#include "range_list.h"

#include "pm4_aux.h"
#include "range_list.h"
#include "../common/ontcns_aux.h"
#include "../common/oc_assert.h"
#include "../klib/ksort.h"

int largest_cover_range(M4Record* m4v,
	const int nm4,
	int* fbgn,
	int* fend,
	const double min_ident_perc,
	const int min_ovlp_size,
	const int min_cov)
{
	interval_list IL, ID;
	init_intv_list(&IL);
	init_intv_list(&ID);

	for (int i = 0; i < nm4; ++i) {
		int tbgn = m4v[i].soff;
		int tend = m4v[i].send;
		add_intv_list(&IL, tbgn, tend - tbgn, 0);
	}

	if (min_cov > 0) {
		interval_list DE;
		init_intv_list(&DE);
		depth_from_intv_list(&DE, &IL);
		size_t it = 0;
		int ib = 0, ie = 0;

		while (it < kv_size(DE.list)) {
			if (intv_list_depth(DE, it) < min_cov) {
				if (ie > ib) add_intv_list(&ID, ib, ie - ib, 0);
				ib = 0;
				ie = 0;
			} else if(ib == 0 && ie == 0) {
				ib = intv_list_lo(DE, it);
				ie = intv_list_hi(DE, it);
			} else if (ie == intv_list_lo(DE, it)) {
				ie = intv_list_hi(DE, it);
			} else {
				if (ie > ib) add_intv_list(&ID, ib, ie - ib, 0);
				ib = intv_list_lo(DE, it);
				ie = intv_list_hi(DE, it);
			}
			++it;
		}
		if (ie > ib) add_intv_list(&ID, ib, ie - ib, 0);
		destroy_intv_list(&DE);
	}

	merge_intv_list(&IL, min_ovlp_size);

	if (min_cov > 0) {
		interval_list FI;
		init_intv_list(&FI);
		size_t li = 0, di = 0;
		while (li < kv_size(IL.list) && di < kv_size(ID.list)) {
			int ll = intv_list_lo(IL, li);
			int lh = intv_list_hi(IL, li);
			int dl = intv_list_lo(ID, di);
			int dh = intv_list_hi(ID, di);
			int nl = 0;
			int nh = 0;

			if (ll <= dl && dl < lh) {
				nl = dl;
				nh = (lh < dh) ? lh : dh;
			}

			if (dl <= ll && ll < dh) {
				nl = ll;
				nh = (lh < dh) ? lh : dh;
			}

			if (nl < nh) add_intv_list(&FI, nl, nh - nl, 0);
			if (lh <= dh) ++li;
			if (dh <= lh) ++di;
		}
		copy_intv_list(&IL, &FI);
		destroy_intv_list(&FI);
	}

	int max_l = 0, max_r = 0;
	for (size_t i = 0; i < kv_size(IL.list); ++i) {
		int l = kv_A(IL.list, i).lo;
		int h = kv_A(IL.list, i).hi;
		if (h - l > max_r - max_l) {
			max_l = l;
			max_r = h;
		}
	}

	*fbgn = max_l;
	*fend = max_r;
	destroy_intv_list(&IL);
	destroy_intv_list(&ID);

	return max_r > 0;
}

typedef struct {
	M4Record* m4v;
	int nm4;
	int* idx_range;
	int nrange;
	int next_range_id;
	pthread_mutex_t range_get_lock;
	ClippedRange* clipped_ranges;
	pthread_mutex_t range_set_lock;
	double min_ident_perc;
	int min_ovlp_size;
	int min_cov;
} LcrData;

LcrData*
new_LcrData(M4Record* m4v,
			int nm4,
			int* idx_range,
			int nrange,
			ClippedRange* clipped_ranges,
		   	double min_ident_perc,
			int min_ovlp_size,
			int min_cov)
{
	LcrData* data = (LcrData*)malloc(sizeof(LcrData));
	data->m4v = m4v;
	data->nm4 = nm4;
	data->idx_range = idx_range;
	data->nrange = nrange;
	data->next_range_id = 0;
	pthread_mutex_init(&data->range_get_lock, NULL);
	data->clipped_ranges = clipped_ranges;
	pthread_mutex_init(&data->range_set_lock, NULL);
	data->min_ident_perc = min_ident_perc;
	data->min_ovlp_size = min_ovlp_size;
	data->min_cov = min_cov;
	return data;
}

LcrData*
free_LcrData(LcrData* data)
{
	free(data);
	return NULL;
}

static void
preprocess_m4v(M4Record* m4v, int* nm4)
{
	int n = *nm4;
	if (n > 300) {
		ks_introsort_m4_ident_gt(n, m4v);
		*nm4 = 300;
	}
}

void*
lcr_worker(void* arg)
{
	LcrData* data = (LcrData*)(arg);
	M4Record* m4v = NULL;
	int nm4;
	int read_id, left, right, size;
	while (1) {
		m4v = get_next_range(data->m4v, 
							 &data->range_get_lock,
							 data->idx_range,
							 &data->next_range_id,
							 data->nrange,
							 &nm4);
		if (!m4v) break;
		preprocess_m4v(m4v, &nm4);
		BOOL r = largest_cover_range(m4v,
									 nm4,
									 &left,
									 &right,
									 data->min_ident_perc,
									 data-> min_ovlp_size,
									 data->min_cov);
		if (r) {
			pthread_mutex_lock(&data->range_set_lock);
			read_id = m4v[0].sid - 1;
			size = m4v[0].ssize;
			data->clipped_ranges[read_id].left = left;
			data->clipped_ranges[read_id].right = right;
			data->clipped_ranges[read_id].size = size;
			pthread_mutex_unlock(&data->range_set_lock);
		}
	}
	return NULL;
}

void
get_largest_cover_range_for_one_partition(const char* m4_path, 
		const int pid, 
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		ClippedRange* clipped_ranges,
		const int num_threads)
{
	new_kvec(vec_m4, m4v);
	new_kvec(vec_int, idx_range);
	load_partition_m4(m4_path, pid, &m4v, &idx_range);
	if (kv_size(m4v) == 0) {
		free_kvec(m4v);
		free_kvec(idx_range);
		return;
	}
	LcrData* lcr_data = new_LcrData(kv_data(m4v),
									kv_size(m4v),
									kv_data(idx_range),
									kv_size(idx_range) - 1,
									clipped_ranges,
									min_ident_perc,
									min_ovlp_size,
									min_cov);
	pthread_t jobs[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		pthread_create(jobs + i, NULL, lcr_worker, (void*)(lcr_data));
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(jobs[i], NULL);
	}
	free_LcrData(lcr_data);
	free_kvec(m4v);
	free_kvec(idx_range);
}
