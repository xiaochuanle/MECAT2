#include "split_reads.h"

#include "../../../corelib/khash.h"
#include "../../../corelib/ksort.h"
#include "../common/m4_aux.h"

#include <assert.h>
#include <math.h>

#define SUBREAD_LOOP_MAX_SIZE 500
#define SUBREAD_LOOP_EXT_SIZE 2000

KHASH_MAP_INIT_INT(32, int)

static int
adjust_offsets(ClippedRange* clr,
	M4Record* m4,
	int* sbgn,
	int* send,
	int* qbgn,
	int* qend)
{
	int qid = m4->qid;
	if (clr[qid].size == 0) return 0;
	int tid = m4->sid;
	assert(m4->sdir == FWD);
	hbn_assert(clr[tid].size == m4->ssize);
	int tovlbgn = m4->soff;
	int tovlend = m4->send;
	assert(clr[qid].size == m4->qsize);
	assert(clr[tid].size == m4->ssize);
	int tclrbgn = clr[tid].left;
	int tclrend = clr[tid].right;

	int qovlbgn;
	int qovlend;
	int qclrbgn;
	int qclrend;
	if (m4->qdir == FWD) {
		qovlbgn = m4->qoff;
		qovlend = m4->qend;
		qclrbgn = clr[qid].left;
		qclrend = clr[qid].right;
	} else {
		assert(m4->qdir == REV);
		qovlbgn = m4->qsize - m4->qend;
		qovlend = m4->qsize - m4->qoff;
		qclrbgn = m4->qsize - clr[qid].right;
		qclrend = m4->qsize - clr[qid].left;
	}

	// overlap doesn't intersect the clear range
	if (qclrend <= qovlbgn ||
		qovlend <= qclrbgn ||
		tclrend <= tovlbgn ||
		tovlend <= tclrbgn) {
		return 0;
	}

	int qlen = qovlend - qovlbgn;
	int tlen = tovlend - tovlbgn;
	double qfbgn = 1.0 * ((qclrbgn < qovlbgn) ? 0 : (qclrbgn - qovlbgn)) / qlen;
	double tfbgn = 1.0 * ((tclrbgn < tovlbgn) ? 0 : (tclrbgn - tovlbgn)) / tlen;
	double qfend = 1.0 * ((qclrend > qovlend) ? 0 : (qovlend - qclrend)) / qlen;
	double tfend = 1.0 * ((tclrend > tovlend) ? 0 : (tovlend - tclrend)) / tlen;
	double maxbgn = (qfbgn > tfbgn) ? qfbgn : tfbgn;
	double maxend = (qfend > tfend) ? qfend : tfend;
	assert(maxbgn < 1.0);
	assert(maxend < 1.0);

	int qadjbgn = (int)round(maxbgn * qlen);
	int tadjbgn = (int)round(maxbgn * tlen);
	int qadjend = (int)round(maxend * qlen);
	int tadjend = (int)round(maxend * tlen);

	qovlbgn += qadjbgn;
	tovlbgn += tadjbgn;
	qovlend -= qadjend;
	tovlend -= tadjend;
	assert(qclrbgn <= qovlbgn);
	assert(tclrbgn <= tovlbgn);
	assert(qovlend <= qclrend);
	assert(tovlend <= tclrend);

	// reset hangs that are close to zero
	if (qovlbgn - qclrbgn < 15) {
		int limit = tovlbgn - tclrbgn;
		int adjust = qovlbgn - qclrbgn;
		if (adjust > limit) adjust = limit;
		assert(tovlbgn >= adjust);
		assert(qovlbgn >= adjust);
		tovlbgn -= adjust;
		qovlbgn -= adjust;
		assert(tclrbgn <= tovlbgn);
		assert(qclrbgn <= qovlbgn);
	}

	if (qclrend - qovlend < 15) {
		int limit = tclrend - tovlend;
		int adjust = qclrend - qovlend;
		if (adjust > limit) adjust = limit;
		tovlend += adjust;
		qovlend += adjust;
		assert(tovlend <= tclrend);
		assert(qovlend <= qclrend);
	}

	*sbgn = tovlbgn;
	*send = tovlend;
	*qbgn = (m4->qdir == FWD) ? qovlbgn : (m4->qsize - qovlend);
	*qend = (m4->qdir == FWD) ? qovlend : (m4->qsize - qovlbgn);
	return 1;
}

void add_and_filter_overlaps(M4Record* m4v,
	const int nm4,
	ClippedRange* clr,
	vec_adjovlp* adjovlp)
{
	AdjustOverlap aov;
	for (int i = 0; i < nm4; ++i) {
		aov.qid = m4v[i].qid;
		if (!adjust_offsets(clr, m4v + i, &aov.tbgn, &aov.tend, &aov.qbgn, &aov.qend)) continue;
		kv_push(AdjustOverlap, *adjovlp, aov);
	}
}

static int interval_overlap(int b1, int e1, int b2, int e2)
{
	int minmax = (e1 < e2) ? e1 : e2;
	int maxmin = (b1 > b2) ? b1 : b2;
	return (minmax > maxmin) ? (minmax - maxmin) : 0;
}

void detect_subread(const int tid, AdjustOverlap* adjovlp, const int nov, vec_bad_region* blist)
{
	khash_t(32)* next_idx = kh_init(32);
	khash_t(32)* num_ovlps = kh_init(32);
	int large_palindrome = 0;
	int r;
	interval_list bad;
	interval_list bad_all;
	init_intv_list(&bad);
	init_intv_list(&bad_all);

	for (int i = 0; i < nov; ++i) {
		int qid = adjovlp[i].qid;
		khiter_t k = kh_put(32, next_idx, qid, &r);
		kh_value(next_idx, k) = i;
		k = kh_put(32, num_ovlps, qid, &r);
		if (r == 1) kh_value(num_ovlps, k) = 1;
		else ++kh_value(num_ovlps, k);
	}

	for (int i = 0; i < nov; ++i) {
		int qid = adjovlp[i].qid;
		khiter_t k = kh_get(32, num_ovlps, qid);
		assert(kh_exist(num_ovlps, k));
		int n = kh_value(num_ovlps, k);
		if (n != 2) continue;
		k = kh_get(32, next_idx, qid);
		assert(kh_exist(next_idx, k));
		int j = kh_value(next_idx, k);
		assert(j < nov);
		if (i == j) continue;
		assert(adjovlp[i].qid == adjovlp[j].qid);

		int tovlp = interval_overlap(adjovlp[i].tbgn, adjovlp[i].tend, adjovlp[j].tbgn, adjovlp[j].tend);
		int qovlp = interval_overlap(adjovlp[i].qbgn, adjovlp[i].qend, adjovlp[j].qbgn, adjovlp[j].qend);
		if (tovlp == 0 && qovlp == 0) continue;
		if (tovlp > 1000 && qovlp > 1000) large_palindrome = 1;
		if (tovlp > 250 || qovlp < 250) continue;

		int badbgn = (adjovlp[i].tbgn < adjovlp[j].tbgn) ? adjovlp[i].tend : adjovlp[j].tend;
		int badend = (adjovlp[i].tbgn < adjovlp[j].tbgn) ? adjovlp[j].tbgn : adjovlp[i].tbgn;
		if (badbgn > badend) { int a = badbgn; badbgn = badend; badend = a; }
		assert(badbgn <= badend);

		if (badend - badbgn <= SUBREAD_LOOP_MAX_SIZE) add_intv_list(&bad, badbgn, badend - badbgn, 0);
		if (badend - badbgn <= SUBREAD_LOOP_EXT_SIZE) add_intv_list(&bad_all, badbgn, badend - badbgn, 0);
	}
	merge_intv_list(&bad, 0);
	merge_intv_list(&bad_all, 0);

	for (size_t bb = 0; bb < kv_size(bad.list); ++bb) {
		int num_span = 0;
		int all_hits = 0;
		int bad_lo = intv_list_lo(bad, bb);
		int bad_hi = intv_list_hi(bad, bb);

		for (size_t aa = 0; aa < kv_size(bad_all.list); ++aa) {
			int bad_all_lo = intv_list_lo(bad_all, aa);
			int bad_all_hi = intv_list_hi(bad_all, aa);
			if (bad_all_lo <= bad_lo && bad_hi <= bad_all_hi) all_hits += intv_list_count(bad_all, aa);
		}
		assert(all_hits > 0);

		for (int i = 0; i < nov; ++i) {
			if (adjovlp[i].tbgn + 100 < bad_lo && bad_hi + 100 < adjovlp[i].tend) num_span += 1;
		}

		if (num_span > 9) continue;
		if (intv_list_count(bad, bb) + all_hits / 4 + large_palindrome < 3) continue;

		BadRegion region;
		region.id = tid;
		region.type = BadType_subread;
		region.bgn = bad_lo;
		region.end = bad_hi;
		kv_push(BadRegion, *blist, region);
	}

	kh_destroy(32, next_idx);
	kh_destroy(32, num_ovlps);
	destroy_intv_list(&bad);
	destroy_intv_list(&bad_all);
}

int trim_bad_interval(vec_bad_region* blist, int min_size, int* clrbgn, int* clrend)
{
	if (kv_size(*blist) == 0) return 1;
	interval_list good_regions;
	init_intv_list(&good_regions);

	for (size_t bb = 0; bb < kv_size(*blist); ++bb) {
		BadRegion region = kv_A(*blist, bb);
		add_intv_list(&good_regions, region.bgn, region.end - region.bgn, 0);
	}

	int cbgn = *clrbgn;
	int cend = *clrend;
	invert_intv_list(&good_regions, cbgn, cend);

	cbgn = 0;
	cend = 0;
	for (size_t rr = 0; rr < kv_size(good_regions.list); ++rr) {
		int lo = intv_list_lo(good_regions, rr);
		int hi = intv_list_hi(good_regions, rr);
		if (cend - cbgn < hi - lo) {
			cbgn = lo;
			cend = hi;
		}
	}

	*clrbgn = cbgn;
	*clrend = cend;
	destroy_intv_list(&good_regions);
	return cend - cbgn >= min_size;
}

typedef struct {
	M4Record* m4v;
	int nm4;
	int* idx_range;
	int nrange;
	int next_range_id;
	pthread_mutex_t range_get_lock;
	ClippedRange* clear_ranges;
	pthread_mutex_t range_set_lock;
	int min_size;
	ClippedRange* split_ranges;
} SplitReadData;

SplitReadData*
new_SplitReadData(M4Record* m4v,
			int nm4,
			int* idx_range,
			int nrange,
			ClippedRange* clear_ranges,
		   	int min_size,
		    ClippedRange* split_ranges)
{
	SplitReadData* data = (SplitReadData*)malloc(sizeof(SplitReadData));
	data->m4v = m4v;
	data->nm4 = nm4;
	data->idx_range = idx_range;
	data->nrange = nrange;
	data->next_range_id = 0;
	pthread_mutex_init(&data->range_get_lock, NULL);
	data->clear_ranges = clear_ranges;
	pthread_mutex_init(&data->range_set_lock, NULL);
	data->min_size = min_size;
	data->split_ranges = split_ranges;
	return data;
}

SplitReadData*
free_SplitReadData(SplitReadData* data)
{
	free(data);
	return NULL;
}

static void
preprocess_m4v(M4Record* m4v, int* nm4)
{
	const int N = 300;
	int n = *nm4;
	if (n <= N) return;
	ks_introsort_m4_ident_gt(n, m4v);
	*nm4 = N;
}

void*
split_read_worker(void* arg)
{
	SplitReadData* data = (SplitReadData*)(arg);
	M4Record* m4v = NULL;
	int nm4;
	int read_id, left, right, size;
	kv_dinit(vec_adjovlp, adjovlp);
	kv_dinit(vec_bad_region, blist);
	while (1) {
		m4v = get_next_range(data->m4v, 
							 &data->range_get_lock, 
							 data->idx_range, 
							 &data->next_range_id, 
							 data->nrange, 
							 &nm4);
		if (!m4v) break;
        if (data->clear_ranges[m4v[0].sid].size == 0) continue;
		hbn_assert(data->clear_ranges[m4v[0].sid].size == m4v[0].ssize);
		preprocess_m4v(m4v, &nm4);
		kv_clear(adjovlp);
		add_and_filter_overlaps(m4v, nm4, data->clear_ranges, &adjovlp);
		kv_clear(blist);
		detect_subread(m4v[0].sid, kv_data(adjovlp), kv_size(adjovlp), &blist);
		left = data->clear_ranges[m4v[0].sid].left;
		right = data->clear_ranges[m4v[0].sid].right;
		int r = trim_bad_interval(&blist, data->min_size, &left, &right);
		if (r) {
			pthread_mutex_lock(&data->range_set_lock);
			read_id = m4v[0].sid;
			size = m4v[0].ssize;
			hbn_assert(data->split_ranges[read_id].size == size);
			data->split_ranges[read_id].left = left;
			data->split_ranges[read_id].right = right;
			data->split_ranges[read_id].size = size;
			pthread_mutex_unlock(&data->range_set_lock);
		}
	}
    kv_destroy(adjovlp);
    kv_destroy(blist);
	return NULL;
}

void
split_reads_for_one_partition(const char* pm4_dir, 
		const int pid, 
		const int min_size,
		ClippedRange* clipped_ranges,
		ClippedRange* split_ranges,
		const int num_threads)
{
	kv_dinit(vec_int, idx_range);
    size_t m4_count;
	M4Record* m4_array = load_partition_m4(pm4_dir, pid, &m4_count, &idx_range);
    if (m4_count == 0) return;
	SplitReadData* sr_data = new_SplitReadData(m4_array,
							 m4_count,
							 kv_data(idx_range),
							 kv_size(idx_range) - 1,
							 clipped_ranges,
							 min_size,
							 split_ranges);
	
	pthread_t jobs[num_threads];
	for (int i = 0; i < num_threads; ++i) {
		pthread_create(jobs + i, NULL, split_read_worker, (void*)(sr_data));
	}
	for (int i = 0; i < num_threads; ++i) {
		pthread_join(jobs[i], NULL);
	}
	free_SplitReadData(sr_data);
	free(m4_array);
	kv_destroy(idx_range);
}