#include "range_list.h"

#include <assert.h>

#include "../../../corelib/hbn_aux.h"

KSORT_INIT(cov_intv_lt, CovInterval, cov_intv_lt)

void init_intv_list(interval_list* list)
{
	list->is_merged = 0;
	list->is_sorted = 0;
	kv_init(list->list);
}

void destroy_intv_list(interval_list* list)
{
	kv_destroy(list->list);
}

void clear_intv_list(interval_list* list)
{
	list->is_merged = 0;
	list->is_sorted = 0;
	kv_clear(list->list);
}

void copy_intv_list(interval_list* dst, interval_list* src)
{
	dst->is_merged = src->is_merged;
	dst->is_sorted = src->is_sorted;
	kv_copy(CovInterval, dst->list, src->list);
}

void add_intv_list(interval_list* list, int position, int length, int val)
{
	CovInterval covi;
	covi.lo = position;
	covi.hi = position + length;
	covi.ct = 1;
	covi.va = val;
	kv_push(CovInterval, list->list, covi);
	list->is_merged = 0;
	list->is_sorted = 0;
}

void sort_intv_list(interval_list* list)
{
	if (!list->is_sorted) ks_introsort_cov_intv_lt(kv_size(list->list), kv_data(list->list));
	list->is_sorted = 1;
}

void merge_intv_list(interval_list* list, int min_ovlp)
{
	if (list->is_merged) return;
	int curr = 0, next = 1;
	sort_intv_list(list);
	int nrange = kv_size(list->list);
	CovInterval* intv_list = kv_data(list->list);

	while (next < nrange) {
		if (cov_intv_is_invalid(intv_list[curr])) {
			intv_list[curr] = intv_list[next];
			invalid_cov_intv(intv_list[next]);
			++next;
		} else {
			int intersect = 0;
			/// contained
			if (intv_list[curr].lo <= intv_list[next].lo && intv_list[next].hi <= intv_list[curr].hi) intersect = 1;
			/// overlap
			if (intv_list[curr].hi - min_ovlp >= intv_list[next].lo) intersect = 1;
			if (intersect) {
				if (intv_list[curr].hi < intv_list[next].hi) intv_list[curr].hi = intv_list[next].hi;
				intv_list[curr].ct += intv_list[next].ct;
				intv_list[curr].va += intv_list[next].va;
				invalid_cov_intv(intv_list[next]);
				++next;
			} else {
				++curr;
				if (curr != next) intv_list[curr] = intv_list[next];
				++next;
			}
		}
	}
	if (curr + 1 < nrange) kv_resize(CovInterval, list->list, curr + 1);
	list->is_merged = 1;
}

typedef struct {
	int pos; // position of the change in depth
	int change; // the value associated with this object; added or subtracted from va
	int open; // is true, the start of a new interval
} IntvDepthRegion;

#define intv_depth_lt(a, b) (((a).pos < (b).pos) || ((a).pos == (b).pos && (a).open > (b).open))
KSORT_INIT(intv_depth_lt, IntvDepthRegion, intv_depth_lt);

static void compute_depth_for_intv_list(interval_list* list, IntvDepthRegion* id, size_t id_len)
{
	clear_intv_list(list);
	if (id_len == 0) return;

	ks_introsort_intv_depth_lt(id_len, id);
	kv_reserve(CovInterval, list->list, id_len);
	CovInterval* intv_list = kv_data(list->list);
	size_t list_len = 0;
	/// init the first interval
	assert(id[0].open == 1);
	intv_list[list_len].lo = id[0].pos;
	intv_list[list_len].hi = id[0].pos;
	intv_list[list_len].ct = 1;
	intv_list[list_len].va = id[0].change;

	int nct;
	int nva;

	for (size_t i = 1; i < id_len; ++i) {
		intv_list[list_len].hi = id[i].pos;
		if (id[i].open) {
			nct = intv_list[list_len].ct + 1;
			nva = intv_list[list_len].va + id[i].change;
		} else {
			nct = intv_list[list_len].ct - 1;
			nva = intv_list[list_len].va - id[i].change;
		}

		int r = (id[i-1].pos != id[i].pos) || (intv_list[list_len].va != nva);
		if (r) r = intv_list[list_len].lo != intv_list[list_len].hi;
		if (r) {
			++list_len;
			intv_list[list_len].lo = id[i].pos;
			intv_list[list_len].ct = intv_list[list_len - 1].ct;
			intv_list[list_len].va = intv_list[list_len - 1].va;
		}

		intv_list[list_len].hi = id[i].pos;
		intv_list[list_len].ct = nct;
		intv_list[list_len].va = nva;

		r = (list_len > 1) &&
			(intv_list[list_len - 1].hi == intv_list[list_len].lo) &&
			(intv_list[list_len - 1].ct == intv_list[list_len].ct) &&
			(intv_list[list_len - 1].va == intv_list[list_len].va);
		if (r) {
			intv_list[list_len - 1].hi = intv_list[list_len].hi;
			--list_len;
		}
	}

	assert(list_len > 0);
	kv_resize(CovInterval, list->list, list_len);
}

void depth_from_intv_list(interval_list* dst, interval_list* src)
{
	size_t id_len = 2 * kv_size(src->list);
	IntvDepthRegion* id = (IntvDepthRegion*)malloc(sizeof(IntvDepthRegion) * id_len);
	for (size_t i = 0; i < kv_size(src->list); ++i) {
		id[2 * i    ].pos		= intv_list_lo(*src, i);
		id[2 * i    ].change	= intv_list_value(*src, i);
		id[2 * i    ].open		= 1;
		id[2 * i + 1].pos		= intv_list_hi(*src, i);
		id[2 * i + 1].change	= intv_list_value(*src, i);
		id[2 * i + 1].open		= 0;
	}

	compute_depth_for_intv_list(dst, id, id_len);
	free(id);
}

void invert_intv_list(interval_list* list, int invlo, int invhi)
{
	merge_intv_list(list, 0);

	int list_len = kv_size(list->list);
	int inv_len = 0;
	int inv_max = list_len + 2;
	CovInterval* inv = (CovInterval*)malloc(sizeof(CovInterval) * inv_max);

	if (list_len == 0) {
		inv[inv_len].lo = invlo;
		inv[inv_len].hi = invhi;
		inv[inv_len].ct = 1;
		inv[inv_len].va = 0;
		++inv_len;
	} else {
		if (invlo < intv_list_lo(*list, 0)) {
			inv[inv_len].lo = invlo;
			inv[inv_len].hi = intv_list_lo(*list, 0);
			inv[inv_len].ct = 1;
			inv[inv_len].va = 0;
			++inv_len;
		}

		for (size_t i = 1; i < kv_size(list->list); ++i) {
			if (intv_list_hi(*list, i-1) < intv_list_lo(*list, i)) {
				inv[inv_len].lo = intv_list_hi(*list, i-1);
				inv[inv_len].hi = intv_list_lo(*list, i);
				inv[inv_len].ct = 1;
				inv[inv_len].va = 0;
				++inv_len;
			}
		}

		if (intv_list_hi(*list, list_len - 1) < invhi) {
			inv[inv_len].lo = intv_list_hi(*list, list_len - 1);
			inv[inv_len].hi = invhi;
			inv[inv_len].ct = 1;
			inv[inv_len].va = 0;
			++inv_len;
		}
	}

	assert(inv_len <= inv_max);
	kv_clear(list->list);
	for (int i = 0; i < inv_len; ++i) {
		kv_push(CovInterval, list->list, inv[i]);
	}
	free(inv);
}