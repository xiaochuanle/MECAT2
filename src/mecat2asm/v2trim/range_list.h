#ifndef RANGE_LIST_H
#define RANGE_LIST_H

#include "../klib/ksort.h"
#include "../klib/kvec.h"

typedef struct {
	int left, right, size;
} ClippedRange;

typedef kvec_t(ClippedRange) vec_ClippedRange;

typedef struct {
	int lo;
	int hi;
} ClearRange;

#define is_deleted_range(clr) ((clr).lo == -1)
typedef  kvec_t(ClearRange) vec_clr_range;

typedef struct {
	int lo;
	int hi;
	int ct; // number of source ranges
	int va; // value;
} CovInterval;

#define cov_intv_lt(lhs, rhs) (((lhs).lo < (rhs).lo) || ((lhs).lo == (rhs).lo && (lhs).hi < (rhs).hi))
#define cov_intv_is_invalid(covi) ((covi).lo == 0 && (covi).hi == 0)
#define invalid_cov_intv(covi) ((covi).lo = 0, (covi).hi = 0)
typedef kvec_t(CovInterval) vec_cov_intv;

typedef struct {
	int				is_sorted;
	int				is_merged;
	vec_cov_intv	list;
} interval_list;

#define intv_list_is_merged(ilist) ((ilist).is_merged)
#define intv_list_is_sorted(ilist) ((ilist).is_sorted)
#define intv_list_lo(ilist, i)		(kv_A((ilist).list, i).lo)
#define intv_list_hi(ilist, i)		(kv_A((ilist).list, i).hi)
#define intv_list_count(ilist, i)	(kv_A((ilist).list, i).ct)
#define intv_list_depth(ilist, i)	(kv_A((ilist).list, i).ct)
#define intv_list_value(ilist, i)	(kv_A((ilist).list, i).va)

#ifdef __cplusplus
extern "C" {
#endif

void init_intv_list(interval_list* list);
void destroy_intv_list(interval_list* list);
void copy_intv_list(interval_list* dst, interval_list* src);
void add_intv_list(interval_list* list, int position, int length, int val);
void sort_intv_list(interval_list* list);
void merge_intv_list(interval_list* list, int min_ovlp);
void depth_from_intv_list(interval_list* dst, interval_list* src);
void invert_intv_list(interval_list* list, int invlo, int invhi);

#ifdef __cplusplus
}
#endif

#endif // RANGE_LIST_H
