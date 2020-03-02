#ifndef __LARGEST_COVER_RANGE_H
#define __LARGEST_COVER_RANGE_H

#include "../common/range_list.h"

#ifdef __cplusplus
extern "C" {
#endif

void
get_largest_cover_range_for_one_partition(const char* pm4_dir, 
		const int pid, 
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		ClippedRange* clipped_ranges,
		const int num_threads);

#ifdef __cplusplus
}
#endif

#endif // __LARGEST_COVER_RANGE_H