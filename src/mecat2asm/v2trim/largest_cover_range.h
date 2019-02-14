#ifndef LARGEST_COVER_RANGE_H
#define LARGEST_COVER_RANGE_H

#include "m4_record.h"
#include "range_list.h"

void
get_largest_cover_range_for_one_partition(const char* m4_path, 
		const int pid, 
		const double min_ident_perc,
		const int min_ovlp_size,
		const int min_cov,
		ClippedRange* clipped_ranges,
		const int num_threads);

#endif // LARGEST_COVER_RANGE_H
