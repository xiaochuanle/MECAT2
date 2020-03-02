#ifndef __SPLIT_READS_H
#define __SPLIT_READS_H

#include "../common/range_list.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int qid;
	int sid;
	int qbgn, qend;
	int tbgn, tend;
} AdjustOverlap;

typedef kvec_t(AdjustOverlap) vec_adjovlp;

typedef struct {
	int id;
	int type;
	int bgn;
	int end;
} BadRegion;

typedef kvec_t(BadRegion) vec_bad_region;

#define BadType_nothing	0
#define BadType_5spur	1
#define BadType_3spur	2
#define BadType_chimera 3
#define BadType_subread	4

void
split_reads_for_one_partition(const char* pm4_dir, 
							  const int pid, 
							  const int min_size,
							  ClippedRange* clipped_ranges,
							  ClippedRange* split_range,
							  const int num_threads);

#ifdef __cplusplus
}
#endif

#endif // __SPLIT_READS_H