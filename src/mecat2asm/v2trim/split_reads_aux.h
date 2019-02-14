#ifndef SPLIT_READS_AUX_H
#define SPLIT_READS_AUX_H

#include "../klib/kvec.h"
#include "range_list.h"

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
split_reads_for_one_partition(const char* m4_path, 
							  const int pid, 
							  const int min_size,
							  ClippedRange* clipped_ranges,
							  ClippedRange* split_range,
							  const int num_threads);

#endif // SPLIT_READS_AUX_H
