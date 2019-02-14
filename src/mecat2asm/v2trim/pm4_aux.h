#ifndef PM4_AUX_H
#define PM4_AUX_H

#include "../klib/kstring.h"
#include "m4_record.h"

void
make_partition_name(const char* prefix, const int pid, kstring_t* name);

void
make_partition_index_name(const char* prefix, kstring_t* name);

int 
load_num_partitions(const char* m4_path);

int
load_num_reads(const char* wrk_dir);

void
dump_num_partitions(const char* m4_path, const int np);

void
load_partition_m4(const char* m4_path, const int pid, vec_m4* m4v, vec_int* idx_range);

M4Record*
get_next_range(M4Record* m4v, pthread_mutex_t* range_get_lock, int* idx_range, int* next_range_id, int nrange, int* nm4);

typedef struct {
	int nw;
	int mnw;
	FILE** file_list;
} AsmM4Writers;

void
pm4_main(const char* wrk_dir, 
		 const char* m4_path, 
		 const int num_dumpped_files, 
		 const int partition_size, 
		 const double min_ident_perc,
		 const int num_threads);

#endif // PCAN_AUX_H
