#ifndef __HBN_TASK_STRUCT_H
#define __HBN_TASK_STRUCT_H

#include "cns_aux.h"

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const HbnProgramOptions* opts;
    int pid;
    HbnConsensusInitHit* cns_hit_array;
    size_t cns_hit_count;
    size_t cns_hit_idx;
    RawReadsReader* raw_reads;
    RawReadCnsInfo* cns_info_array;
    int cns_info_count;
    int cns_info_idx;
    pthread_mutex_t cns_info_lock;
    kstring_t output_buf;
    pthread_mutex_t out_lock;
    FILE* out;
    CnsThreadData** thread_data_array;
} hbn_task_struct;

hbn_task_struct*
hbn_task_struct_new(const HbnProgramOptions* opts);

hbn_task_struct*
hbn_task_struct_free(hbn_task_struct* ht_struct);

void
hbn_task_struct_load_partition_info(hbn_task_struct* ht_struct, const int pid);

BOOL
hbn_task_struct_load_batch_info(hbn_task_struct* ht_struct);

void
hbn_task_struct_dump_results(hbn_task_struct* ht_struct);

#ifdef __cplusplus
}
#endif

#endif // __HBN_TASK_STRUCT_H