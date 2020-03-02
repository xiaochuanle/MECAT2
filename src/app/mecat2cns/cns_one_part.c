#include "cns_one_part.h"

#include "cns_one_read.h"

#include <pthread.h>

static void*
cns_thread_worker(void* params)
{
    CnsThreadData* data = (CnsThreadData*)(params);
    while (1) {
        int cns_info_idx = -1;
        pthread_mutex_lock(data->cns_info_idx_lock);
        cns_info_idx = *data->cns_info_idx;
        ++(*data->cns_info_idx);
        pthread_mutex_unlock(data->cns_info_idx_lock);
        if (cns_info_idx >= data->cns_info_count) break;
        consensus_one_read(data, cns_info_idx);
        //exit(0);
    }
    return NULL;
}

void
cns_one_part(hbn_task_struct* ht_struct, const int pid)
{
    hbn_task_struct_load_partition_info(ht_struct, pid);
    const int num_threads = ht_struct->opts->num_threads;
    pthread_t jobid_array[num_threads];
    char job_name[1000];
    char buf1[64];
    char buf2[64];
    while (hbn_task_struct_load_batch_info(ht_struct)) {
        int min_subject_id = ht_struct->cns_info_array[0].oid;
        int max_subject_id = ht_struct->cns_info_array[ht_struct->cns_info_count-1].oid + 1;
        u64_to_fixed_width_string_r(min_subject_id, buf1, HBN_DIGIT_WIDTH);
        u64_to_fixed_width_string_r(max_subject_id, buf2, HBN_DIGIT_WIDTH);
        sprintf(job_name, "correcting subject %s --- %s", buf1, buf2);
        hbn_timing_begin(job_name);
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobid_array + i, NULL, cns_thread_worker, ht_struct->thread_data_array[i]);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobid_array[i], NULL);
        }
        hbn_task_struct_dump_results(ht_struct);
        hbn_timing_end(job_name);
        //break;
    }
}