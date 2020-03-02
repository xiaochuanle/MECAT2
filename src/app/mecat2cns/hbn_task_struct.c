#include "hbn_task_struct.h"

#include "../../corelib/partition_aux.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

hbn_task_struct*
hbn_task_struct_new(const HbnProgramOptions* opts)
{
    hbn_task_struct* ht_struct = (hbn_task_struct*)calloc(1, sizeof(hbn_task_struct));
    ht_struct->opts = opts;
    ht_struct->raw_reads = RawReadsReaderNew(opts->db_dir, opts->db_title, opts->use_batch_mode);
    ht_struct->cns_info_array = (RawReadCnsInfo*)calloc(opts->batch_size, sizeof(RawReadCnsInfo));
    pthread_mutex_init(&ht_struct->cns_info_lock, NULL);
    ks_init(ht_struct->output_buf);
    pthread_mutex_init(&ht_struct->out_lock, NULL);
    ht_struct->thread_data_array = (CnsThreadData**)calloc(opts->num_threads, sizeof(CnsThreadData*));
    for (int i = 0; i < opts->num_threads; ++i) {
        ht_struct->thread_data_array[i] = CnsThreadDataNew(opts,
                                                ht_struct->raw_reads,
                                                ht_struct->cns_info_array,
                                                &ht_struct->cns_info_idx,
                                                &ht_struct->cns_info_lock,
                                                &ht_struct->output_buf,
                                                &ht_struct->out_lock);
    }
    return ht_struct;
}

hbn_task_struct*
hbn_task_struct_free(hbn_task_struct* ht_struct)
{
    RawReadsReaderFree(ht_struct->raw_reads);
    free(ht_struct->cns_info_array);
    ks_destroy(ht_struct->output_buf);
    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->thread_data_array[i] = CnsThreadDataFree(ht_struct->thread_data_array[i]);
    }
    free(ht_struct->thread_data_array);
    if (ht_struct->out) hbn_fclose(ht_struct->out);
    free(ht_struct);
    return NULL;
}

void
hbn_task_struct_load_partition_info(hbn_task_struct* ht_struct, const int pid)
{
    if (ht_struct->cns_hit_array) sfree(ht_struct->cns_hit_array);
    ht_struct->cns_hit_array = load_and_sort_cns_hits(ht_struct->opts->can_dir,
        pid,
        ht_struct->opts->use_batch_mode,
        &ht_struct->cns_hit_count);
    ht_struct->cns_hit_idx = 0;
    ht_struct->pid = pid;
    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->thread_data_array[i]->cns_hit_array = ht_struct->cns_hit_array;
        ht_struct->thread_data_array[i]->cns_hit_count = ht_struct->cns_hit_count;
    }

    if (ht_struct->out) hbn_fclose(ht_struct->out);
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(ht_struct->opts->can_dir,
        DEFAULT_PART_PREFIX,
        ht_struct->pid,
        path);
    strcat(path, ".cns.fasta");
    hbn_fopen(ht_struct->out, path, "w");
}

BOOL
hbn_task_struct_load_batch_info(hbn_task_struct* ht_struct)
{
    BOOL r = set_next_raw_read_batch_info(ht_struct->cns_info_array,
                &ht_struct->cns_info_count,
                ht_struct->cns_hit_array,
                ht_struct->cns_hit_count,
                &ht_struct->cns_hit_idx,
                ht_struct->raw_reads,
                ht_struct->opts->batch_size);

    if (!r) return r;
    ht_struct->cns_info_idx = 0;
    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->thread_data_array[i]->cns_info_count = ht_struct->cns_info_count;
    }
    return TRUE;
}

void
hbn_task_struct_dump_results(hbn_task_struct* ht_struct)
{
    for (int i = 0; i < ht_struct->cns_info_count; ++i) {
        RawReadCnsInfo* cns_info = ht_struct->cns_info_array + i;
        if (cns_info->cns_fasta_size == 0) continue;
        hbn_assert(cns_info->cns_fasta_offset < ks_size(ht_struct->output_buf));
        const char* s = ks_s(ht_struct->output_buf) + cns_info->cns_fasta_offset;
        hbn_fwrite(s, 1, cns_info->cns_fasta_size, ht_struct->out);
    }
    ks_clear(ht_struct->output_buf);
}