#include "hbn_task_struct.h"

hbn_task_struct*
hbn_task_struct_new(const HbnProgramOptions* opts)
{
    hbn_task_struct* ht_struct = (hbn_task_struct*)calloc(1, sizeof(hbn_task_struct));

    ht_struct->qi_vs_sj_out = NULL;
    ht_struct->out = NULL;
    pthread_mutex_init(&ht_struct->out_lock , NULL);

    const int query_and_subject_are_the_same = strcmp(opts->query, opts->subject) == 0;

    ht_struct->query_db_title = INIT_QUERY_DB_TITLE;
    ht_struct->query_vol_index = -1;
    ht_struct->query_vol = NULL;

    ht_struct->subject_db_title = (query_and_subject_are_the_same) ? INIT_QUERY_DB_TITLE : INIT_SUBJECT_DB_TITLE;
    ht_struct->subject_vol_index = -1;
    ht_struct->subject_vol = NULL;

    ht_struct->opts = opts;
    ht_struct->query_and_subject_are_the_same = query_and_subject_are_the_same;
    ht_struct->lktbl = NULL;
    hbn_assert(ht_struct->opts->num_threads > 0);
    ht_struct->word_data_array = (WordFindData**)calloc(ht_struct->opts->num_threads, sizeof(WordFindData*));
    ht_struct->hit_extn_data_array = (HbnSubseqHitExtnData**)calloc(ht_struct->opts->num_threads, sizeof(HbnSubseqHitExtnData*));
    ht_struct->results_array = (HbnHSPResults**)calloc(ht_struct->opts->num_threads, sizeof(HbnHSPResults*));
    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->results_array[i] = HbnHSPResultsNew(HBN_QUERY_CHUNK_SIZE);
        ht_struct->hit_extn_data_array[i] = HbnSubseqHitExtnDataNew(opts);
    }

    hbn_fopen(ht_struct->out, ht_struct->opts->output, "w");

    return ht_struct;
}

hbn_task_struct*
hbn_task_struct_free(hbn_task_struct* ht_struct)
{
    if (ht_struct->out) hbn_fclose(ht_struct->out);
    ht_struct->out = NULL;

    hbn_task_struct_destroy_query_vol_context(ht_struct);
    hbn_task_struct_destroy_subject_vol_context(ht_struct);

    if (!ht_struct->opts->keep_db) {
        char cmd[HBN_MAX_PATH_LEN];
        sprintf(cmd, "rm -r %s", ht_struct->opts->db_dir);
        hbn_system(cmd);
    }

    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->hit_extn_data_array[i] = HbnSubseqHitExtnDataFree(ht_struct->hit_extn_data_array[i]);
        ht_struct->results_array[i] = HbnHSPResultsFree(ht_struct->results_array[i]);
    }
    free(ht_struct->hit_extn_data_array);
    free(ht_struct->results_array);
    free(ht_struct->word_data_array);
    free(ht_struct);
    return NULL;
}

void
hbn_task_struct_destroy_query_vol_context(hbn_task_struct* ht_struct)
{
    if (ht_struct->query_vol) {
        hbn_assert(ht_struct->query_vol_index >= 0);
        hbn_assert(ht_struct->qi_vs_sj_out);
        hbn_fclose(ht_struct->qi_vs_sj_out);
        CSeqDBFree(ht_struct->query_vol);
    }
    ht_struct->query_vol = NULL;
    ht_struct->query_vol_index = -1;
    ht_struct->qi_vs_sj_out = NULL;
}

void
hbn_task_struct_build_query_vol_context(hbn_task_struct* ht_struct, int query_vol_index)
{
    hbn_task_struct_destroy_query_vol_context(ht_struct);
    ht_struct->query_vol = seqdb_load(ht_struct->opts->db_dir, ht_struct->query_db_title, query_vol_index);
    ht_struct->query_vol_index = query_vol_index;
    hbn_assert(ht_struct->subject_vol);
    hbn_assert(ht_struct->subject_vol_index >= 0);
    ht_struct->qi_vs_sj_out = open_qi_vs_sj_results_file(ht_struct->opts->db_dir,
                                kBackupAlignResultsDir,
                                query_vol_index,
                                ht_struct->subject_vol_index);
}

void
hbn_task_struct_destroy_subject_vol_context(hbn_task_struct* ht_struct)
{
    if (ht_struct->subject_vol) {
        hbn_assert(ht_struct->subject_vol_index >= 0);
        hbn_assert(ht_struct->lktbl);
        hbn_assert(ht_struct->word_data_array);
        CSeqDBFree(ht_struct->subject_vol);
        destroy_lookup_table(ht_struct->lktbl);
        for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
            ht_struct->word_data_array[i] = WordFindDataFree(ht_struct->word_data_array[i]);
        }
    }
    ht_struct->subject_vol = NULL;
    ht_struct->subject_vol_index = -1;
    ht_struct->lktbl = NULL;
}

void
hbn_task_struct_build_subject_vol_context(hbn_task_struct* ht_struct, int subject_vol_index)
{
    hbn_task_struct_destroy_subject_vol_context(ht_struct);
    ht_struct->subject_vol_index = subject_vol_index;
    ht_struct->subject_vol = seqdb_load_unpacked(ht_struct->opts->db_dir, ht_struct->subject_db_title, subject_vol_index);
    ht_struct->lktbl = build_lookup_table(ht_struct->subject_vol,
                            ht_struct->opts->kmer_size,
                            ht_struct->opts->kmer_window_size,
                            ht_struct->opts->max_kmer_occ,
                            ht_struct->opts->num_threads);
    set_kmer_block_size_info(ht_struct->opts->block_size);
    for (int i = 0; i < ht_struct->opts->num_threads; ++i) {
        ht_struct->word_data_array[i] = WordFindDataNew(ht_struct->subject_vol, 
                                    ht_struct->lktbl, 
                                    ht_struct->opts->kmer_size, 
                                    1, 
                                    ht_struct->opts->min_ddfs, 
                                    ht_struct->query_and_subject_are_the_same);
    }
}