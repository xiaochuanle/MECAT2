#ifndef __HBN_TASK_STRUCT_H
#define __HBN_TASK_STRUCT_H

#include "cmdline_args.h"
#include "hbn_extend_subseq_hit.h"
#include "hbn_job_control.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/build_db.h"
#include "../../algo/hbn_lookup_table.h"
#include "../../algo/word_finder.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    FILE*               qi_vs_sj_out;
    FILE*               out;
    pthread_mutex_t     out_lock;

    const char*         query_db_title;
    int                 query_vol_index;
    CSeqDB*             query_vol;
    
    const char*         subject_db_title;
    int                 subject_vol_index;
    CSeqDB*             subject_vol;

    const HbnProgramOptions*  opts;

    BOOL                query_and_subject_are_the_same;
    LookupTable*        lktbl;
    WordFindData**      word_data_array;
    HbnSubseqHitExtnData** hit_extn_data_array;
    HbnHSPResults**     results_array;
} hbn_task_struct;

hbn_task_struct*
hbn_task_struct_new(const HbnProgramOptions* opts);

hbn_task_struct*
hbn_task_struct_free(hbn_task_struct* ht_struct);

void
hbn_task_struct_destroy_query_vol_context(hbn_task_struct* ht_struct);

void
hbn_task_struct_build_query_vol_context(hbn_task_struct* ht_struct, int query_vol_index);

void
hbn_task_struct_destroy_subject_vol_context(hbn_task_struct* ht_struct);

void
hbn_task_struct_build_subject_vol_context(hbn_task_struct* ht_struct, int subject_vol_index);

#ifdef __cplusplus
}
#endif

#endif // __HBN_TASK_STRUCT_H