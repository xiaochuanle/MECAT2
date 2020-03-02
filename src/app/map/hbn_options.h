#ifndef __HBN_OPTIONS_H
#define __HBN_OPTIONS_H

#include "../../ncbi_blast/setup/blast_options.h"
#include "../../ncbi_blast/setup/blast_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define HBN_QUERY_CHUNK_SIZE    20

typedef enum {
    eSeqid,
    eSeqidx,
    eSubseq,
    eSubseqx,
    eM4,
    eM4x,
    ePaf,
    eSAM,
    eInvalidFmt
} EOutputFormat;

extern const char* outfmt_names[eInvalidFmt];

EOutputFormat
string_to_outfmt(const char* str);

typedef enum {
    eHbnTask_pm,
    eHbnTask_rm,
    eHbnInvalidTask
} EHbnTask;

extern const char* hbn_task_names[eHbnInvalidTask];

EHbnTask
string_to_task(const char* str);

typedef struct {
    /// genearl search options
    EHbnTask        align_task;
    double          expect_value;
    int             penalty;
    int             reward;
    int             gap_open;
    int             gap_extend;

    /// input query options
    int             strand;

    /// database options
    const char*     db_dir;
    int             keep_db;
    int             min_query_size;
    int             max_query_vol_seqs;
    size_t          max_query_vol_res;
    int             min_subject_size;
    int             max_subject_vol_seqs;
    size_t          max_subject_vol_res;

    /// ddf scoring options
    int             kmer_size;
    int             kmer_window_size;
    int             max_kmer_occ;
    int             block_size;
    int             min_ddfs;

    /// mem scoring options
    int             memsc_kmer_size;
    int             memsc_kmer_window;
    int             memsc_mem_size;
    int             memsc_score;

    /// formatting options
    EOutputFormat   outfmt;
    int             dump_cigar;
    int             dump_md;

    /// query filtering options
    int             use_dust_masker;
    int             dust_masker_level;
    int             dust_masker_window;
    int             dust_masker_linker;
    int             use_lower_case_masker;

    /// restrict search or results
    double          perc_identity;
    double          query_cov_hsp_perc;
    int             query_cov_hsp_res;
    int             max_hsps_per_subject;
    int             hitlist_size;
    int             keep_best_hsp_per_subject;
    int             skip_overhang;

    /// statistical options
    i64             searchsp_eff;
    i64             db_length;
    int             dbseq_num;

    /// misc options
    int             num_threads;
    int             node_id;
    int             num_nodes;

    const char*     query;
    const char*     subject;
    const char*     output;
} HbnProgramOptions;

HbnProgramOptions*
HbnProgramOptionsFree(HbnProgramOptions* opts);

#ifdef __cplusplus
}
#endif

#endif // __HBN_OPTIONS_H