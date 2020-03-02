#ifndef __CNS_AUX_H
#define __CNS_AUX_H

#include "cns_options.h"
#include "raw_reads_reader.h"
#include "fccns.h"
#include "../../algo/chain_dp.h"
#include "../../algo/diff_gapalign.h"
#include "../../corelib/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_CNS_OVLPS   60
#define MAX_CNS_COV 15

typedef struct {
    int oid;
    int raw_read_size;
    int cns_from;
    int cns_to;
    int cns_read_size;
    size_t can_from;
    size_t can_to;
    size_t cns_fasta_offset;
    size_t cns_fasta_size;
} RawReadCnsInfo;

typedef struct {
    int start, end;
} MappingRange;

typedef struct {
    const HbnProgramOptions* opts;
    RawReadsReader* raw_reads;
    RawReadCnsInfo* cns_info_array;
    int cns_info_count;
    int* cns_info_idx;
    pthread_mutex_t* cns_info_idx_lock;
    HbnConsensusInitHit* cns_hit_array;
    size_t cns_hit_count;
    kstring_t qaux;
    kstring_t saux;
    vec_u8 fwd_read;
    vec_u8 rev_read;
    vec_u8 fwd_subject;
    vec_u8 rev_subject;
    kstring_t* cns_out;
    pthread_mutex_t* cns_out_lock;
    vec_u8 cov_stats;
    FCCnsData* cns_data;
    DiffGapAlignData* diff_data;
    Ksw2Data* ksw;
} CnsThreadData;

CnsThreadData*
CnsThreadDataNew(const HbnProgramOptions* opts,
    RawReadsReader* raw_reads,
    RawReadCnsInfo* cns_info_array,
    int* cns_info_idx,
    pthread_mutex_t* cns_info_idx_lock,
    kstring_t* cns_out,
    pthread_mutex_t* out_lock);

CnsThreadData*
CnsThreadDataFree(CnsThreadData* data);

void
CnsThreadDataInit(CnsThreadData* data, 
    const u8* fwd_subject,
    const u8* rev_subject,
    const int subject_size);

void 
normalize_gaps(const char* qstr, 
    const char* tstr, 
    const int aln_size, 
    kstring_t* qnorm, 
    kstring_t* tnorm, 
    const BOOL push);
    
HbnConsensusInitHit*
load_and_sort_cns_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    size_t* hit_count);

BOOL
set_next_raw_read_batch_info(RawReadCnsInfo* cns_info_array,
    int* cns_info_count,
    HbnConsensusInitHit* cns_hit_array,
    size_t cns_hit_count,
    size_t* cns_hit_idx,
    RawReadsReader* raw_reads,
    const int batch_size);

#ifdef __cplusplus
}
#endif

#endif // __CNS_AUX_H