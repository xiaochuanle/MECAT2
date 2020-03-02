#ifndef __HBN_EXTEND_SUBSEQ_HIT_H
#define __HBN_EXTEND_SUBSEQ_HIT_H

#include "hbn_options.h"
#include "hbn_subseq_hit.h"
#include "../../algo/diff_gapalign.h"
#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vec_subseq_hit fwd_sbjct_subseq_list;
    vec_subseq_hit rev_sbjct_subseq_list;
    vec_subseq_hit sbjct_subseq_list;
    DiffGapAlignData* diff_data;
} HbnSubseqHitExtnData;

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataNew(const HbnProgramOptions* opts);

HbnSubseqHitExtnData*
HbnSubseqHitExtnDataFree(HbnSubseqHitExtnData* data);

void
hbn_extend_query_subseq_hit_list(HbnSubseqHit* subseq_hit_array,
    int subseq_hit_count,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_id,
    const int query_length,
    const text_t* db,
    const HbnProgramOptions* opts,
    HbnSubseqHitExtnData* data,
    BlastHitList* hit_list,
    HbnHSPResults* results);

#ifdef __cplusplus
}
#endif

#endif // __HBN_EXTEND_SUBSEQ_HIT_H