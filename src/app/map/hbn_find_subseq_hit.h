#ifndef __HBN_FIND_SUBSEQ_HIT_H
#define __HBN_FIND_SUBSEQ_HIT_H

#include "hbn_options.h"
#include "hbn_subseq_hit.h"
#include "../../algo/word_finder.h"

#ifdef __cplusplus
extern "C" {
#endif

void
find_candidate_subject_subseqs(HbnInitHit* init_hit_array,
    int init_hit_count,
    const int query_id,
    const int query_size,
    const HbnProgramOptions* opts,
    const text_t* db,
    vec_subseq_hit* fwd_sbjct_subseq_list,
    vec_subseq_hit* rev_sbjct_subseq_list,
    vec_subseq_hit* sbjct_subseq_list);

#ifdef __cplusplus
}
#endif

#endif // __HBN_FIND_SUBSEQ_HIT_H