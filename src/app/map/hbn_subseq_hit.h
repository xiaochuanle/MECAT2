#ifndef __HBN_SUBSEQ_HIT_H
#define __HBN_SUBSEQ_HIT_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    int qoff;
    int sid;
    int soff;
    size_t sfrom;
    size_t sto;
    int score;
    int chain_seed_offset;
    int chain_seed_count;
} HbnSubseqHit;

#define dump_subseq_hit(output_func, stream, hit) \
    output_func(stream, "qdir = %d, qoff = %d, sid = %d, soff = %d, sfrom = %zu, sto = %zu, score = %d\n", \
        (hit).qdir, \
        (hit).qoff, \
        (hit).sid, \
        (hit).soff, \
        (hit).sfrom, \
        (hit).sto, \
        (hit).score)

typedef kvec_t(HbnSubseqHit) vec_subseq_hit;

void ks_introsort_subseq_hit_score_gt(size_t n, HbnSubseqHit* a);

void ks_introsort_subseq_hit_sid_lt(size_t n, HbnSubseqHit* a);

void ks_introsort_subseq_hit_sfrom_lt(size_t n, HbnSubseqHit* a);

#ifdef __cplusplus
}
#endif

#endif // __HBN_SUBSEQ_HIT_H