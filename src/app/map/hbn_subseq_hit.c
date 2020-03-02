#include "hbn_subseq_hit.h"

#include "../../corelib/ksort.h"

#define subseq_hit_score_gt(a, b) ((a).score > (b).score)
KSORT_INIT(subseq_hit_score_gt, HbnSubseqHit, subseq_hit_score_gt);

#define subseq_hit_sfrom_lt(a, b) (((a).sfrom < (b).sfrom) || ((a).sfrom == (b).sfrom && (a).sto > (b).sto))
KSORT_INIT(subseq_hit_sfrom_lt, HbnSubseqHit, subseq_hit_sfrom_lt);