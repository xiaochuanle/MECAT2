#include "gapped_candidate.h"

#include "ksort.h"

#define init_hit_score_gt(a, b) ((a).score > (b).score)
KSORT_INIT(init_hit_score_gt, HbnInitHit, init_hit_score_gt);

#define init_hit_qid_lt(a, b) ((a).qid < (b).qid)
KSORT_INIT(init_hit_qid_lt, HbnInitHit, init_hit_qid_lt);

#define init_hit_sid_lt(a, b) ((a).sid < (b).sid)
KSORT_INIT(init_hit_sid_lt, HbnInitHit, init_hit_sid_lt);

///////////////

#define cns_hit_sid_lt(a, b) ((a).sid < (b).sid)
KSORT_INIT(cns_hit_sid_lt, HbnConsensusInitHit, cns_hit_sid_lt);

#define cns_hit_score_gt(a, b) ( \
    ((a).score > (b).score) \
    || \
    ((a).score == (b).score && (a).sid < (b).sid) \
)
KSORT_INIT(cns_hit_score_gt, HbnConsensusInitHit, cns_hit_score_gt);