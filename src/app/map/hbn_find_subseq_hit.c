#include "hbn_find_subseq_hit.h"

#include "../../corelib/khash.h"
#include "../../corelib/ksort.h"

KHASH_MAP_INIT_INT(map_int_int, int);

typedef struct {
    int oid;
    int offset;
    int count;
    int score;
} SubjectInitHitInfo;

typedef kvec_t(SubjectInitHitInfo) vec_sbjct_hit_info;

#define sbjct_init_hit_info_score_gt(lhs, rhs) ((lhs).score > (rhs).score)
KSORT_INIT(sbjct_init_hit_info_score_gt, SubjectInitHitInfo, sbjct_init_hit_info_score_gt);

static int
group_subject_init_hits(HbnInitHit* init_hit_array, int init_hit_count, const int hitlist_size_max)
{
    ks_introsort_init_hit_sid_lt(init_hit_count, init_hit_array);
    kv_dinit(vec_sbjct_hit_info, sbjct_hit_info_list);
    SubjectInitHitInfo sbjct_hit_info;
    int i = 0;
    while (i < init_hit_count) {
        int j = i + 1;
        while (j < init_hit_count && init_hit_array[i].sid == init_hit_array[j].sid) ++j;
        sbjct_hit_info.offset = i;
        sbjct_hit_info.count = j - i;
        sbjct_hit_info.score = 0;
        for (int k = i; k < j; ++k) sbjct_hit_info.score = hbn_max(sbjct_hit_info.score, init_hit_array[k].score);
        kv_push(SubjectInitHitInfo, sbjct_hit_info_list, sbjct_hit_info);
        i = j;
    }

    ks_introsort_sbjct_init_hit_info_score_gt(kv_size(sbjct_hit_info_list), kv_data(sbjct_hit_info_list));
    HbnInitHit* tmp_init_hit_array = (HbnInitHit*)malloc(sizeof(HbnInitHit) * init_hit_count);
    int init_hit_idx = 0;
    for (i = 0; i < hitlist_size_max && i < kv_size(sbjct_hit_info_list); ++i) {
        HbnInitHit* c = init_hit_array + kv_A(sbjct_hit_info_list, i).offset;
        int n = kv_A(sbjct_hit_info_list, i).count;
        memcpy(tmp_init_hit_array + init_hit_idx, c, sizeof(HbnInitHit) * n);
        init_hit_idx += n;
    }
    memcpy(init_hit_array, tmp_init_hit_array, sizeof(HbnInitHit) * init_hit_idx);
    free(tmp_init_hit_array);
    kv_destroy(sbjct_hit_info_list);
    return init_hit_idx;
}

static void
adjust_init_hit_subject_offset(HbnInitHit* hit, const EHbnTask task, const int query_size, const size_t subject_size, size_t* sfrom, size_t* sto)
{
    if (task == eHbnTask_pm
        ||
        subject_size <= query_size
        ||
        subject_size - query_size <= 2000) {
        *sfrom = 0;
        *sto = subject_size;
        return;
    }

    size_t ql = hit->qoff;
    size_t sl = hit->soff;
    size_t ll = hbn_min(ql, sl);
    size_t from = hit->soff - ll;
    int E = ll * 0.3;
    if (E > 2000) E = 2000;
    from = (from <= E) ? 0 : (from - E);

    size_t qr = query_size - hit->qoff;
    size_t sr = subject_size - hit->soff;
    size_t rl = hbn_min(qr, sr);
    size_t to = hit->soff + rl;
    E = rl * 0.3;
    if (E > 2000) E = 2000;
    to = (subject_size - to <= E) ? subject_size : (to + E);

    *sfrom = from;
    *sto = to;
}

void
setup_candidates_for_one_subject(HbnInitHit* init_hit_array,
    int init_hit_count,
    const int query_id,
    const int query_size,
    const size_t subject_size,
    const HbnProgramOptions* opts,
    const text_t* db,
    vec_subseq_hit* fwd_sbjct_subseq_list,
    vec_subseq_hit* rev_sbjct_subseq_list,
    vec_subseq_hit* sbjct_subseq_list)
{
    kv_clear(*fwd_sbjct_subseq_list);
    kv_clear(*rev_sbjct_subseq_list);

    for (int i = 0; i < init_hit_count; ++i) {
        HbnInitHit* hit = init_hit_array + i;
        //HBN_LOG("process\t");
       // dump_init_hit(fprintf, stderr, *hit);
        hbn_assert(hit->sid == init_hit_array[0].sid);
        hbn_assert(hit->sdir == FWD);
        size_t sfrom = 0, sto = 0;
        adjust_init_hit_subject_offset(hit, opts->align_task, query_size, subject_size, &sfrom, &sto);
        HbnSubseqHit subseq;
        subseq.qid = query_id;
        subseq.qdir = hit->qdir;
        subseq.qoff = hit->qoff;
        subseq.sid = hit->sid;
        subseq.soff = hit->soff - sfrom;
        subseq.sfrom = sfrom;
        subseq.sto = sto;
        subseq.score = hit->score;
        //fprintf(stderr, "adding\t");
        //dump_can_seg(fprintf, stderr, subseq);
        vec_subseq_hit* sbject_subseq_list = (subseq.qdir == FWD) ? fwd_sbjct_subseq_list : rev_sbjct_subseq_list;
        kv_push(HbnSubseqHit, *sbject_subseq_list, subseq);
        //dump_can_sbjct_subseq(fprintf, stderr, subseq);
    }

    //merge_adjacent_subject_subseqs(fwd_sbjct_subseq_list, opts->align_task);
    //merge_adjacent_subject_subseqs(rev_sbjct_subseq_list, opts->align_task);
    //kv_clear(*sbjct_subseq_list);
    for (size_t i = 0; i < kv_size(*fwd_sbjct_subseq_list); ++i) {
        HbnSubseqHit subseq = kv_A(*fwd_sbjct_subseq_list, i);
        kv_push(HbnSubseqHit, *sbjct_subseq_list, subseq);
    }
    for (size_t i = 0; i < kv_size(*rev_sbjct_subseq_list); ++i) {
        HbnSubseqHit subseq = kv_A(*rev_sbjct_subseq_list, i);
        kv_push(HbnSubseqHit, *sbjct_subseq_list, subseq);
    }
}

void
find_candidate_subject_subseqs(HbnInitHit* init_hit_array,
    int init_hit_count,
    const int query_id,
    const int query_size,
    const HbnProgramOptions* opts,
    const text_t* db,
    vec_subseq_hit* fwd_sbjct_subseq_list,
    vec_subseq_hit* rev_sbjct_subseq_list,
    vec_subseq_hit* sbjct_subseq_list)
{
    kv_clear(*sbjct_subseq_list);
    if (init_hit_count == 0) return;
    init_hit_count = group_subject_init_hits(init_hit_array, init_hit_count, opts->hitlist_size);

    int i = 0;
    int hitlist_count = 0;
    while (i < init_hit_count && hitlist_count < opts->hitlist_size) {
        int j = i + 1;
        while (j < init_hit_count && init_hit_array[i].sid == init_hit_array[j].sid) ++j;
        int subject_id = init_hit_array[i].sid;
        size_t subject_size = seqdb_seq_size(db, subject_id);
        setup_candidates_for_one_subject(init_hit_array + i,
            j - i,
            query_id,
            query_size,
            subject_size,
            opts,
            db,
            fwd_sbjct_subseq_list,
            rev_sbjct_subseq_list,
            sbjct_subseq_list);
        i = j;
    }
}