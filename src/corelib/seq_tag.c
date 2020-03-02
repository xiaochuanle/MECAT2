#include "seq_tag.h"

#include "string2hsp.h"
#include "ksort.h"

static const int kDumpRepeatHsp = 0;
static int RepeatSeq = 0;
static const int kDumpCircularHsp = 0;
static int CircularSeq = 0;
static const int kDumpTandemRepeatHsp = 0;
static int TandemRepeatSeq = 0;
static const int kDumpRcDup = 1;
static int RcDupSeq = 0;

static const int max_e = 200;
static const int MAX_E = 200;

SeqTag* g_seq_tag = NULL;

const char* sSeqTagNames[] = {
    "perfect",
    "unperfect",
    "repeat",
    "chimeric",
    "circular",
    "tandem-repeat",
    "rcdup"
};

const char* GetSeqTagName(int tag)
{
    return sSeqTagNames[tag];
}

int SeqNameToTag(const char* name)
{
    int tag = SEQ_TAG_MAX;
    for (int i = 0; i < SEQ_TAG_MAX; ++i) {
        if (strcmp(name, sSeqTagNames[i]) == 0) {
            tag = i;
            break;
        }
    }
    if (tag == SEQ_TAG_MAX) HBN_ERR("Invalid sequence tag: %s", name);
    return tag;
}

//////////////////

BOOL is_contained_hsp(const HbnHSP* hsp)
{
    size_t qoff = hsp->qoff;
    size_t qend = hsp->qend;
    size_t soff = hsp->soff;
    size_t send = hsp->send;
    
    BOOL r = (qoff <= max_e) && (hsp->qsize - qend <= max_e);
    if (r) return r;
    r = (soff <= max_e) && (hsp->ssize - send <= max_e);
    return r;
}

static seq_tag_t
examine_perfect_and_repeat_read(HbnHSP* hsp_array, int hsp_count)
{
    HbnHSP perfect_hsp_array[hsp_count];
    int perfect_hsp_count = 0;
    for (int i = 0; i < hsp_count; ++i) {
        if (is_contained_hsp(hsp_array + i)) {
            perfect_hsp_array[perfect_hsp_count++] = hsp_array[i];
        }
    }
    if (perfect_hsp_count == 0) return SEQ_TAG_MAX;
    if (perfect_hsp_count == 1) {
        if (g_seq_tag) {
            g_seq_tag->tag = SEQ_TAG_PERFECT;
            g_seq_tag->ident_perc = hsp_array[0].ident_perc;
            g_seq_tag->sid = hsp_array[0].sid;
            g_seq_tag->qoff = 0;
            g_seq_tag->qend = hsp_array[0].qsize;
            g_seq_tag->soff = hsp_array[0].soff;
            g_seq_tag->send = hsp_array[0].send;
        }
        return SEQ_TAG_PERFECT;
    }

    ks_introsort_hbn_hsp_ident_gt(perfect_hsp_count, perfect_hsp_array);
    HbnHSP* hsp1 = &perfect_hsp_array[0];
    HbnHSP* hsp2 = &perfect_hsp_array[1];
    if (hsp1->ident_perc - hsp2->ident_perc > 8.0) return SEQ_TAG_PERFECT;

    if (kDumpRepeatHsp) {
        HBN_LOG("find repeat hsp");
        ks_dinit(hspstr);
        for (int i = 0; i < perfect_hsp_count; ++i) {
            hsp2string(&perfect_hsp_array[i], &hspstr);
            fprintf(stderr, "%s\n", ks_s(hspstr));
        }
        ks_destroy(hspstr);
    }
    return SEQ_TAG_REPEAT;
}

static seq_tag_t
s_examine_circular_hsp_one_subject(HbnHSP* hsp_array, int hsp_count)
{
    if (hsp_count < 2) return SEQ_TAG_MAX;

    HbnHSP lhspv[hsp_count];
    int lhsp = 0;
    HbnHSP rhspv[hsp_count];
    int rhsp = 0;
    for (int i = 0; i < hsp_count; ++i) {
        HbnHSP* hsp = hsp_array + i;
        if (hsp->soff <= max_e && hsp->qsize - hsp->qend <= max_e) {
            lhspv[lhsp++] = *hsp;
        }
        if (hsp->ssize - hsp->send <= max_e && hsp->qoff <= max_e) {
            rhspv[rhsp++] = *hsp;
        }
    }
    if (lhsp == 0 || rhsp == 0) return SEQ_TAG_MAX;

    for (int i = 0; i < lhsp; ++i) {
        HbnHSP* lh = &lhspv[i];
        hbn_assert(lh->sdir == FWD);
        for (int j = 0; j < rhsp; ++j) {
            HbnHSP* rh = &rhspv[j];
            hbn_assert(rh->sdir == FWD);
            if (lh->qdir != rh->qdir) continue;
            if (lh->qoff + max_e >= rh->qend || rh->qend + max_e >= lh->qoff) {
                if (kDumpCircularHsp) {
                    HBN_LOG("find circular hsp");
                    ks_dinit(hspstr);
                    hsp2string(lh, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    hsp2string(rh, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                }
                if (g_seq_tag) {
                    g_seq_tag->sid = lh->sid;
                    g_seq_tag->tag = SEQ_TAG_CIRCULAR;
                    g_seq_tag->qoff = rh->qoff;
                    g_seq_tag->qend = rh->qend;
                    g_seq_tag->soff = rh->soff;
                    g_seq_tag->send = rh->send;
                    if (rh->qdir == REV) {
                        g_seq_tag->qoff = rh->qsize - rh->qend;
                        g_seq_tag->qend = rh->qsize - rh->qoff;
                    }
                }
                return SEQ_TAG_CIRCULAR;
            }
        }
    }

    return SEQ_TAG_MAX;
}

static seq_tag_t
examine_circular_hsp(HbnHSP* hsp_array, int hsp_count)
{
    ks_introsort_hbn_hsp_sid_lt(hsp_count, hsp_array);
    int i = 0;
    while (i < hsp_count) {
        int j = i + 1;
        while (j < hsp_count && hsp_array[i].sid == hsp_array[j].sid) ++j;
        seq_tag_t tag = s_examine_circular_hsp_one_subject(hsp_array + i, j - i);
        if (tag == SEQ_TAG_CIRCULAR) return tag;
        i = j;
    }

    return SEQ_TAG_MAX;
}

#define hbn_hsp_send_gt(a, b) ((a).send > (b).send)
KSORT_INIT(hbn_hsp_send_gt, HbnHSP, hbn_hsp_send_gt);

static seq_tag_t
s_examine_tandem_repeat_hsp_qdir_one_subject(HbnHSP* hsp_array, int hsp_count)
{
    if (hsp_count < 2) return SEQ_TAG_MAX;
    ks_introsort_hbn_hsp_soff_lt(hsp_count, hsp_array);
    int i = 0;
    int* cov = (int*)malloc(sizeof(int) * hsp_array[0].qsize);
    while (i < hsp_count) {
        size_t soff = hsp_array[i].soff;
        int j = i + 1;
        while (j < hsp_count && soff + max_e >= hsp_array[j].soff) ++j;
        if (j - i > 1) {
            memset(cov, 0, sizeof(int) * hsp_array[0].qsize);
            for (int k = i; k < j; ++k) {
                for (size_t pos = hsp_array[k].qoff; pos < hsp_array[k].qend; ++pos) cov[pos] = 1;
            }
            size_t n = 0;
            for (size_t pos = 0; pos < hsp_array[0].qsize; ++pos) n += cov[pos];
            if (hsp_array[0].qsize - n <= max_e) {
                if (kDumpTandemRepeatHsp) {
                    HBN_LOG("find tandem repeat hsp");
                    ks_dinit(hspstr);
                    for (int k = i; k < j; ++k) {
                        hsp2string(hsp_array + k, &hspstr);
                        fprintf(stderr, "%s\n", ks_s(hspstr));
                    }
                    ks_destroy(hspstr);
                }
                if (j - i == 2) {
                    HbnHSP* h1 = hsp_array + i;
                    HbnHSP* h2 = hsp_array + i + 1;
                    HbnHSP* hsp = (h1->send - h1->soff > h2->send - h2->soff) ? h1 : h2;

#if 0
                    HBN_LOG("two fold tandem repeat");
                    ks_dinit(hspstr);
                    hsp2string(h1, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    hsp2string(h2, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    ks_destroy(hspstr);
#endif

                    g_seq_tag->tag = SEQ_TAG_TANDEM_REPEAT;
                    g_seq_tag->ident_perc = hsp->ident_perc;
                    g_seq_tag->sid = hsp->sid;
                    g_seq_tag->qoff = hsp->qoff;
                    g_seq_tag->qend = hsp->qend;
                    g_seq_tag->soff = hsp->soff;
                    g_seq_tag->send = hsp->send;
                    if (hsp->qdir == REV) {
                        g_seq_tag->qoff = hsp->qsize - hsp->qend;
                        g_seq_tag->qend = hsp->qsize - hsp->qoff;
                    }
                }
                ++TandemRepeatSeq;
                return SEQ_TAG_TANDEM_REPEAT;
            }
        }
        i = j;
    }

    ks_introsort_hbn_hsp_send_gt(hsp_count, hsp_array);
    i = 0;
    while (i < hsp_count) {
        size_t send = hsp_array[i].send;
        int j = i + 1;
        while (j < hsp_count && send <= hsp_array[j].send + max_e) ++j;
        if (j - i > 1) {
            memset(cov, 0, sizeof(int) * hsp_array[0].qsize);
            for (int k = i; k < j; ++k) {
                for (size_t pos = hsp_array[k].qoff; pos < hsp_array[k].qend; ++pos) cov[pos] = 1;
            }
            size_t n = 0;
            for (size_t pos = 0; pos < hsp_array[0].qsize; ++pos) n += cov[pos];
            if (hsp_array[0].qsize - n <= max_e) {
                if (kDumpTandemRepeatHsp) {
                    HBN_LOG("find tandem repeat hsp %d", TandemRepeatSeq);
                    ks_dinit(hspstr);
                    for (int k = i; k < j; ++k) {
                        hsp2string(hsp_array + k, &hspstr);
                        fprintf(stderr, "%s\n", ks_s(hspstr));
                    }
                }
                ++TandemRepeatSeq;
                if (j - i == 2) {
                    HbnHSP* h1 = hsp_array + i;
                    HbnHSP* h2 = hsp_array + i + 1;
                    HbnHSP* hsp = (h1->send - h1->soff > h2->send - h2->soff) ? h1 : h2;
#if 0
                    HBN_LOG("two fold tandem repeat");
                    ks_dinit(hspstr);
                    hsp2string(h1, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    hsp2string(h2, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    ks_destroy(hspstr);
#endif
                    
                    g_seq_tag->tag = SEQ_TAG_TANDEM_REPEAT;
                    g_seq_tag->ident_perc = hsp->ident_perc;
                    g_seq_tag->sid = hsp->sid;
                    g_seq_tag->qoff = hsp->qoff;
                    g_seq_tag->qend = hsp->qend;
                    g_seq_tag->soff = hsp->soff;
                    g_seq_tag->send = hsp->send;
                    if (hsp->qdir == REV) {
                        g_seq_tag->qoff = hsp->qsize - hsp->qend;
                        g_seq_tag->qend = hsp->qsize - hsp->qoff;
                    }
                }
                return SEQ_TAG_TANDEM_REPEAT;
            }
        }
        i = j;
    }   

    free(cov);
    return SEQ_TAG_MAX;
}

static seq_tag_t
s_examine_tandem_repeat_hsp_one_subject(HbnHSP* hsp_array, int hsp_count)
{
    HbnHSP t_hsp_array[hsp_count];
    ks_introsort_hbn_hsp_sid_lt(hsp_count, hsp_array);
    int i = 0;
    while (i < hsp_count) {
        int j = i + 1;
        while (j < hsp_count && hsp_array[i].sid == hsp_array[j].sid) ++j;
        int n = 0;
        for (int k = i; k < j; ++k) {
            if (hsp_array[k].qdir == FWD) t_hsp_array[n++] = hsp_array[k];
        }
        seq_tag_t tag = s_examine_tandem_repeat_hsp_qdir_one_subject(t_hsp_array, n);
        if (tag == SEQ_TAG_TANDEM_REPEAT) return tag;
        n = 0;
        for (int k = i; k < j; ++k) {
            if (hsp_array[k].qdir == REV) t_hsp_array[n++] = hsp_array[k];
        }
        tag = s_examine_tandem_repeat_hsp_qdir_one_subject(t_hsp_array, n);
        if (tag == SEQ_TAG_TANDEM_REPEAT) return tag;

        i = j;
    }
    return SEQ_TAG_MAX;
}

static seq_tag_t
examine_tandem_repeat_hsp(HbnHSP* hsp_array, int hsp_count)
{
    ks_introsort_hbn_hsp_sid_lt(hsp_count, hsp_array);
    int i = 0;
    while (i < hsp_count) {
        int j = i + 1;
        while (j < hsp_count && hsp_array[i].sid == hsp_array[j].sid) ++j;
        seq_tag_t tag = s_examine_tandem_repeat_hsp_one_subject(hsp_array + i, j - i);
        if (tag == SEQ_TAG_TANDEM_REPEAT) return SEQ_TAG_TANDEM_REPEAT;
        i = j;
    }
    return SEQ_TAG_MAX;
}

static seq_tag_t
s_examine_rcdup_hsp_one_subject(HbnHSP* hsp_array, int hsp_count)
{
    if (hsp_count < 2) return SEQ_TAG_MAX;

    HbnHSP fhspv[hsp_count];
    HbnHSP rhspv[hsp_count];
    int fhsp = 0;
    int rhsp = 0;
    for (int i = 0; i < hsp_count; ++i) {
        if (hsp_array[i].qdir == FWD) fhspv[fhsp++] = hsp_array[i];
        else rhspv[rhsp++] = hsp_array[i];
    }
    if (fhsp == 0 || rhsp == 0) return SEQ_TAG_MAX;

    int qsize = hsp_array[0].qsize;
    int* cov = (int*)malloc(sizeof(int) * qsize);
    for (int i = 0; i < fhsp; ++i) {
        HbnHSP* fh = fhspv + i;
        for (int j = 0; j < rhsp; ++j) {
            HbnHSP* rh = rhspv + j;
            int r = (fh->soff <= rh->soff + max_e && fh->send + max_e >= rh->send)
                    ||
                    (fh->soff + max_e >= rh->soff && fh->send <= rh->send + max_e);
            if (!r) continue;
            memset(cov, 0, sizeof(int) * qsize);
            for (size_t pos = fh->qoff; pos < fh->qend; ++pos) cov[pos] = 1;
            size_t qoff = qsize - rh->qend;
            size_t qend = qsize - rh->qoff;
            for (size_t pos = qoff; pos < qend; ++pos) cov[pos] = 1;
            int n = 0;
            for (size_t pos = 0; pos < qsize; ++pos) n += cov[pos];
            if (qsize - n <= max_e) {
                if (kDumpRcDup) {
                    HBN_LOG("find rcdup %d", RcDupSeq);
                    ks_dinit(hspstr);
                    hsp2string(fh, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));
                    hsp2string(rh, &hspstr);
                    fprintf(stderr, "%s\n", ks_s(hspstr));                    
                }
                if (g_seq_tag) {
                    HbnHSP* hsp = (rh->send - rh->soff > fh->send - fh->soff) ? rh : fh;
                    g_seq_tag->tag = SEQ_TAG_RCDUP;
                    g_seq_tag->ident_perc = hsp->ident_perc;
                    g_seq_tag->qoff = hsp->qoff;
                    g_seq_tag->qend = hsp->qend;
                    g_seq_tag->sid = hsp->sid;
                    g_seq_tag->soff = hsp->soff;
                    g_seq_tag->send = hsp->send;
                    if (hsp->qdir == REV) {
                        g_seq_tag->qoff = hsp->qsize - hsp->qend;
                        g_seq_tag->qend = hsp->qsize - hsp->qoff;
                    }
                }
                
                ++RcDupSeq;
                return SEQ_TAG_RCDUP;
            }
        }
    }
    free(cov);
    return SEQ_TAG_MAX;
}

static seq_tag_t
examine_rcdup_hsp(HbnHSP* hsp_array, int hsp_count)
{
    ks_introsort_hbn_hsp_sid_lt(hsp_count, hsp_array);
    int i = 0;
    while (i < hsp_count) {
        int j = i + 1;
        while (j < hsp_count && hsp_array[i].sid == hsp_array[j].sid) ++j;
        seq_tag_t tag = s_examine_rcdup_hsp_one_subject(hsp_array + i, j - i);
        if (tag == SEQ_TAG_RCDUP) return tag;
        i = j;
    }
    return SEQ_TAG_MAX;    
}

seq_tag_t TagSeqFromHSPArray(HbnHSP* hsp_array, int hsp_count, double min_ident_perc, SeqTag* seq_tag)
{
    int n = 0;
    for (int i = 0; i < hsp_count; ++i) {
        if (hsp_array[i].ident_perc + 8.0 < min_ident_perc) continue;
        HbnHSP* hsp = hsp_array + i;
        normalise_hbn_hsp_sdir(hsp);
        hbn_assert(hsp->sdir == FWD);
        if (hsp->qdir == REV) {
            size_t qoff = hsp->qsize - hsp->qend;
            size_t qend = hsp->qsize - hsp->qoff;
            hsp->qoff = qoff;
            hsp->qend = qend;
        }
        hsp_array[n++] = hsp_array[i];
    }
    hsp_count = n;

    g_seq_tag = seq_tag;
    seq_tag_t tag;
    tag = examine_perfect_and_repeat_read(hsp_array, hsp_count);
    if (tag == SEQ_TAG_PERFECT || tag == SEQ_TAG_REPEAT) return tag;

    tag = examine_circular_hsp(hsp_array, hsp_count);
    if (tag == SEQ_TAG_CIRCULAR) return tag;

    tag = examine_rcdup_hsp(hsp_array, hsp_count);
    if (tag == SEQ_TAG_RCDUP) return tag;

    tag = examine_tandem_repeat_hsp(hsp_array, hsp_count);
    if (tag == SEQ_TAG_TANDEM_REPEAT) return tag;

    return SEQ_TAG_UNPERFECT;
}

////////////

static BOOL
examine_ovlp_quality_normal(int qoff,
    int qend,
    int soff,
    int send,
    double min_ovlp_frac,
    const int read_size,
    const int subject_size)
{
    const int MAX_E = 200;
    BOOL r = (qend - qoff >= min_ovlp_frac / 100.0)
             ||
             (send - soff >= min_ovlp_frac / 100.0);
    if (r) return r;
    r = (read_size - qend <= MAX_E && soff <= MAX_E)
        ||
        (subject_size - send <= MAX_E && soff <= MAX_E);
    return r;
}

static BOOL
examine_ovlp_quality_perfect(int qoff,
    int qend,
    int soff,
    int send,
    const int read_size,
    const int subject_size)
{
    BOOL r = (qoff <= MAX_E || read_size - qend <= MAX_E)
             ||
             (soff <= MAX_E || subject_size - send <= MAX_E);
    if (r) return r;
    r = (read_size - qend <= MAX_E && soff <= MAX_E)
        ||
        (subject_size - send <= MAX_E && soff <= MAX_E);
    if (!r) return r;
    return r;
}

static BOOL
examine_ovlp_quality_unperfect(int qoff,
    int qend,
    int soff,
    int send,
    const double min_ovlp_frac,
    const int read_size,
    const int subject_size)
{
    BOOL r = (qend - qoff >= read_size * min_ovlp_frac / 100.0)
             ||
             (send - soff >= subject_size * min_ovlp_frac / 100.0);
    return r;
}

BOOL examine_ovlp_quality(const SeqTag* tag_array,
    int qoff,
    int qend,
    int soff,
    int send,
    double min_ovlp_frac,
    const int read_id,
    const int read_size,
    const int subject_id,
    const int subject_size)
{
    if (!tag_array) return examine_ovlp_quality_normal(qoff, qend, soff, send, min_ovlp_frac, read_size, subject_size);
    SeqTag stag = tag_array[subject_id];
    SeqTag qtag = tag_array[read_id];
    if (stag.tag == SEQ_TAG_UNPERFECT) {
        return examine_ovlp_quality_unperfect(qoff, qend, soff, send, 0.8, read_size, subject_size);
    } else if (stag.tag == SEQ_TAG_PERFECT && qtag.tag == SEQ_TAG_UNPERFECT) {
        return examine_ovlp_quality_unperfect(qoff, qend, soff, send, 0.8, read_size, subject_size);
    }
    return TRUE;
}

/////////////

static BOOL
s_examine_perfect_subject_hit(const SeqTag* tag_array, const HbnConsensusInitHit* hit)
{
    SeqTag qtag = tag_array[hit->qid];
    if (qtag.tag == SEQ_TAG_UNPERFECT) return TRUE;

    SeqTag stag = tag_array[hit->sid];
    if (qtag.tag == SEQ_TAG_PERFECT) {
        if (stag.sid != qtag.sid) return FALSE;
        if (stag.soff + MAX_E >= qtag.soff && stag.send <= qtag.send + MAX_E) return TRUE;
        if (qtag.soff + MAX_E >= stag.soff && qtag.send <= stag.send + MAX_E) return TRUE;
        if (qtag.soff <= stag.soff && qtag.send >= stag.soff) return TRUE;
        if (stag.soff <= qtag.soff && stag.send >= qtag.soff) return TRUE;
        return FALSE;
    }

    if (qtag.tag == SEQ_TAG_CIRCULAR) {
        if (stag.send <= qtag.qsize - qtag.qend + MAX_E) return TRUE;
        if (stag.send + MAX_E >= qtag.soff) return TRUE;
        return FALSE;
    }

    return FALSE;
}

static BOOL
s_examine_circular_subject_hit(const SeqTag* tag_array, const HbnConsensusInitHit* hit)
{
    SeqTag qtag = tag_array[hit->qid];
    if (qtag.tag == SEQ_TAG_UNPERFECT) return FALSE;
    if (qtag.tag == SEQ_TAG_CIRCULAR) return TRUE;

    if (qtag.tag == SEQ_TAG_PERFECT) {
        SeqTag stag = tag_array[hit->sid];
        if (qtag.send <= stag.qsize - stag.qend + MAX_E) return TRUE;
        if (qtag.send + MAX_E >= stag.soff) return TRUE;
        return FALSE;
    }

    return FALSE;
}

static BOOL
s_examine_unperfect_subject_hit(const SeqTag* tag_array, const HbnConsensusInitHit* hit)
{
    return TRUE;
}

BOOL
query_tag_is_valid(const SeqTag* tag_array, const HbnConsensusInitHit* hit)
{
    if (!tag_array) return TRUE;

    SeqTag stag = tag_array[hit->sid];
    if (stag.tag == SEQ_TAG_PERFECT) return s_examine_perfect_subject_hit(tag_array, hit);
    if (stag.tag == SEQ_TAG_CIRCULAR) return s_examine_circular_subject_hit(tag_array, hit);
    if (stag.tag == SEQ_TAG_UNPERFECT) return s_examine_unperfect_subject_hit(tag_array, hit);
    hbn_assert(0);
    return FALSE;
}