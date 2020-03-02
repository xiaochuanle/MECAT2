#include "hbn_hit.h"

#include "ksort.h"

#define hbn_hsp_soff_lt(a, b) ( \
    ((a).soff < (b).soff) \
    || \
    ((a).soff == (b).soff && (a).send > (b).send) \
)
KSORT_INIT(hbn_hsp_soff_lt, HbnHSP, hbn_hsp_soff_lt);

#define hbn_hsp_ident_gt(a, b) ((a).ident_perc > (b).ident_perc)
KSORT_INIT(hbn_hsp_ident_gt, HbnHSP, hbn_hsp_ident_gt);

#define hbn_hsp_score_gt(a, b) ((a).score > (b).score)
KSORT_INIT(hbn_hsp_score_gt, HbnHSP, hbn_hsp_score_gt);

#define hbn_hsp_list_score_gt(a, b) ( \
    ((a).best_score > (b).best_score) \
    || \
    ((a).best_score == (b).best_score && (a).oid < (b).oid) \
)
KSORT_INIT(hbn_hsp_list_score_gt, HbnHSPList, hbn_hsp_list_score_gt);

#define hbn_hsp_sid_lt(a, b) ((a).sid < (b).sid)
KSORT_INIT(hbn_hsp_sid_lt, HbnHSP, hbn_hsp_sid_lt);

void normalise_hbn_hsp_offsets(HbnHSP* hsp)
{
    if (hsp->use_normalised_offset) return;

    if (hsp->qdir == REV) {
        size_t off = hsp->qsize - hsp->qend;
        size_t end = hsp->qsize - hsp->qoff;
        hsp->qoff = off;
        hsp->qend = end;
    }

    if (hsp->sdir == REV) {
        size_t off = hsp->ssize - hsp->send;
        size_t end = hsp->ssize - hsp->soff;
        hsp->soff = off;
        hsp->send = end;
    }

    hsp->use_normalised_offset = 1;
}

void normalise_hbn_hsp_sdir(HbnHSP* hsp)
{
    if (hsp->sdir == REV) {
        hsp->sdir = CMP_STRAND(hsp->sdir);
        hsp->qdir = CMP_STRAND(hsp->qdir);
    } else {
        hbn_assert(hsp->sdir == FWD);
    }
}

HbnAlignResults*
HbnAlignResultsNew(const int hitlist_max)
{
    HbnAlignResults* retval = (HbnAlignResults*)calloc(1, sizeof(HbnAlignResults));
    retval->num_queries = 0;
    retval->hitlist_max = hitlist_max;
    retval->hitlist_array = (HbnHitList*)calloc(hitlist_max, sizeof(HbnHitList));
    retval->hsp_alloc = SmallObjectAllocNew(sizeof(HbnHSP));
    retval->hsplist_alloc = SmallObjectAllocNew(sizeof(HbnHSPList));
    for (int i = 0; i < hitlist_max; ++i) {
        retval->hitlist_array[i].hsplist_array = NULL;
        retval->hitlist_array[i].hsplist_count = 0;
    }
    ks_init(retval->aligned_strings);
    kv_init(retval->cns_hit_list);
    ks_init(retval->output_buf);
    return retval;    
}

HbnAlignResults*
HbnAlignResultsFree(HbnAlignResults* results)
{
    free(results->hitlist_array);
    SmallObjectAllocFree(results->hsp_alloc);
    SmallObjectAllocFree(results->hsplist_alloc);
    ks_destroy(results->aligned_strings);
    kv_destroy(results->cns_hit_list);
    ks_destroy(results->output_buf);
    free(results);
    return NULL;
}

void
HbnAlignResultsClear(HbnAlignResults* results, int num_queries)
{
    for (int i = 0; i < results->num_queries; ++i) results->hitlist_array[i].hsplist_count = 0;
    hbn_assert(num_queries <= results->hitlist_max);
    results->num_queries = num_queries;
    SmallObjectAllocClear(results->hsp_alloc);
    SmallObjectAllocClear(results->hsplist_alloc);
    ks_clear(results->aligned_strings);
    kv_clear(results->cns_hit_list);
    ks_clear(results->output_buf);
}