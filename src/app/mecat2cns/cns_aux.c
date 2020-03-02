#include "cns_aux.h"

#include "../../corelib/ksort.h"
#include "../../corelib/partition_aux.h"
#include "../../ncbi_blast/c_ncbi_blast_aux.h"

void 
normalize_gaps(const char* qstr, 
    const char* tstr, 
    const int aln_size, 
    kstring_t* qnorm, 
    kstring_t* tnorm, 
    const BOOL push)
{
    ks_clear(*qnorm);
    ks_clear(*tnorm);

#ifndef NDEBUG
    int qcnt = 0, tcnt = 0;
    for (int i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != GAP_CHAR) ++qcnt;
        if (tc != GAP_CHAR) ++tcnt;
    }
#endif

    // convert mismatches to indels
    for (int i = 0; i < aln_size; ++i)
    {
        const char qc = qstr[i];
        const char tc = tstr[i];
        if (qc != tc && qc != GAP_CHAR && tc != GAP_CHAR) {
            kputc(GAP_CHAR, qnorm);
            kputc(qc, qnorm);
            kputc(tc, tnorm);
            kputc(GAP_CHAR, tnorm);
        } else { 
            kputc(qc, qnorm);
            kputc(tc, tnorm);
        }
    }

    // push gaps to the right, but not pass the end
    if (push)
    {
        int qlen = ks_size(*qnorm);
        char* qn = ks_s(*qnorm);
        int tlen = ks_size(*tnorm);
        char* tn = ks_s(*tnorm);
        for (int i = 0; i < qlen - 1; ++i)
        {
            // push target gaps
            if (tn[i] == GAP_CHAR)
            {
                int j = i;
                while (1)
                {
                    const char c = tn[++j];
                    if (c != GAP_CHAR || j > qlen - 1)
                    {
                        if (c == qn[i]) { tn[i] = c; tn[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
            // push query gaps
            if (qn[i] == GAP_CHAR)
            {
                int j = i;
                while (1)
                {
                    const char c = qn[++j];
                    if (c != GAP_CHAR || j > tlen - 1)
                    {
                        if (c == tn[i]) { qn[i] = c; qn[j] = GAP_CHAR; }
                        break;
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    int qcnt2 = 0, tcnt2 = 0;
    for (size_t i = 0; i < ks_size(*qnorm); ++i) 
        if (ks_A(*qnorm, i) != GAP_CHAR) ++qcnt2;
    for (size_t i = 0; i < ks_size(*tnorm); ++i)
        if (ks_A(*tnorm, i) != GAP_CHAR) ++tcnt2;
    hbn_assert(qcnt == qcnt2);
    hbn_assert(tcnt == tcnt2);
#endif
}

HbnConsensusInitHit*
load_and_sort_cns_hits(const char* can_dir, 
    const int pid, 
    const BOOL use_batch_mode,
    size_t* hit_count)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    HbnConsensusInitHit* hit_array = load_part_records(path, sizeof(HbnConsensusInitHit), hit_count);
    size_t n = *hit_count;
    if (n == 0) return NULL;
    ks_introsort_cns_hit_sid_lt(n, hit_array);
    if (!use_batch_mode) return hit_array;

    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && hit_array[i].sid == hit_array[j].sid) ++j;
        int cnt = j - i;
        if (cnt > MAX_CNS_OVLPS) {
            for (size_t k = i + MAX_CNS_OVLPS; k < j; ++k) hit_array[k].sid = -1;
        }
        i = j;
    }

    i = 0;
    for (size_t k = 0; k < n; ++k) if (hit_array[k].sid >= 0) hit_array[i++] = hit_array[k];
    *hit_count = i;
    return hit_array;
}

BOOL
set_next_raw_read_batch_info(RawReadCnsInfo* cns_info_array,
    int* cns_info_count,
    HbnConsensusInitHit* cns_hit_array,
    size_t cns_hit_count,
    size_t* cns_hit_idx,
    RawReadsReader* raw_reads,
    const int batch_size)
{
    size_t i = *cns_hit_idx;
    if (i >= cns_hit_count) return FALSE;
    int p = 0;
    memset(cns_info_array, 0, sizeof(RawReadCnsInfo) * batch_size);
    int hit_cnt = 0;
    while (i < cns_hit_count) {
        size_t j = i + 1;
        while (j < cns_hit_count && cns_hit_array[i].sid == cns_hit_array[j].sid) ++j;
        cns_info_array[p].oid = cns_hit_array[i].sid;
        cns_info_array[p].can_from = i;
        cns_info_array[p].can_to = j;
        hit_cnt += j - i;
        i = j;
        ++p;
        if (p == batch_size) break;
    }

    *cns_info_count = p;
    *cns_hit_idx = i;
    RawReadsReaderLoadRawReads(cns_hit_array + cns_info_array[0].can_from,
        hit_cnt,
        raw_reads);
    return TRUE;
}

CnsThreadData*
CnsThreadDataNew(const HbnProgramOptions* opts,
    RawReadsReader* raw_reads,
    RawReadCnsInfo* cns_info_array,
    int* cns_info_idx,
    pthread_mutex_t* cns_info_idx_lock,
    kstring_t* cns_out,
    pthread_mutex_t* out_lock)
{
    CnsThreadData* data = (CnsThreadData*)calloc(1, sizeof(CnsThreadData));
    data->opts = opts;
    data->raw_reads = raw_reads;
    data->cns_info_array = cns_info_array;
    data->cns_info_idx = cns_info_idx;
    data->cns_info_idx_lock = cns_info_idx_lock;
    ks_init(data->qaux);
    ks_init(data->saux);
    kv_init(data->fwd_read);
    kv_init(data->rev_read);
    kv_init(data->fwd_subject);
    kv_init(data->rev_subject);
    data->cns_out = cns_out;
    data->cns_out_lock = out_lock;
    kv_init(data->cov_stats);
    data->cns_data = FCCnsDataNew();
    data->diff_data = DiffGapAlignDataNew();
    data->ksw = Ksw2DataNew();
    ksw2_extd2_set_params(data->ksw);
    return data;
}

CnsThreadData*
CnsThreadDataFree(CnsThreadData* data)
{
    ks_destroy(data->qaux);
    ks_destroy(data->saux);
    kv_destroy(data->fwd_read);
    kv_destroy(data->rev_read);
    kv_destroy(data->fwd_subject);
    kv_destroy(data->rev_subject);
    kv_destroy(data->cov_stats);
    FCCnsDataFree(data->cns_data);
    DiffGapAlignDataFree(data->diff_data);
    Ksw2DataFree(data->ksw);
    free(data);
    return NULL;
}

void
CnsThreadDataInit(CnsThreadData* data, 
    const u8* fwd_subject,
    const u8* rev_subject,
    const int subject_size)
{
    FCCnsDataClear(data->cns_data);
    kv_resize(u8, data->cov_stats, subject_size);
    kv_zero(u8, data->cov_stats);
}