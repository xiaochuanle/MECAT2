#include "hbn_align_one_volume.h"

#include "hbn_find_subseq_hit.h"
#include "hbn_extend_subseq_hit.h"
#include "mecat_results.h"
#include "../../corelib/m4_record.h"
#include "../../ncbi_blast/setup/hsp2string.h"

#include <pthread.h>

static hbn_task_struct* g_task_struct = NULL;
static int g_thread_index = 0;
static pthread_mutex_t g_thread_index_lock;
static int g_query_index = 0;
static pthread_mutex_t g_query_index_lock;

static void
init_global_values(hbn_task_struct* task_struct)
{
    g_task_struct = task_struct;
    g_thread_index = 0;
    pthread_mutex_init(&g_thread_index_lock, NULL);
    g_query_index = 0;
    pthread_mutex_init(&g_query_index_lock, NULL);
}

static int
get_next_query_chunk(
    const text_t* queries, 
    int* next_query_id,
    pthread_mutex_t* query_id_lock,
    BLAST_SequenceBlk* query_blk, 
    BlastQueryInfo* query_info)
{
    int from = 0, to = 0;
    pthread_mutex_lock(query_id_lock);
    from = *next_query_id;
    *next_query_id += HBN_QUERY_CHUNK_SIZE;
    pthread_mutex_unlock(query_id_lock);
    if (from >= queries->dbinfo.num_seqs) return 0;
    to = hbn_min(from + HBN_QUERY_CHUNK_SIZE, queries->dbinfo.num_seqs);
    int num_queries = to - from;

    int length = 0;
    for (int i = from; i < to; ++i) length += seqdb_seq_size(queries, i);
    length *= 2;

    BlastContextInfo ctx_info;
    int ctx_idx = 0;
    int seq_idx = 0;
    int max_length = 0;
    int min_length = I32_MAX;
    query_blk->sequence = (Uint1*)realloc(query_blk->sequence, length);
    query_blk->length = length;

    for (int i = from; i < to; ++i) {
        ctx_info.query_length = seqdb_seq_size(queries, i);
        ctx_info.eff_searchsp = 0;
        ctx_info.length_adjustment = 0;
        ctx_info.query_index = i;
        ctx_info.frame = 0;
        ctx_info.is_valid = TRUE;
        ctx_info.segment_flags = 0;

        max_length = hbn_max(max_length, ctx_info.query_length);
        min_length = hbn_min(min_length, ctx_info.query_length);

        size_t seq_from = seqdb_seq_offset(queries, i);
        size_t seq_to = seq_from + seqdb_seq_size(queries, i);
        ctx_info.query_offset = seq_idx;
        query_info->contexts[ctx_idx++] = ctx_info;
        for (size_t k = seq_from; k < seq_to; ++k) {
            Uint1 c = _get_pac(queries->packed_seq, k);
            query_blk->sequence[seq_idx++] = c;
        }

        ctx_info.query_offset = seq_idx;
        query_info->contexts[ctx_idx++] = ctx_info;
        size_t k = seq_to;
        while (k > seq_from) {
            --k;
            Uint1 c = _get_pac(queries->packed_seq, k);
            c = c ^ 3;
            query_blk->sequence[seq_idx++] = c;
        }
    }
    hbn_assert(seq_idx == length);

    query_info->first_context = 0;
    query_info->last_context = ctx_idx - 1;
    query_info->num_queries = num_queries;
    query_info->max_length = max_length;
    query_info->min_length = min_length;  

    return num_queries;  
}

static void
dump_subseq_hits(HbnSubseqHit* subseq_array,
    const int subseq_count,
    const text_t* query_vol,
    const text_t* subject_vol,
    EOutputFormat outfmt,
    kstring_t* out_buf)
{
    ks_clear(*out_buf);
    const int query_start_id = query_vol->dbinfo.seq_start_id;
    const int subject_start_id = subject_vol->dbinfo.seq_start_id;

    char line[1024];
    for (int i = 0; i < subseq_count; ++i) {
        HbnSubseqHit* subseq = subseq_array + i;
        if (outfmt == eSeqid || outfmt == eSeqidx) {
            HbnConsensusInitHit cns_hit;
            int query_length = seqdb_seq_size(query_vol, subseq->qid);
            cns_hit.qid = subseq->qid + query_start_id;
            cns_hit.qoff = (subseq->qdir == FWD) ? subseq->qoff : (query_length - 1 - subseq->qoff);
            cns_hit.sid = subseq->sid + subject_start_id;
            cns_hit.soff = subseq->soff;
            cns_hit.score = subseq->score;
            cns_hit.strand = (subseq->qdir == FWD);
            if (outfmt == eSeqid) {
                dump_cns_hit(sprintf, line, cns_hit);
                kputsn(line, strlen(line), out_buf);
            } else {
                const char* input = (const char*)(&cns_hit);
                const int input_len = sizeof(HbnConsensusInitHit);
                kputsn(input, input_len, out_buf);
            }
            continue;
        }
    }
}

static void
dump_m4_hits(const text_t* query_vol,
    const text_t* subject_vol,
    HbnHSPResults* results,
    const HbnProgramOptions* opts)
{
    EOutputFormat outfmt = opts->outfmt;
    ks_dinit(line);
    ks_clear(results->output_buf);
    for (int i = 0; i < results->num_queries; ++i) {
        BlastHitList* hit_list = results->hitlist_array + i;
        if (hit_list->hsplist_array == 0) continue;
        for (int j = 0; j < hit_list->hsplist_count; ++j) {
            BlastHSPList* hsp_list = hit_list->hsplist_array[j];
            for (int k = 0; k < hsp_list->hspcnt; ++k) {
                BlastHSP* hsp = hsp_list->hsp_array[k];
                if (hsp->hbn_query.strand == REV) {
                    int offset = hsp->hbn_query.seq_size - hsp->hbn_query.end;
                    int end = hsp->hbn_query.seq_size - hsp->hbn_query.offset;
                    hsp->hbn_query.offset = offset;
                    hsp->hbn_query.end = end;
                }
                if (hsp->hbn_subject.strand == REV) {
                    int offset = hsp->hbn_subject.seq_size - hsp->hbn_subject.end;
                    int end = hsp->hbn_subject.seq_size - hsp->hbn_subject.offset;
                    hsp->hbn_subject.offset = offset;
                    hsp->hbn_subject.end = end;
                }
                const char* qname = seqdb_seq_name(query_vol, hsp->hbn_query.oid);
                const char* sname = seqdb_seq_name(subject_vol, hsp->hbn_subject.oid);
                hsp->hbn_query.oid = query_vol->dbinfo.seq_start_id + hsp->hbn_query.oid;
                hsp->hbn_subject.oid = subject_vol->dbinfo.seq_start_id + hsp->hbn_subject.oid;

                if (outfmt >= eM4 && outfmt <= eM4x) {
                    print_one_m4_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_cigar,
                        opts->dump_md,
                        outfmt == eM4x,
                        &line,
                        &results->output_buf);
                } else if (outfmt == ePaf) {
                    print_one_paf_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_cigar,
                        opts->dump_md,
                        &results->output_buf);
                } else if (outfmt == eSAM) {
                    print_one_sam_result(hsp,
                        &results->aligned_strings,
                        qname,
                        sname,
                        opts->dump_md,
                        &results->output_buf);
                }
            }
        }
    }
    ks_destroy(line);
}

static void
align_one_block(hbn_task_struct* task_struct, const int thread_id, BLAST_SequenceBlk* query_blk, BlastQueryInfo* query_info)
{
    hbn_assert(thread_id < task_struct->opts->num_threads);
    CSeqDB* query_vol = task_struct->query_vol;
    CSeqDB* subject_vol = task_struct->subject_vol;
    WordFindData* word_data = task_struct->word_data_array[thread_id];
    HbnHSPResults* results = task_struct->results_array[thread_id];
    vec_init_hit* init_hit_list = &word_data->init_hit_list;
    const HbnProgramOptions* opts = task_struct->opts;
    HbnSubseqHitExtnData* extn_data = task_struct->hit_extn_data_array[thread_id];
    HbnHSPResultsClear(results, query_info->num_queries);
    vec_int_pair* seeding_subseq_list = &word_data->seeding_subseqs;
    kv_resize(IntPair, *seeding_subseq_list, 1);
    vec_subseq_hit* fwd_sbjct_subseq_list = &extn_data->fwd_sbjct_subseq_list;
    vec_subseq_hit* rev_sbjct_subseq_list = &extn_data->rev_sbjct_subseq_list;
    vec_subseq_hit* sbjct_subseq_list = &extn_data->sbjct_subseq_list;
    kv_dinit(vec_subseq_hit, subseq_hit_sink);
    
    for (int i = 0; i < query_info->num_queries; ++i) {
        kv_clear(*init_hit_list);
        int ctx_id = i << 1;
        BlastContextInfo ctx_info = query_info->contexts[ctx_id];
        int query_id = ctx_info.query_index;
        int query_strand = FWD;
        int query_length = ctx_info.query_length;
        u8* fwd_query = query_blk->sequence + ctx_info.query_offset;
        IntPair ip = { 0, query_length };
        kv_front(*seeding_subseq_list) = ip;
        ddfs_find_candidates(word_data, fwd_query, query_id, query_vol->dbinfo.seq_start_id, query_strand, query_length);

        ctx_id += 1;
        ctx_info = query_info->contexts[ctx_id];
        query_strand = REV;
        u8* rev_query = query_blk->sequence + ctx_info.query_offset;
        ddfs_find_candidates(word_data, rev_query, query_id, query_vol->dbinfo.seq_start_id, query_strand, query_length);

        HbnInitHit* init_hit_array = kv_data(*init_hit_list);
        int init_hit_count = kv_size(*init_hit_list);
        //if (query_id == 4924) HBN_LOG("number of init hits: %d", init_hit_count);

        find_candidate_subject_subseqs(init_hit_array,
            init_hit_count,
            query_id,
            query_length,
            opts,
            subject_vol,
            fwd_sbjct_subseq_list,
            rev_sbjct_subseq_list,
            sbjct_subseq_list);

        if (opts->outfmt >= eSeqid && opts->outfmt <= eSubseqx) {
            kv_push_v(HbnSubseqHit, 
                subseq_hit_sink, 
                kv_data(*sbjct_subseq_list), 
                kv_size(*sbjct_subseq_list));
            continue;
        }

        hbn_extend_query_subseq_hit_list(kv_data(*sbjct_subseq_list),
            kv_size(*sbjct_subseq_list),
            fwd_query,
            rev_query,
            query_id,
            query_length,
            subject_vol,
            opts,
            extn_data,
            results->hitlist_array + i,
            results);
    }

    if (opts->outfmt >= eSeqid && opts->outfmt <= eSubseqx) {
        dump_subseq_hits(kv_data(subseq_hit_sink),
            kv_size(subseq_hit_sink),
            query_vol,
            subject_vol,
            opts->outfmt,
            &results->output_buf);
    } else if (opts->outfmt >= eM4 && opts->outfmt <= eSAM) {
        dump_m4_hits(query_vol, subject_vol, results, opts);
    }

    pthread_mutex_lock(&task_struct->out_lock);
    if (task_struct->out) {
        hbn_fwrite(ks_s(results->output_buf), 1, ks_size(results->output_buf), task_struct->out);
    }
    if (task_struct->qi_vs_sj_out) {
        hbn_fwrite(ks_s(results->output_buf), 1, ks_size(results->output_buf), task_struct->qi_vs_sj_out);
    }
    pthread_mutex_unlock(&task_struct->out_lock);

    kv_destroy(subseq_hit_sink);
    int gid = query_info->contexts[0].query_index;
    if (gid && (gid % 1000 == 0)) HBN_LOG("%8d queries processed", gid);
}

static void*
hbn_align_worker()
{
    int thread_id = 0;
    pthread_mutex_lock(&g_thread_index_lock);
    thread_id = g_thread_index++;
    pthread_mutex_unlock(&g_thread_index_lock);
    hbn_assert(thread_id < g_task_struct->opts->num_threads);
    BLAST_SequenceBlk* query_blk = BLAST_SequenceBlkNew();
    BlastQueryInfo* query_info = BlastQueryInfoNew(HBN_QUERY_CHUNK_SIZE * 2);

    while (get_next_query_chunk(g_task_struct->query_vol,
                &g_query_index,
                &g_query_index_lock,
                query_blk,
                query_info)) {
        //HBN_LOG("mapping 20 reads");
        align_one_block(g_task_struct, thread_id, query_blk, query_info);
        //break;
    }
    BLAST_SequenceBlkFree(query_blk);
    BlastQueryInfoFree(query_info);
    return NULL;
}

void
hbn_align_one_volume(hbn_task_struct* task_struct)
{
    init_global_values(task_struct);
    const int num_threads = task_struct->opts->num_threads;
    pthread_t job_ids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(job_ids + i, NULL, hbn_align_worker, NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(job_ids[i], NULL);
    }
}