#include "hbn_format.h"

#include "sam_writter.h"

const char* OutFmtStrings[] = {
    "can", "canbin", "m4", "m4bin", "m4x", "m4y", "sam", "UnkwownFmt"
};

const char*
outfmt_2_string(EOutputFormat outfmt)
{
    return OutFmtStrings[outfmt];
}

EOutputFormat
string_2_outfmt(const char* str)
{
    int id = 7;
    if (strcmp(str, OutFmtStrings[0]) == 0) id = 0;
    if (strcmp(str, OutFmtStrings[1]) == 0) id = 1;
    if (strcmp(str, OutFmtStrings[2]) == 0) id = 2;
    if (strcmp(str, OutFmtStrings[3]) == 0) id = 3;
    if (strcmp(str, OutFmtStrings[4]) == 0) id = 4;
    if (strcmp(str, OutFmtStrings[5]) == 0) id = 5;
    if (strcmp(str, OutFmtStrings[6]) == 0) id = 6;
    return id;
}

void
dump_cns_hits(HbnConsensusInitHit* cns_hit_array,
    const int cns_hit_count,
    const int query_start_id,
    const int subject_start_id,
    const EOutputFormat outfmt,
    kstring_t* out_buf)
{
    for (int i = 0; i < cns_hit_count; ++i) {
        cns_hit_array[i].qid += query_start_id;
        cns_hit_array[i].sid += subject_start_id;
    }

    if (outfmt == eOutFmtCanBin) {
        const char* input = (const char*)(cns_hit_array);
        const int len = sizeof(HbnConsensusInitHit) * cns_hit_count;
        kputsn(input, len, out_buf);
    } else {
        hbn_assert(outfmt == eOutFmtCan);
        for (int i = 0; i < cns_hit_count; ++i) dump_cns_hit(ksprintf, out_buf, cns_hit_array[i]);
    }
}

void
dump_hbn_hsps(const text_t* queries,
    const text_t* database,
    HbnHSP* hsp_array,
    const int hsp_count,
    kstring_t* aligned_strings,
    const EOutputFormat outfmt,
    kstring_t* out)
{
    ks_introsort_hbn_hsp_score_gt(hsp_count, hsp_array);
    for (int i = 0; i < hsp_count; ++i) {
        HbnHSP* hsp = hsp_array + i;
        normalise_hbn_hsp_offsets(hsp);
        if (outfmt_is_sam(outfmt)) {
            print_sam_result(hsp,
                queries,
                database,
                aligned_strings,
                out);
            continue;
        }

        if (outfmt == eOutFmtM4Bin) {
            hsp->qid += queries->dbinfo.seq_start_id;
            if (hsp->sid == -1) hsp->sid = hsp->qid;
            hsp->sid += database->dbinfo.seq_start_id;
            kputsn((const char*)(hsp), sizeof(HbnHSP), out);
            continue;
        } 

        hbn_assert(hsp->qid < seqdb_num_seqs(queries));
        hbn_assert(hsp->sid < seqdb_num_seqs(database));
        const char* qhdr = seqdb_seq_name(queries, hsp->qid);
        const char* shdr = (hsp->sid == -1) ? qhdr : seqdb_seq_name(database, hsp->sid);
        dump_hbn_hsp_gi(ksprintf, out, *hsp);
        if (outfmt == eOutFmtM4x) {
            ks_back(*out) = '\t';
            print_cigar(hsp, aligned_strings, 0, out);
            kputc('\n', out);
        } else if (outfmt == eOutFmtM4y) {
            int n = hsp->saln_offset - hsp->qaln_offset;
            const char* q = ks_s(*aligned_strings) + hsp->qaln_offset;
            const char* s = ks_s(*aligned_strings) + hsp->saln_offset;
            kputsn(q, n, out);
            kputc('\n', out);
            for (int u = 0; u < n; ++u) {
                if (q[u] == s[u]) kputc('|', out);
                else kputc('*', out);
            }
            kputc('\n', out);
            kputsn(s, n, out);
            kputc('\n', out);
        }        
    }
}

void
dump_hit_list(const text_t* queries, 
    const text_t* database, 
    HbnHitList* hitlist, 
    kstring_t* aligned_strings,
    const EOutputFormat outfmt,
    kstring_t* out)
{
    ks_introsort_hbn_hsp_list_score_gt(hitlist->hsplist_count, hitlist->hsplist_array);
    for (int i = 0; i < hitlist->hsplist_count; ++i) {
        dump_hbn_hsps(queries,
            database,
            hitlist->hsplist_array[i].hsp_array,
            hitlist->hsplist_array[i].hsp_count,
            aligned_strings,
            outfmt,
            out);     
    }
}

void
dump_hbn_align_results(const text_t* queries,
    const text_t* database,
    HbnAlignResults* results,
    EOutputFormat outfmt)
{
    ks_clear(results->output_buf);
    if (outfmt_is_can(outfmt)) {
        dump_cns_hits(kv_data(results->cns_hit_list),
            kv_size(results->cns_hit_list),
            queries->dbinfo.seq_start_id,
            database->dbinfo.seq_start_id,
            outfmt,
            &results->output_buf);
        return;
    }

    for (int i = 0; i < results->num_queries; ++i) {
        if (!results->hitlist_array[i].hsplist_count) continue;
        dump_hit_list(queries,
            database,
            &results->hitlist_array[i],
            &results->aligned_strings,
            outfmt,
            &results->output_buf);
    }
}