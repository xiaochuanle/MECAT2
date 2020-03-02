#include "../../corelib/line_reader.h"
#include "../../corelib/khash.h"
#include "../../corelib/hbn_hit.h"
#include "../../corelib/fasta.h"
#include "../../corelib/seq_tag_report.h"
#include "../../corelib/string2hsp.h"
#include "../../corelib/ksort.h"
#include "../../corelib/cstr_util.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <ctype.h>

//>NW_003315967.2 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_1_CTG1_1
//>NT_187628.1 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_8_CTG1_1
//>NT_187627.1 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_6_CTG1_1
//>NW_003315968.2 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_2_CTG1_1
//>NW_003315969.2 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_3_CTG1_1
//>NW_003315970.2 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_4_CTG1_1
//>NT_187626.1 Homo sapiens chromosome 21 genomic scaffold, GRCh38.p7 alternate locus group ALT_REF_LOCI_1 HSCHR21_5_CTG2

int is_chr21(const char* name)
{
#if 0
    const char* chr21_names[7] = {
        "NW_003315967.2",
        "NT_187628.1",
        "NT_187627.1",
        "NW_003315968.2",
        "NW_003315969.2",
        "NW_003315970.2",
        "NT_187626.1"
    };
    for (int i = 0; i < 7; ++i) {
        if (strcmp(chr21_names[i], name) == 0) return 1;
    }
    return 0;
#endif

    const char* chr21_name = "NC_000021.9";
    return strcmp(chr21_name, name) == 0;
}

const int min_seq_size = 1000;
const int max_e = 50;
const double min_r = 0.6;
const double min_ident_perc = 80.0;
int seq_cnt = 0;
int chr21_seq_cnt = 0;
size_t chr21_seq_res = 0;

int is_perfect_hsp(HbnHSP* hsp)
{
    int r = (hsp->qoff <= max_e && hsp->qsize - hsp->qend <= max_e)
            ||
            (hsp->soff <= max_e && hsp->ssize - hsp->send <= max_e);
    return r;
}

int is_true_overlap_hsp(HbnHSP* hsp)
{
    int r = (hsp->soff <= max_e && hsp->qsize - hsp->qend <= max_e)
            ||
            (hsp->qoff <= max_e && hsp->ssize - hsp->send <= max_e);
    if (!r) return r;

    r = (hsp->qend - hsp->qoff >= hsp->qsize * min_r)
        ||
        (hsp->send - hsp->soff >= hsp->ssize * min_r);
    return r;
}

int examine_one_hsp_list(HbnHSP* hsp_array, int hsp_count, NameToIdMap* subject_name2id_map)
{
    if (hsp_array[0].qsize < min_seq_size) return 0;

    int n_perfect = 0;
    int n_ovlp = 0;
    HbnHSP hsp;
    for (int i = 0; i < hsp_count; ++i) {
        if (is_perfect_hsp(hsp_array + i)) {
            ++n_perfect;
            hsp = hsp_array[i];
        }
    }
    if (n_perfect > 1) return 0;
    if (n_perfect == 1) {
        if (hsp.ident_perc >= min_ident_perc) goto examine_name;
        return 0;
    }

    for (int i = 0; i < hsp_count; ++i) {
        if (is_true_overlap_hsp(hsp_array + i)) {
            ++n_ovlp;
            hsp = hsp_array[i];
        }
    }
    if (n_ovlp > 1) return 0;
    if (n_ovlp == 1) {
        if (hsp.ident_perc >= min_ident_perc) goto examine_name;
        return 0;
    }

#if 1
    for (int i = 0; i < hsp_count; ++i) {
        hsp = hsp_array[i];
        if (hsp.ident_perc < min_ident_perc) continue;
        int r = (hsp.send - hsp.soff >= hsp.ssize * min_r
                ||
                hsp.qend - hsp.qoff >= hsp.qsize * min_r);
        if (!r) continue;
        if (is_chr21(Name2IdMap_id2name(subject_name2id_map, hsp.sid))) return 1;
    }
#endif
    return 0;

examine_name:
    if (is_chr21(Name2IdMap_id2name(subject_name2id_map, hsp.sid))) {
        return 1;
    }
    return 0;
}

void dump_chr21_queries(const char* seqdb_dir_path, int* chr21_query_flag_array)
{
    int nvol = seqdb_load_num_volumes(seqdb_dir_path, INIT_QUERY_TITLE);
    int qid = -1;
    kv_dinit(vec_u8, query);
    hbn_dfopen(out, "chr21_queries.fasta", "w");
    int seq = 0;
    size_t res =0;
    for (int i = 0; i < nvol; ++i) {
        CSeqDB* vol = seqdb_load(seqdb_dir_path, INIT_QUERY_TITLE, i);
        for (int k = 0; k < seqdb_num_seqs(vol); ++k) {
            ++qid;
            if (!chr21_query_flag_array[qid]) continue;
            seqdb_extract_sequence(vol, k, FWD, &query);
            for (size_t pos = 0; pos < kv_size(query); ++pos) {
                int c = kv_A(query, pos);
                c = DECODE_RESIDUE(c);
                kv_A(query, pos) = c;
            }
            fprintf(out, ">%s\n", seqdb_seq_name(vol, k));
            hbn_fwrite(kv_data(query), 1, kv_size(query), out);
            fprintf(out, "\n");
            ++seq;
            res += kv_size(query);
        }
        CSeqDBFree(vol);
    }

    hbn_fclose(out);
    kv_destroy(query);
    HBN_LOG("dump %d queries, %zu residues", seq, res);
}

void dump_chr21(const char* seqdb_dir_path)
{
    int nvol = seqdb_load_num_volumes(seqdb_dir_path, INIT_DB_TITLE);
    kv_dinit(vec_u8, subject);
    hbn_dfopen(out, "chr21_subject.fasta", "w");
    int seq = 0;
    size_t res =0;
    for (int i = 0; i < nvol; ++i) {
        CSeqDB* vol = seqdb_load(seqdb_dir_path, INIT_DB_TITLE, i);
        for (int k = 0; k < seqdb_num_seqs(vol); ++k) {
            if (!is_chr21(seqdb_seq_name(vol, k))) continue;
            seqdb_extract_sequence(vol, k, FWD, &subject);
            for (size_t pos = 0; pos < kv_size(subject); ++pos) {
                int c = kv_A(subject, pos);
                c = DECODE_RESIDUE(c);
                kv_A(subject, pos) = c;
            }
            fprintf(out, ">%s\n", seqdb_seq_name(vol, k));
            hbn_fwrite(kv_data(subject), 1, kv_size(subject), out);
            fprintf(out, "\n");
            ++seq;
            res += kv_size(subject);
        }
        CSeqDBFree(vol);
    }

    hbn_fclose(out);
    kv_destroy(subject);
    HBN_LOG("dump %d queries, %zu residues", seq, res);    
}

#define seqinfo_size_gt(a, b) ((a).seq_size > (b).seq_size)
KSORT_INIT(seqinfo_size_gt, CSeqInfo, seqinfo_size_gt);

static void
print_seq_info(CSeqInfo* seqinfo_array, int num_seqs, const char* names)
{
    ks_introsort_seqinfo_size_gt(num_seqs, seqinfo_array);
    char size_buf[256];
    for (int i = 0; i < num_seqs; ++i) {
        char* name = names + seqinfo_array[i].hdr_offset;
        char* size = u64_to_string_datasize(seqinfo_array[i].seq_size, size_buf);
        fprintf(stderr, "%s    %s\n", name, size);
    }
}

#if 1
int main(int argc, char* argv[])
{
    hbn_assert(argc == 3);
    const char* seqdb_dir_path = argv[1];
    const char* hsp_path = argv[2];

    CSeqDBInfo db_info = seqdb_load_volume_info(seqdb_dir_path, INIT_DB_TITLE, 0);
    hbn_assert(db_info.hdr_offset_from == 0);
    char* subject_names = load_seq_headers(seqdb_dir_path, INIT_DB_TITLE, 0, db_info.hdr_offset_to);
    CSeqInfo* seqinfo_array = load_seq_infos(seqdb_dir_path, INIT_DB_TITLE, 0, db_info.num_seqs);
    print_seq_info(seqinfo_array, db_info.num_seqs, subject_names);
    return 0;

    CSeqDBInfo query_info = seqdb_load_volume_info(seqdb_dir_path, INIT_QUERY_TITLE, 0);
    hbn_assert(query_info.hdr_offset_from == 0);
    char* query_names = load_seq_headers(seqdb_dir_path, INIT_QUERY_TITLE, 0, query_info.hdr_offset_to);
    int* chr21_query_flag_array = (int*)calloc(query_info.num_seqs, sizeof(int));

    HbnHspReader* hspreader = HbnHspReaderNew(hsp_path, query_names, query_info.num_seqs, subject_names, db_info.num_seqs);
    kv_dinit(vec_hsp, hsp_list);
    HbnHSP hsp;

    while (1) {
        if (!HbnHspReaderGet(hspreader, &hsp)) break;
        kv_clear(hsp_list);
        kv_push(HbnHSP, hsp_list, hsp);
        int qid = hsp.qid;

        while (1) {
            if (!HbnHspReaderGet(hspreader, &hsp)) break;
            if (hsp.qid != qid) {
                HbnHspReaderUnget(hspreader);
                break;
            }
            kv_push(HbnHSP, hsp_list, hsp);
        }

        ++seq_cnt;
        HbnHSP* hsp_array = kv_data(hsp_list);
        int hsp_count = kv_size(hsp_list);
        for (int i = 0; i < hsp_count; ++i) {
            normalise_hbn_hsp_sdir(hsp_array + i);
            if (hsp_array[i].qdir == REV) {
                size_t qoff = hsp_array[i].qsize - hsp_array[i].qend;
                size_t qend = hsp_array[i].qsize - hsp_array[i].qoff;
                hsp_array[i].qoff = qoff;
                hsp_array[i].qend = qend;
            }
        }
        if (examine_one_hsp_list(hsp_array, hsp_count, hspreader->subject_name2id_map)) {
            ++chr21_seq_cnt;
            chr21_seq_res += hsp_array[0].qsize;
            chr21_query_flag_array[qid] = 1;
        }
    }
    HbnHspReaderFree(hspreader);
    free(query_names);
    free(subject_names);
    kv_destroy(hsp_list);

    HBN_LOG("total queries: %d, number of queries: %d, residues: %zu", seq_cnt, chr21_seq_cnt, chr21_seq_res);

    dump_chr21_queries(seqdb_dir_path, chr21_query_flag_array);
    dump_chr21(seqdb_dir_path);
}
#else

extern int find_chr21_name(const char* name);

int main(int argc, char* argv[])
{
    HbnFastaReader* reader = HbnFastaReaderNew(argv[1]);
    size_t total_res = 0;
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnFastaReaderReadOneSeq(reader);
        kputc('\0', &reader->name);
        const char* name = ks_s(reader->name);
        kputc('\0', &reader->comment);
        const char* comment = ks_s(reader->comment);
        //HBN_LOG("testing %s", name);
        if (find_chr21_name(comment)) {
            HBN_LOG("%s -- %zu", name, ks_size(reader->sequence));
            total_res += ks_size(reader->sequence);
        }
    }
    HBN_LOG("total res: %zu", total_res);
}
#endif