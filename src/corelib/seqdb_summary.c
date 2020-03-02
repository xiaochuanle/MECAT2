#include "seqdb_summary.h"

#include "db_format.h"
#include "line_reader.h"
#include "ksort.h"
#include "fasta.h"

#include <ctype.h>

#define size_t_gt(a, b) ((a) > (b))
KSORT_INIT(size_t_gt, size_t, size_t_gt);

typedef struct {
    double percentile;
    const char* Nx;
    const char* Lx;
} DbLenPercentile;

static const DbLenPercentile LenStatPoints[DbLenStatN] = {
    { 0.1,  "N10", "L10" },
    { 0.2,  "N20", "L20" },
    { 0.25, "N25", "L25" },
    { 0.3,  "N30", "L30" },
    { 0.4,  "N40", "L40" },
    { 0.5,  "N50", "L50" },
    { 0.6,  "N60", "L60" },
    { 0.7,  "N70", "L70" },
    { 0.75, "N75", "L75" },
    { 0.8,  "N80", "L80" },
    { 0.9,  "N90", "L90" }
};

static void
calc_percentile_info(const size_t* len_array, const int len_count, const size_t target, int* nx, int* lx)
{
    size_t sum = 0;
    for (int i = 0; i < len_count; ++i) {
        sum += len_array[i];
        if (sum >= target) {
            *nx = len_array[i];
            *lx = i + 1;
            return;
        }
    }
    *nx = len_array[len_count - 1];
    *lx = len_count;
}

static void
fill_summary_len_info(size_t* len_array, const int len_count, CSeqDBSummary* summary)
{
    if (len_count == 0) return;
    ks_introsort_size_t_gt(len_count, len_array);
    size_t sum = 0;
    for (int i = 0; i < len_count; ++i) sum += len_array[i];
    summary->seq_count = len_count;
    summary->res_count = sum;
    summary->max_len = len_array[0];
    summary->min_len = len_array[len_count-1];
    summary->avg_len = sum / len_count;
    summary->med_len = len_array[len_count/2];
    for (int i = 0; i < DbLenStatN; ++i) {
        size_t target = sum * LenStatPoints[i].percentile;
        calc_percentile_info(len_array, len_count, target, &summary->len_stats[i].first, &summary->len_stats[i].second);
    }
}

static CSeqDBSummary*
load_summary_from_seqdb(const char* seqdb_dir, const char* seqdb_title, const int min_seq_size)
{
    CSeqDBSummary* summary = (CSeqDBSummary*)calloc(1, sizeof(CSeqDBSummary));
    CSeqDBInfo dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, 0);
    if (dbinfo.num_seqs == 0) return NULL;
    CSeqInfo* seqinfo_array = load_seq_infos(seqdb_dir, seqdb_title, 0, dbinfo.num_seqs);
    size_t* len_array = (size_t*)malloc(sizeof(size_t) * dbinfo.num_seqs);
    int n = 0;
    for (int i = 0; i < dbinfo.num_seqs; ++i) {
        if (seqinfo_array[i].seq_size >= min_seq_size) len_array[n++] = seqinfo_array[i].seq_size;
    }
    free(seqinfo_array);
    CAmbigSubseq* ambig_array = load_ambig_subseqs(seqdb_dir, seqdb_title, dbinfo.ambig_offset_from, dbinfo.ambig_offset_to);
    summary->ambig_subseq_count = dbinfo.ambig_offset_to - dbinfo.ambig_offset_from;
    summary->ambig_res_count = 0;
    for (size_t i = dbinfo.ambig_offset_from; i < dbinfo.ambig_offset_to; ++i) {
        summary->ambig_res_count += ambig_array[i].count;
    }
    free(ambig_array);
    fill_summary_len_info(len_array, n, summary);
    free(len_array);

    return summary;
}

static void
add_seq_ambig_info(const char* seq, const size_t len, size_t* ambig_subseq_count, size_t* ambig_res_count)
{
    size_t subseqs = 0;
    size_t res = 0;
    size_t i = 0;
    while (i < len) {
        int c = seq[i];
        c = nst_nt16_table[c];
        if (c > 3) {
            ++subseqs;
            size_t j = i;
            while (j < len) {
                c = seq[j];
                c = nst_nt16_table[c];
                if (c < 4) break;
                ++res;
                ++j;
            }
            i = j;
        } else {
            ++i;
        }
    }
    *ambig_subseq_count += subseqs;
    *ambig_res_count += res;
}

static void
add_len_info_from_one_file(const char* file_path, const int min_seq_size, vec_size_t* len_list, size_t* ambig_subseq_count, size_t* ambig_res_count)
{
    HbnFastaReader* reader = HbnFastaReaderNew(file_path);
    HbnFastaReaderSkipErrorFormatedSequences(reader);
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        if (!HbnFastaReaderReadOneSeq(reader)) continue;
        if (ks_size(reader->sequence) < min_seq_size) continue;
        kv_push(size_t, *len_list, ks_size(reader->sequence));
        add_seq_ambig_info(ks_s(reader->sequence), ks_size(reader->sequence), ambig_subseq_count, ambig_res_count);
    }
    HbnFastaReaderFree(reader);
}

CSeqDBSummary*
CSeqDBSummaryBuild(int argc, char* argv[])
{
    if (argc < 2) return NULL;
    
    int min_seq_size = 0;
    if (argc > 2) {
        if (argc > 4) {
            HBN_LOG("too many arguments");
            return NULL;
        }
        if (!string_is_valid_number(argv[2])) {
            const char* seqdb_dir = argv[1];
            const char* seqdb_title = argv[2];
            if (argc == 4) {
                if (!string_is_valid_number(argv[3])) {
                    HBN_LOG("'%s' seems not a plausible number");
                    return NULL;
                }
                min_seq_size = atoi(argv[3]);
            }
            return load_summary_from_seqdb(seqdb_dir, seqdb_title, min_seq_size);
        }
        min_seq_size = atoi(argv[argc-1]);
    }

    kv_dinit(vec_size_t, len_list);
    size_t ambig_subseq_count = 0;
    size_t ambig_res_count = 0;
    const char* input = argv[1];
    EDbFormat fmt = hbn_guess_db_format(input);
    if (fmt == eDbFormatUnknown) {
        HbnLineReader* reader = HbnLineReaderNew(input);
        while (!HbnLineReaderAtEof(reader)) {
            HbnLineReaderReadOneLine(reader);
            kstring_t* line = &reader->line;
            kputc('\0', line);
            add_len_info_from_one_file(ks_s(*line), min_seq_size, &len_list, &ambig_subseq_count, &ambig_res_count);
        }
        HbnLineReaderFree(reader);
    } else {
        add_len_info_from_one_file(input, min_seq_size, &len_list, &ambig_subseq_count, &ambig_res_count);
    }  

    CSeqDBSummary* summary = (CSeqDBSummary*)calloc(1, sizeof(CSeqDBSummary));
    summary->ambig_subseq_count = ambig_subseq_count;
    summary->ambig_res_count = ambig_res_count;
    fill_summary_len_info(kv_data(len_list), kv_size(len_list), summary);
    kv_destroy(len_list);
    return summary;
}

void
CSeqDBSummaryView(const CSeqDBSummary* summary)
{
    const int kPrologLen = 24;
    FILE* out = stderr;
    print_fixed_width_string(out, "sequences:", kPrologLen);
    print_digit_with_comma(out, summary->seq_count);
    fprintf(out, "\n");
    print_fixed_width_string(out, "residues:", kPrologLen);
    print_digit_with_comma(out, summary->res_count);
    fprintf(out, "\n");
    if (summary->ambig_subseq_count) {
        print_fixed_width_string(out, "ambig subseqs:", kPrologLen);
        print_digit_with_comma(out, summary->ambig_subseq_count);
        fprintf(out, "\n");
        print_fixed_width_string(out, "ambig residues:", kPrologLen);
        print_digit_with_comma(out, summary->ambig_res_count);
        double res_perc = 100.0 * summary->ambig_res_count / summary->res_count;
        fprintf(out, "(%.2lf%%)", res_perc);
        fprintf(out, "\n");
    }
    print_fixed_width_string(out, "max:", kPrologLen);
    print_digit_with_comma(out, summary->max_len);
    fprintf(out, "\n");
    print_fixed_width_string(out, "min:", kPrologLen);
    print_digit_with_comma(out, summary->min_len);
    fprintf(out, "\n");
    print_fixed_width_string(out, "avg:", kPrologLen);
    print_digit_with_comma(out, summary->avg_len);
    fprintf(out, "\n");
    print_fixed_width_string(out, "median:", kPrologLen);
    print_digit_with_comma(out, summary->med_len);
    fprintf(out, "\n");
    fprintf(out, "\n");
    for (int i = 0; i < DbLenStatN; ++i) {
        char buf[64], buf1[64];
        sprintf(buf, "(%s, %s):", LenStatPoints[i].Nx, LenStatPoints[i].Lx);
        print_fixed_width_string(out, buf, kPrologLen);
        i64_to_string_with_comma_r(summary->len_stats[i].first, buf);
        print_fixed_width_string(out, buf, 10);
        i64_to_string_with_comma_r(summary->len_stats[i].second, buf1);
        print_fixed_width_string(out, buf1, 10);
        fprintf(out, "\n");
        //fprintf(out, "(%s, %s)\n", buf, buf1);
    }
}

CSeqDBSummary*
CSeqDBSummaryFree(CSeqDBSummary* summary)
{
    free(summary);
    return NULL;
}