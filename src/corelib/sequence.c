#include "sequence.h"

HbnSequence* 
hbn_sequence_new(const char* path)
{
    HbnSequence* seq = (HbnSequence*)malloc(sizeof(HbnSequence));
    ks_init(seq->sequence);
    ks_init(seq->name);
    seq->reader = path ? hbn_line_reader_new(path) : NULL;
    return seq;
}

HbnSequence*
hbn_sequence_free(HbnSequence* seq)
{
    ks_destroy(seq->sequence);
    ks_destroy(seq->name);
    if (seq->reader) hbn_line_reader_free(seq->reader);
    free(seq);
    return NULL;
}

void
hbn_sequence_reset(HbnSequence* seq)
{
    ks_clear(seq->sequence);
    ks_clear(seq->name);
}

static void
parse_def_line(HbnSequence* seq, kstring_t* line, const size_t line_number)
{
    hbn_assert(ks_A(*line, 0) == '>' || ks_A(*line, 0) == '@');
    size_t i = 1;
    for (; i < ks_size(*line); ++i) {
        if (isspace(ks_A(*line, i))) break;
    }
    if (i == 1) HBN_ERR("FastaReader: Defline lacks a proper ID around line %lu", line_number);
    kputsn(ks_s(*line) + 1, i - 1, &seq->name);
}

#define ascii_is_alpha(c) (((c) >= 'a' && (c) <= 'z') || ((c) >= 'A' && (c) <= 'Z'))
#define ascii_is_ambig_nuc(c) (nst_nt16_table[c] > 3 && nst_nt16_table[c] < 16)

static void
check_data_line(HbnSequence* seq, kstring_t* line)
{
    // make sure the first data line has at least SOME resemblance to actual sequence data
    if (!ks_empty(seq->sequence)) return;
    size_t good = 0, bad = 0;
    // in case the data has huge sequences all on the first line we do need
    // a cutoff and "70" seems reasonable since it's the default width of
    // CFastaOstream (as of 2017-03-09)
    size_t len_to_check = hbn_min(ks_size(*line), 70);
    size_t ambig_nuc = 0;
    for (size_t pos = 0; pos < len_to_check; ++pos) {
        unsigned char c = ks_A(*line, pos);
        if (ascii_is_alpha(c) || c == '*') {
            ++good;
            if (ascii_is_ambig_nuc(c)) ++ambig_nuc;
        } else if (c == '-') {
            // hyphens
        } else if (isspace(c) || (c >= '0' && c <= '9')) {
            // treat whitespace and digits as neutral
        } else if (c == ';') {
            break; // comment --- ignore rest of line
        } else {
            ++bad;
        }
    }

    if (bad >= good / 3 &&
        (len_to_check > 3 || good == 0 || bad > good)) {
        HBN_ERR("FastaReader: Near line %lu, there's a line that doesn't look like plausible data, but it's not marked as defline or comment", seq->reader->line_number);
    }
    // warn if more than a certain percentage is ambiguous nucleotides
    const static size_t kWarnPercentAmbiguous = 40; // e.g. "40" means "40%"
    const size_t percent_ambig = (good == 0) ? 100 : (ambig_nuc * 100 / good);
    if (len_to_check > 3 && percent_ambig > kWarnPercentAmbiguous) {
        char* msg = (char*)malloc(ks_size(seq->name) + 200);
        kputc('\0', &seq->name);
        sprintf(msg, "FastaReader: Start of first data line in sequence with header '%s' is about %lu%% ambiguous nucleotides"
                     " (shouldn't be over %lu%%)", ks_s(seq->name),  percent_ambig, kWarnPercentAmbiguous);
        HBN_WARN("%s", msg);
        free(msg);
        ks_pop_back(seq->name);
    }
}

static void
parse_data_line(HbnSequence* seq, kstring_t* line)
{
    check_data_line(seq, line);
    // most lines won't have a comment (';') 
    // so optimize for that case as much as possible
    const size_t s_len = ks_size(*line);
    ks_reserve(&seq->sequence, ks_size(seq->sequence) + s_len);
    size_t seq_pos = ks_size(seq->sequence);
    for (size_t pos = 0; pos < s_len; ++pos) {
        const unsigned char c = ks_A(*line, pos);
        if (c == ';') break;
        if (c == '\t' || c == '\v' || c == '\f' || c == '\r' || c == ' ') continue;
        if (nst_nt16_table[c] == 16) {
            char msg[512];
            sprintf(msg, "FastaReader: There are invalid nucleotide residue(s) in input sequence around line %lu and column %lu", seq->reader->line_number, pos);
            HBN_ERR("%s", msg);
        }
        ks_A(seq->sequence, seq_pos) = c;
        ++seq_pos;
    }
    ks_set_size(&seq->sequence, seq_pos);
    //HBN_LOG("seq_pos = %lu, s_len = %lu", seq_pos, s_len);
}

static void
parse_fastq_sequence(HbnSequence* seq)
{
    if (hbn_line_reader_at_eof(seq->reader)) return;

    // header
    hbn_line_reader_read_one_line(seq->reader);
    parse_def_line(seq, &seq->reader->line, seq->reader->line_number);
    // sequence
    if (hbn_line_reader_at_eof(seq->reader)) {
        HBN_LOG("FastaReader: Unexpected end-of-file around line %lu", seq->reader->line_number);
    }
    hbn_line_reader_read_one_line(seq->reader);
    parse_data_line(seq, &seq->reader->line);
    // plus
    if (hbn_line_reader_at_eof(seq->reader)) {
        HBN_LOG("FastaReader: Unexpected end-of-file around line %lu", seq->reader->line_number);
    }
    hbn_line_reader_read_one_line(seq->reader);
    if (ks_empty(seq->reader->line) || ks_A(seq->reader->line, 0) != '+') {
        HBN_ERR("FastaReader: This line (around line %lu) seems not a proper third line of a fastq format sequence, which should starts with '+'", seq->reader->line_number);
    }
    // quality
    if (hbn_line_reader_at_eof(seq->reader)) {
        HBN_LOG("FastaReader: Unexpected end-of-file around line %lu", seq->reader->line_number);
    }
    hbn_line_reader_read_one_line(seq->reader);
    if (ks_size(seq->reader->line) != ks_size(seq->sequence)) {
        HBN_LOG("%lu -- %lu", ks_size(seq->reader->line), ks_size(seq->sequence));
        HBN_ERR("FastaReader: Quality score line and the sequence line do not have the same length around line %lu", seq->reader->line_number);
    }
}

void hbn_sequence_read_one_seq(HbnSequence* seq)
{
    hbn_sequence_reset(seq);
    int need_defline = 1;
    while (!hbn_line_reader_at_eof(seq->reader)) {
        char c = hbn_line_reader_peek_char(seq->reader);
        if (hbn_line_reader_at_eof(seq->reader)) {
            HBN_LOG("FastaReader: Unexpected end-of-file around line %lu", seq->reader->line_number);
        }
        if (c == '>') {
            hbn_line_reader_read_one_line(seq->reader);
            if (need_defline) {
                parse_def_line(seq, &seq->reader->line, seq->reader->line_number);
                need_defline = 0;
                continue;
            } else {
                // start of the next sequence
                hbn_line_reader_unget_line(seq->reader);
                break;
            }
        }
        if (c == '@') {
            // a fastq file, process seperately
            parse_fastq_sequence(seq);
            need_defline = 0;
            break;
        }

        hbn_line_reader_read_one_line(seq->reader);
        kstring_t* line = &seq->reader->line;
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue; // ignore lines containing only whitespace
        c = ks_front(*line);
        if (c == '!' || c == '#' || c == ';') {
            // no content, just a comment or blank line
            continue;
        } else if (need_defline) {
            if (c == '>') {
                parse_def_line(seq, line, seq->reader->line_number);
                need_defline = 0;
                continue;
            } else if (c == '@') {
                hbn_line_reader_unget_line(seq->reader);
                parse_fastq_sequence(seq);
                need_defline = 0;
                break;
            } else {
                HBN_ERR("FastaReader: Input doesn't start with a defline or comment around line %lu", seq->reader->line_number);
            }
        }

        parse_data_line(seq, line);
    }

    if (need_defline && hbn_line_reader_at_eof(seq->reader)) {
        HBN_ERR("FastaReader: Expected defline around line %lu", seq->reader->line_number);
    }
}