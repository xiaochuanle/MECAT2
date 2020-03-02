#include "fasta.h"

#include "khash.h"

#include <ctype.h>

#ifndef bool
typedef BOOL bool;
#define true TRUE
#define false FALSE
#endif

// The FASTA reader uses these heavily, but the standard versions
// aren't inlined on as many configurations as one might hope, and we
// don't necessarily want locale-dependent behavior anyway.

bool s_ASCII_IsUpper(unsigned char c)
{
    return c >= 'A'  &&  c <= 'Z';
}

bool s_ASCII_IsLower(unsigned char c)
{
    return c >= 'a'  &&  c <= 'z';
}

bool s_ASCII_IsAlpha(unsigned char c)
{
    return s_ASCII_IsUpper(c)  ||  s_ASCII_IsLower(c);
}

// The arg *must* be a lowercase letter or this won't work
unsigned char s_ASCII_MustBeLowerToUpper(unsigned char c)
{
    return c + ('A' - 'a');
}

#if 1
bool s_ASCII_IsAmbigNuc(unsigned char c)
{
    switch(c) {
    case 'U': case 'u':
    case 'R': case 'r':
    case 'Y': case 'y':
    case 'S': case 's':
    case 'W': case 'w':
    case 'K': case 'k':
    case 'M': case 'm':
    case 'B': case 'b':
    case 'D': case 'd':
    case 'H': case 'h':
    case 'V': case 'v':
    case 'N': case 'n':
        return true;
    default:
        return false;
    }
}
#endif

#if 0
inline static bool s_ASCII_IsUnAmbigNuc(unsigned char c)
{
    switch( c ) {
    case 'A': case 'C': case 'G': case 'T':
    case 'a': case 'c': case 'g': case 't':
        return true;
    default:
        return false;
    }
}
#endif

KHASH_SET_INIT_STR(SeqIdSet);

typedef char* charptr;
typedef kvec_t(charptr) vec_charptr;

typedef struct {
    vec_charptr name_list;
    khash_t(SeqIdSet)* name_set;
} CachedSeqNames;

CachedSeqNames*
CachedSeqNamesNew()
{
    CachedSeqNames* names = (CachedSeqNames*)calloc(1, sizeof(CachedSeqNames));
    kv_init(names->name_list);
    names->name_set = kh_init_SeqIdSet();
    return names;
}

CachedSeqNames*
CachedSeqNamesFree(CachedSeqNames* names)
{
    for (size_t i = 0; i < kv_size(names->name_list); ++i) {
        char* name = kv_A(names->name_list, i);
        free(name);
    }
    kv_destroy(names->name_list);
    kh_destroy_SeqIdSet(names->name_set);
    free(names);
    return NULL;
}

static BOOL 
name_is_cached(CachedSeqNames* names, const char* new_name)
{
    return kh_get_SeqIdSet(names->name_set, new_name) != kh_end(names->name_set);
}

void
cache_new_name(CachedSeqNames* names, const char* new_name)
{
    int n = strlen(new_name) + 1;
    char* s = (char*)malloc(n);
    memcpy(s, new_name, n);
    kv_push(charptr, names->name_list, s);
    int r;
    kh_put_SeqIdSet(names->name_set, s, &r);
    hbn_assert(r == 1);
}

HbnFastaReader*
HbnFastaReaderNew(const char* filename)
{
    HbnFastaReader* reader = (HbnFastaReader*)calloc(1, sizeof(HbnFastaReader));
    reader->filename = filename;
    reader->line_reader = HbnLineReaderNew(filename);
    ks_init(reader->name);
    ks_init(reader->comment);
    ks_init(reader->sequence);
    ks_init(reader->plus);
    ks_init(reader->qual);
    reader->skip_error_formated_sequences = FALSE;
    reader->cached_names = CachedSeqNamesNew();
    return reader;
}

HbnFastaReader*
HbnFastaReaderFree(HbnFastaReader* reader)
{
    reader->line_reader = HbnLineReaderFree(reader->line_reader);
    ks_destroy(reader->name);
    ks_destroy(reader->comment);
    ks_destroy(reader->sequence);
    ks_destroy(reader->plus);
    ks_destroy(reader->qual);
    reader->cached_names = CachedSeqNamesFree(reader->cached_names);
    free(reader);
    return NULL;
}

void
HbnFastaReaderClear(HbnFastaReader* reader)
{
    ks_clear(reader->name);
    ks_clear(reader->comment);
    ks_clear(reader->sequence);
    ks_clear(reader->plus);
    ks_clear(reader->qual);
}

kstring_t* HbnFastaReaderName(HbnFastaReader* reader)
{
    return &reader->name;
}

kstring_t* HbnFastaReaderComment(HbnFastaReader* reader)
{
    return &reader->comment;
}

kstring_t* HbnFastaReaderSequence(HbnFastaReader* reader)
{
    return &reader->sequence;
}

kstring_t* HbnFastaReaderPlusLine(HbnFastaReader* reader)
{
    return &reader->plus;
}

kstring_t* HbnFastaReaderQualityScoresLine(HbnFastaReader* reader)
{
    return &reader->qual;
}

void HbnFastaReaderSkipErrorFormatedSequences(HbnFastaReader* reader)
{
    reader->skip_error_formated_sequences = TRUE;
}

size_t HbnFastaReaderLineNumber(const HbnFastaReader* reader)
{
    return HbnLineReaderLineNumber(reader->line_reader);
}

static BOOL
s_parse_defline(kstring_t* line, HbnFastaReader* reader)
{
    const char* s = ks_s(*line);
    size_t n = ks_size(*line);
    hbn_assert(s[0] == '>' || s[0] == '@');
    if (n == 1) {
        kputc('\0', line);
        HBN_WARN("HbnFastaReader: Defline '%s' around (%zu, %s) lacks a proper ID",
            ks_s(*line),
            HbnFastaReaderLineNumber(reader),
            reader->filename);
        return FALSE;
    }
    size_t i = 1;
    while (i < n && (!isspace(s[i]))) ++i;
    kputsn(&s[1], i-1, &reader->name);
    while (i < n && isspace(s[i])) ++i;
    if (i < n) kputsn(&s[i], n-i, &reader->comment);
    if (reader->cached_names == NULL) return TRUE;

    kputc('\0', &reader->name);
    CachedSeqNames* cached_names = (CachedSeqNames*)reader->cached_names;
    if (name_is_cached(cached_names, ks_s(reader->name))) {
        HBN_WARN("HbnFastaReader: Seq-id '%s' is a duplicate around (%zu, %s).",
            ks_s(reader->name),
            HbnFastaReaderLineNumber(reader),
            reader->filename);
        return FALSE;
    }
    cache_new_name(cached_names, ks_s(reader->name));
    ks_pop_back(reader->name);

    return TRUE;
}

static BOOL
s_check_data_line(HbnFastaReader* reader, kstring_t* line)
{
    // make sure the first data line has at least SOME resemblance to
    // actual sequence data
    if (!ks_empty(reader->sequence)) return TRUE;

    size_t good = 0, bad = 0;
    // in case the data has huge sequences all on the first line we do need
    // a cutoff and "70" seems reasonable since it's the default width of
    // CFastaOstream (as of 2017-03-09)
    size_t len_to_check = hbn_min(ks_size(*line), 70);
    size_t ambig_nuc = 0;
    const char* s = ks_s(*line);
    for (size_t pos = 0; pos < len_to_check; ++pos) {
        unsigned char c = s[pos];
        if (s_ASCII_IsAlpha(c) || c == '*') {
            ++good;
            if (s_ASCII_IsAmbigNuc(c)) {
                ++ambig_nuc;
            }
        } else if (c == '-') {
            ++good;
        } else if (isspace(c) || (c >= '0' && c <= '9')) {
            // treat whilespace and digits as neutral
        } else if (c == ';') {
            break; // comment - ignore rest of line
        } else {
            ++bad;
        }
    }
    if (bad >= good / 3 &&
        (len_to_check > 3 || good == 0 || bad > good)) {
        HBN_WARN("HbnFastaReader: Near (%zu, %s), there's a line that doesn't look like"
            " plausible data, but it's not marked as defline or comment.",
            HbnLineReaderLineNumber(reader->line_reader),
            reader->filename);
        return FALSE;        
    }
    // warn if more than a certain percentage is ambiguous nucleotides
    const static size_t kWarnPercentAmbiguous = 40; // e.g. "40" means "40%"
    const size_t percent_ambig = (good == 0)?100:((ambig_nuc * 100) / good);
    if( len_to_check > 3 && percent_ambig > kWarnPercentAmbiguous ) {
        HBN_WARN( 
            "HbnFastaReader: Start of first data line (%zu, %s) in seq is about "
            "%zu%% ambiguous nucleotides (shouldn't be over %zu%%)",
            HbnLineReaderLineNumber(reader->line_reader),
            reader->filename,
            percent_ambig,
            kWarnPercentAmbiguous);
    }
    return TRUE;
}

static BOOL 
s_parse_data_line(HbnFastaReader* reader, kstring_t* line)
{
    if (!s_check_data_line(reader, line)) return FALSE;
    // most lines won't have a comment (';') 
    // so optimize for that case as much as possible
    const size_t s_len = ks_size(*line);
    ks_reserve(&reader->sequence, ks_size(reader->sequence) + s_len);
    size_t seq_pos = ks_size(reader->sequence);
    for (size_t pos = 0; pos < s_len; ++pos) {
        const unsigned char c = ks_A(*line, pos);
        if (c == ';') break;
        if (c == '\t' || c == '\v' || c == '\f' || c == '\r' || c == ' ') continue;
        if (nst_nt16_table[c] == 16) {
            char msg[512];
            if (isprint(c)) {
                sprintf(msg, "HbnFastaReader: "
                    "There are invalid nucleotide residue(s) '%c' in input sequence around (%zu, %zu, %s)", 
                    c, pos + 1, HbnLineReaderLineNumber(reader->line_reader), reader->filename);
            } else {
                sprintf(msg, "HbnFastaReader: "
                    "There are invalid nucleotide residue(s) (ascii code: %d) in input sequence around (%zu, %zu, %s)", 
                    c, pos + 1, HbnLineReaderLineNumber(reader->line_reader), reader->filename);
            }
            HBN_WARN(msg);
            return FALSE;
        }
        ks_A(reader->sequence, seq_pos) = c;
        ++seq_pos;
    }
    ks_set_size(&reader->sequence, seq_pos);
    return TRUE;
}

static void
s_advance_to_next_header(HbnFastaReader* reader)
{
    kstring_t* line = HbnLineReaderLine(reader->line_reader);
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue;
        int c = ks_front(*line);
        if (c == '>' || c == '@') {
            HbnLineReaderUngetline(reader->line_reader);
            break;
        }
    }
}

static BOOL
s_parse_fastq_sequence(HbnFastaReader* reader)
{
    /// header
    int find_header = 0;
    kstring_t* line = HbnLineReaderLine(reader->line_reader);
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue;
        if (ks_front(*line) == '@') {
            find_header = 1;
            break;
        }
    }
    
    if (!find_header) return TRUE;
    s_parse_defline(line, reader);

    /// sequence
    int find_sequence = 0;
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue;
        find_sequence = 1;
        break;
    }
    if (!find_sequence) {
        HBN_WARN("HbnFastaReader: sequence data is missing from FASTQ sequence at around (%zu, %s)",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        return FALSE;
    }
    if (!s_parse_data_line(reader, line)) return FALSE;

    /// '+'
    int find_plus = 0;
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue;
        find_plus = 1;
        break;
    }
    if (!find_plus) {
        HBN_WARN("HbnFastaReader: plus data line (which shall start with '+') is missing from FASTQ sequence at around (%zu, %s)",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        return FALSE;        
    }
    if (ks_front(*line) != '+') {
        HBN_WARN("HbnFastaReader: detect an invalid plus data (which shall start with '+') "
            "line from FASTQ sequence at around (%s, %d).",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        return FALSE;
    }

    /// qual line
    int find_qual = 0;
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue;
        find_qual = 1;
        break;
    }
    if (!find_qual) {
        HBN_WARN("HbnFastaReader: "
            "quality scores line is missing from FASTQ sequence at around (%zu, %s)",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        return FALSE;
    }
    if (ks_size(*line) != ks_size(reader->sequence)) {
        HBN_WARN("HbnFastaReader: "
            "quality line and sequence line do not have the same length at around (%s, %d).",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        return FALSE;            
    }    

    return TRUE;
}

#define handle_read_error do { \
    if (reader->skip_error_formated_sequences) { \
        HBN_WARN("We skip this sequence and restart scanning from the next one."); \
        s_advance_to_next_header(reader); \
    } else { \
        HBN_ERR("We terminate the scanning process and exit with EXIT_FAILURE."); \
    } \
} while(0)

BOOL HbnFastaReaderReadOneSeq(HbnFastaReader* reader)
{
    HbnFastaReaderClear(reader);
    BOOL need_defline = TRUE;
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        HbnLineReaderReadOneLine(reader->line_reader);
        kstring_t* line = HbnLineReaderLine(reader->line_reader);
        truncate_end_spaces(line);
        if (ks_empty(*line)) continue; //  ignore lines containing only whitespace
        const char c = ks_front(*line);
        if (c == '>')  {
            if (need_defline) {
                if (!s_parse_defline(line, reader)) {
                    handle_read_error;
                    need_defline = TRUE;
                    continue;
                }
                need_defline = FALSE;
                continue;
            } else {
                HbnLineReaderUngetline(reader->line_reader);
                // start of the next sequence;
                break;
            }
        } else if (c == '@') {
            HbnLineReaderUngetline(reader->line_reader);
            if (!s_parse_fastq_sequence(reader)) {
                handle_read_error;
                need_defline = TRUE;
                continue;
            }
            need_defline = FALSE;
            break;
        }

        if (c == '!' || c == '#' || c == ';') {
            // no content, just a comment or blank line
            continue;
        } else if (need_defline) {
            HbnLineReaderUngetline(reader->line_reader);
            HBN_WARN("HbnFastaReader: Input doesn't start with a defline or comment around (%zu, %s)",
                HbnLineReaderLineNumber(reader->line_reader),
                reader->filename);
            handle_read_error;
            need_defline = TRUE;
        }

        if (!s_parse_data_line(reader, line)) {
            handle_read_error;
            need_defline = TRUE;
        }
    }

    if (need_defline && HbnLineReaderAtEof(reader->line_reader)) {
        HBN_WARN("HbnFastaReader: Expected defline around (%zu, %s)",
            HbnLineReaderLineNumber(reader->line_reader), reader->filename);
        if (reader->skip_error_formated_sequences) {
            HBN_WARN("We skipt this error.");
            return FALSE;
        } else {
            HBN_ERR("We terminate the scanning process and exit with EXIT_FAILURE.");
        }
    }

    return TRUE;
}