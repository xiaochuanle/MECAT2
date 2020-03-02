#include "mecat_results.h"

#include "../../corelib/m4_record.h"
#include "../../ncbi_blast/setup/hsp2string.h"

const char* kSamVersion = "1.4";

void
print_sam_prolog(FILE* out, const char* sam_version, const char* prog_version, int argc, char* argv[])
{
    fprintf(out, "@HD\t");
    fprintf(out, "VN:%s\t", sam_version);
    fprintf(out, "SO:unknown\t");
    fprintf(out, "GO:query\n");

    fprintf(out, "@PG\t");
    fprintf(out, "ID:0\t");
    fprintf(out, "VN:%s\t", prog_version);
    fprintf(out, "CL:");
    for (int i = 0; i < argc - 1; ++i) fprintf(out, "%s ", argv[i]);
    fprintf(out, "%s", argv[argc-1]);
    fprintf(out, "\t");
    fprintf(out, "PN:mecat2map\n");
}

static void
print_cigar_core(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    kstring_t* out)
{
    int aln_size = hsp->hsp_info.subject_align_offset - hsp->hsp_info.query_align_offset;
    const char* q = ks_s(*aligned_strings) + hsp->hsp_info.subject_align_offset;
    const char* s = ks_s(*aligned_strings) + hsp->hsp_info.query_align_offset;
    int i = 0;
    while (i < aln_size) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == GAP_CHAR) {  // delete from subject (gap in query)
            type = 'D';
            while (i < aln_size && q[i] == GAP_CHAR) {
                hbn_assert(s[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (s[i] == GAP_CHAR) { // insert into subject (gap in subject)
            type = 'I';
            while (i < aln_size && s[i] == GAP_CHAR) {
                hbn_assert(q[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else { // substitution
            type = 'M';
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (i < aln_size && q[i] != GAP_CHAR && s[i] != GAP_CHAR) {
                ++i;
                ++cnt;
            }
        }
        ksprintf(out, "%d%c", cnt, type);
    }
}

void
print_sam_cigar(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    kstring_t* out)
{
    if (hsp->hbn_query.offset) ksprintf(out, "%zuS", hsp->hbn_query.offset);
    print_cigar_core(hsp, aligned_strings, out);
    if (hsp->hbn_query.end < hsp->hbn_query.seq_size) 
        ksprintf(out, "%zuS", hsp->hbn_query.seq_size - hsp->hbn_query.end);
}

void
print_paf_cigar(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    kstring_t* out)
{
    ksprintf(out, "cg:Z:");
    print_cigar_core(hsp, aligned_strings, out);
}

void
print_md(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    kstring_t* out)
{
    ksprintf(out, "MD:Z:");
    int aln_size = hsp->hsp_info.subject_align_offset - hsp->hsp_info.query_align_offset;
    const char* q = ks_s(*aligned_strings) + hsp->hsp_info.subject_align_offset;
    const char* s = ks_s(*aligned_strings) + hsp->hsp_info.query_align_offset;
    int i = 0, l_md = 0;
    while (i < aln_size) {
        int j = i;
        if (s[i] == GAP_CHAR) { // insert into subject (gap in  query)
            while (j < aln_size && s[j] == GAP_CHAR) {
                hbn_assert(q[j] != GAP_CHAR);
                ++j;
            }
        } else if (q[i] == GAP_CHAR) { // delete from subject (gap in subject)
            while (j < aln_size && q[j] == GAP_CHAR) {
                hbn_assert(s[j] != GAP_CHAR);
                ++j;
            }
            ksprintf(out, "%d^", l_md);
            kputsn(&s[i], j - i, out);
            l_md = 0;
        } else { // substitution
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (j < aln_size && q[j] != GAP_CHAR && s[j] != GAP_CHAR) ++j;
            for (int k = i; k < j; ++k) {
                if (q[k] != s[k]) {
                    ksprintf(out, "%d%c", l_md, s[k]);
                    l_md = 0;
                } else {
                    ++l_md;
                }
            }
        }
        i = j;
    }   
    if (l_md > 0) ksprintf(out, "%d", l_md); 
}

#if 0
                                         1 | string | Query sequence name                                     |
                                      |  2 |  int   | Query sequence length                                   |
                                      |  3 |  int   | Query start coordinate (0-based)                        |
                                      |  4 |  int   | Query end coordinate (0-based)                          |
                                      |  5 |  char  | `+' if query/target on the same strand; `-' if opposite |
                                      |  6 | string | Target sequence name                                    |
                                      |  7 |  int   | Target sequence length                                  |
                                      |  8 |  int   | Target start coordinate on the original strand          |
                                      |  9 |  int   | Target end coordinate on the original strand            |
                                      | 10 |  int   | Number of matching bases in the mapping                 |
                                      | 11 |  int   | Number bases, including gaps, in the mapping            |
                                      | 12 |  int   | Mapping quality (0-255 with 255 for missing)  
#endif 

void
print_one_paf_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out)
{
    const char tab = '\t';
    ksprintf(out, "%s", qname); /// 1) query name
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_query.seq_size); /// 2) query length
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_query.offset); /// 3) query start coordinate (0-based)
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_query.end); /// 4) query end coordinate (0-based)
    kputc(tab, out);
    /// 5) '+' if query and subject on the same strand; '-' if opposite
    if (hsp->hbn_query.strand == hsp->hbn_subject.strand) {
        kputc('+', out);
    } else {
        kputc('-', out);
    }
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 6) subject name
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_subject.seq_size); /// 7) subject length
    kputc(tab,out);
    ksprintf(out, "%zu", hsp->hbn_subject.offset); /// 8) subject start coordinate
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_subject.end); /// 9) subject end coordinate
    kputc(tab, out);
    
    int num_ident = 0;
    int aln_size = hsp->hsp_info.subject_align_offset - hsp->hsp_info.query_align_offset;
    const char* q = ks_s(*aligned_strings) + hsp->hsp_info.subject_align_offset;
    const char* s = ks_s(*aligned_strings) + hsp->hsp_info.query_align_offset;
    for (int i = 0; i < aln_size; ++i) if (q[i] == s[i]) ++num_ident;
    ksprintf(out, "%d", num_ident); /// 10) number of matching bases in the alignment
    kputc(tab, out);
    ksprintf(out, "%d", aln_size); /// 11) number of bases in the alignment (including gaps and mismatch bases)
    kputc(tab, out);
    ksprintf(out, "%d", 60); /// 12) mapq
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", hsp->hsp_info.chain_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", hsp->hsp_info.raw_score); ///  dp score
    if (dump_cigar) {
        kputc(tab, out);
        print_paf_cigar(hsp, aligned_strings, out);
    }
    if (dump_md) {
        kputc(tab, out);
        print_md(hsp, aligned_strings, out);
    }
    kputc('\n', out);
}

void
print_one_sam_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_md,
    kstring_t* out)
{
    const char tab = '\t';
    int flag = 0;
    if (hsp->hbn_query.strand ==  REV) flag |= 0x10; // reverse query strand
    ksprintf(out, "%s", qname); /// 1) query name
    kputc(tab, out);
    ksprintf(out, "%d", flag); /// 2) flag
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 3) subject name
    kputc(tab, out);
    ksprintf(out, "%zu", hsp->hbn_subject.offset + 1); /// 4) left most subject position (1-based)
    kputc(tab, out);
    ksprintf(out, "%d", 60); /// 5) mapq
    kputc(tab, out);
    print_sam_cigar(hsp, aligned_strings, out); /// 6) cigar
    kputc(tab, out);
    ksprintf(out, "*"); /// 7) rnext
    kputc(tab, out);
    ksprintf(out, "*"); /// 8) pnext
    kputc(tab, out);
    ksprintf(out, "%d", 0); /// 9) subject length
    kputc(tab, out);
    /// 10) aligned subsequence
    int num_ident = 0;
    int aln_size = hsp->hsp_info.subject_align_offset - hsp->hsp_info.query_align_offset;
    const char* q = ks_s(*aligned_strings) + hsp->hsp_info.subject_align_offset;
    const char* s = ks_s(*aligned_strings) + hsp->hsp_info.query_align_offset;
    for (int i = 0; i < aln_size; ++i) {
        if (q[i] == s[i]) ++num_ident;
        if (q[i] != GAP_CHAR) kputc(q[i], out);
    }
    kputc(tab, out);
    ksprintf(out, "*"); /// 11) quality score
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", hsp->hsp_info.chain_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", hsp->hsp_info.raw_score); ///  dp score
    if (dump_md) {
        kputc(tab, out);
        print_md(hsp, aligned_strings, out);
    }
    kputc('\n', out);
}

void
print_one_m4_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_cigar,
    const BOOL dump_md,
    const BOOL binary,
    kstring_t* line,
    kstring_t* out)
{
    if (!binary) {
        ks_clear(*line);
        blasthsp_to_string(hsp, line, qname, sname);
        if (dump_cigar || dump_md) ks_pop_back(*line);
        if (dump_cigar) {
            kputc('\t', line);
            print_paf_cigar(hsp, aligned_strings, line);
        }
        if (dump_md) {
            kputc('\t', line);
            print_md(hsp, aligned_strings, line);
        }
        if (dump_cigar || dump_md) kputc('\n', line);
        kputsn(ks_s(*line), ks_size(*line), out);
    } else {
        M4Record m4;
        m4.qid = hsp->hbn_query.oid;
        m4.qdir = hsp->hbn_query.strand;
        m4.qoff = hsp->hbn_query.offset;
        m4.qend = hsp->hbn_query.end;
        m4.qsize = hsp->hbn_query.seq_size;
        m4.sid = hsp->hbn_subject.oid;
        m4.sdir = hsp->hbn_subject.strand;
        m4.soff = hsp->hbn_subject.offset;
        m4.send = hsp->hbn_subject.end;
        m4.ssize = hsp->hbn_subject.seq_size;
        m4.ident_perc = hsp->hsp_info.perc_identity;
        m4.score = hsp->hsp_info.raw_score;
        const char* input = (const char*)(&m4);
        const int input_len = sizeof(M4Record);
        kputsn(input, input_len, out);
    }
}