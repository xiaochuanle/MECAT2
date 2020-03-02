#include "sam_writter.h"

#include <stdio.h>

#include "hbn_package_version.h"

void
print_sam_header(FILE* out)
{
    fprintf(out, "@HD\t");
    fprintf(out, "VN:1.4\tSO:unknown\tGO:query\n");
}

void
print_sam_program(int argc, char* argv[], FILE* out)
{
	fprintf(out, "@PG\tID:0\tVN:%s\tCL:", HBN_PACKAGE_VERSION);
	int i;
	for (i = 0; i < argc; ++i) fprintf(out, "%s ",argv[i]);
	fprintf(out, "\t");
	fprintf(out, "PN:%s\n", HBN_PACKAGE_NAME);
}

void
print_cigar(HbnHSP* hsp, kstring_t* aligned_strings, int print_clip_info, kstring_t* out_buf)
{
    if (print_clip_info && hsp->qoff) ksprintf(out_buf, "%zuH", hsp->qoff);
    int n = hsp->saln_offset - hsp->qaln_offset;
    const char* q = ks_s(*aligned_strings) + hsp->qaln_offset;
    const char* s = ks_s(*aligned_strings) + hsp->saln_offset;
    int i = 0;
    while (i < n) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == '-') { // delete from subject
            type = 'D';
            while (i < n && q[i] == '-') {
                hbn_assert(s[i] != '-');
                ++i;
                ++cnt;
            }
        } else if (s[i] == '-') { // insert into reference
            type = 'I';
            while (i < n && s[i] == '-') {
                hbn_assert(q[i] != '-');
                ++i;
                ++cnt;
            }
        } else { // match or mismatch
            hbn_assert(q[i] != '-' && s[i] != '-');
            type = 'M';
            while (i < n && q[i] != '-' && s[i] != '-') {
                ++i;
                ++cnt;
            }
        }
        ksprintf(out_buf, "%d%c", cnt, type);
    }
    if (print_clip_info && hsp->qend != hsp->qsize) ksprintf(out_buf, "%zuH", hsp->qsize - hsp->qend);
}

void
print_sam_result(HbnHSP* hsp,
    const text_t* queries,
    const text_t* database,
    kstring_t* aligned_strings,
    kstring_t* out_buf)
{
    normalise_hbn_hsp_offsets(hsp);
    normalise_hbn_hsp_sdir(hsp);
    int flag = 0;
    if (hsp->qdir == REV) flag = 0x10; // seq being reverse complemented
    ksprintf(out_buf, "%s\t", seqdb_seq_name(queries, hsp->qid)); // 1) query name
    ksprintf(out_buf, "%d\t", flag); // 2) flag
    ksprintf(out_buf, "%s\t", seqdb_seq_name(database, hsp->sid)); // 3) subject name
    ksprintf(out_buf, "%zu\t", hsp->soff + 1); // 4) 1-based left most position
    ksprintf(out_buf, "255\t"); // 5) mapq
    print_cigar(hsp, aligned_strings, 1, out_buf); // 6) cigar
    ksprintf(out_buf, "\t");
    ksprintf(out_buf, "*\t"); // 7) rnext
    ksprintf(out_buf, "*\t"); // 8) pnext
    ksprintf(out_buf, "0\t"); // 9) tlen
    // 10) seq
    int n = hsp->saln_offset - hsp->qaln_offset;
    const char* q = ks_s(*aligned_strings) + hsp->qaln_offset;
    for (int i = 0; i < n; ++i) if (q[i] != '-') ksprintf(out_buf, "%c", q[i]);
    ksprintf(out_buf, "\t");
    ksprintf(out_buf, "*\n"); // 11) qual
}