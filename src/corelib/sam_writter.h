#ifndef __SAM_WRITTER_H
#define __SAM_WRITTER_H

#include "seqdb.h"
#include "hbn_hit.h"

#ifdef __cplusplus
extern "C" {
#endif

void
print_sam_header(FILE* out);

void
print_sam_program(int argc, char* argv[], FILE* out);

void
print_cigar(HbnHSP* hsp, kstring_t* aligned_strings, int print_clip_info, kstring_t* out_buf);

void
print_sam_result(HbnHSP* hsp,
    const text_t* queries,
    const text_t* database,
    kstring_t* aligned_strings,
    kstring_t* out_buf);

#ifdef __cplusplus
}
#endif

#endif // __SAM_WRITTER_H