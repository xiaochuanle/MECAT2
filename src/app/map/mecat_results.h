#ifndef __MECAT_RESULTS_H
#define __MECAT_RESULTS_H

#include "hbn_options.h"
#include "../../ncbi_blast/setup/blast_hits.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const char* kSamVersion;

void
print_sam_prolog(FILE* out, const char* sam_version, const char* prog_version, int argc, char* argv[]);

void
print_one_paf_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out);

void
print_one_sam_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_md,
    kstring_t* out);

void
print_one_m4_result(const BlastHSP* hsp, 
    const kstring_t* aligned_strings, 
    const char* qname,
    const char* sname,
    const BOOL dump_cigar,
    const BOOL dump_md,
    const BOOL binary,
    kstring_t* line,
    kstring_t* out);

#ifdef __cplusplus
}
#endif

#endif // __MECAT_RESULTS_H