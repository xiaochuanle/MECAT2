#ifndef __HSP2STRING_H
#define __HSP2STRING_H

#include "blast_hits.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/name2id_map.h"
#include "../../corelib/line_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

const char* blasthsp_to_string_names(const BlastHSP* hsp, ...);

const char* blasthsp_to_string_name2idmap(const BlastHSP* hsp, ...);

const char* blasthsp_to_string_ids(const BlastHSP* hsp, ...);

#define blasthsp_to_string(hsp_ptr, ...) \
( \
    (HBN_MACRO_ARG_CNT(__VA_ARGS__) == 3) \
    ? \
    blasthsp_to_string_names(hsp_ptr, __VA_ARGS__) \
    : \
    ( \
        (HBN_MACRO_ARG_CNT(__VA_ARGS__) == 2) \
        ? \
        blasthsp_to_string_name2idmap(hsp_ptr, __VA_ARGS__) \
        : \
        blasthsp_to_string_ids(hsp_ptr, __VA_ARGS__) \
    ) \
)

#define dump_blasthsp(output_func, stream, hsp) do { \
    ks_dinit(hspstr_); \
    blasthsp_to_string(&(hsp), &hspstr_); \
    output_func(stream, "%s", ks_s(hspstr_)); \
    ks_destroy(hspstr_); \
} while(0)

BlastHSP* string_to_blasthsp_ids(const char* str, ...);

BlastHSP* string_to_blasthsp_name2idmap(const char* str, ...);

#define string2hsp(str, ...) \
( \
    (HBN_MACRO_ARG_CNT(__VA_ARGS__) == 1) \
    ? \
    string_to_blasthsp_ids(str, __VA_ARGS__) \
    : \
    string_to_blasthsp_name2idmap(str, __VA_ARGS__) \
)

typedef struct {
    NameToIdMap* query_name2id_map;
    NameToIdMap* subject_name2id_map;
    HbnLineReader* line_reader;
    BlastHSP hsp;
    int unget;
} BlastHSPReader;

BlastHSPReader*
BlastHSPReaderNew(const char* hsp_path, 
    const char* query_names,
    const int num_qureis,
    const char* subject_names,
    const int num_subjects);

BlastHSPReader*
BlastHSPReaderFree(BlastHSPReader* reader);

int
BlastHSPReaderGet(BlastHSPReader* reader, BlastHSP* hsp);

void
BlastHSPReaderUnget(BlastHSPReader* reader);

#ifdef __cplusplus
}
#endif

#endif // __HSP2STRING_H