#ifndef __STRING2HSP_H
#define __STRING2HSP_H

#include "hbn_hit.h"
#include "line_reader.h"
#include "name2id_map.h"

#ifdef __cplusplus
extern "C" {
#endif

const char* hsp2string(const HbnHSP* hsp, kstring_t* hspstr);

HbnHSP* string2hsp(const char* str, HbnHSP* hsp);

HbnHSP* string2hsp_names(const char* str, 
            NameToIdMap* query_name2id_map, 
            NameToIdMap* subject_name2id_map,
            HbnHSP* hsp);

const char* hsp2string_names(const HbnHSP* hsp, 
    NameToIdMap* query_name2id_map, 
    NameToIdMap* subject_name2id_map,
    kstring_t* hspstr);

typedef struct {
    NameToIdMap* query_name2id_map;
    NameToIdMap* subject_name2id_map;
    HbnLineReader* line_reader;
    HbnHSP hsp;
    int unget;
} HbnHspReader;

HbnHspReader*
HbnHspReaderNew(const char* hsp_path, 
    const char* query_names,
    const int num_qureis,
    const char* subject_names,
    const int num_subjects);

HbnHspReader*
HbnHspReaderFree(HbnHspReader* reader);

int
HbnHspReaderGet(HbnHspReader* reader, HbnHSP* hsp);

void
HbnHspReaderUnget(HbnHspReader* reader);

#ifdef __cplusplus
}
#endif

#endif // __STRING2HSP_H