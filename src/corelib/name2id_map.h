#ifndef __NAME2ID_MAP_H
#define __NAME2ID_MAP_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void* name2id_map;
    const char** id2name_map;
    int num_seqs;
} NameToIdMap;

const char* Name2IdMap_id2name(NameToIdMap* map, const int id);

int Name2IdMap_name2id(NameToIdMap* map, const char* name);

void
NameToIdMapSet(const char* names, const int num_seqs, NameToIdMap* map);

NameToIdMap*
NameToIdMapFree(NameToIdMap* map);

NameToIdMap*
NameToIdMapNew();

#ifdef __cplusplus
}
#endif

#endif // __NAME2ID_MAP_H