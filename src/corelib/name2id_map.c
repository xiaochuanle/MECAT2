#include "name2id_map.h"

#include "khash.h"
#include "hbn_aux.h"

KHASH_MAP_INIT_STR(name2id_map_t, int);

NameToIdMap*
NameToIdMapNew()
{
    NameToIdMap* map = (NameToIdMap*)calloc(1, sizeof(NameToIdMap));
    map->name2id_map = kh_init_name2id_map_t();
    map->id2name_map = NULL;
    map->num_seqs = 0;
    return map;
}

NameToIdMap*
NameToIdMapFree(NameToIdMap* map)
{
    if (map->name2id_map) kh_destroy(name2id_map_t, map->name2id_map);
    if (map->id2name_map) free(map->id2name_map);
    free(map);
    return NULL;
}

void
NameToIdMapSet(const char* names, const int num_seqs, NameToIdMap* map)
{
    if (map->id2name_map) free(map->id2name_map);
    map->id2name_map = (const char**)malloc(sizeof(const char*) * num_seqs);
    khash_t(name2id_map_t)* name2id_map = (khash_t(name2id_map_t)*)(map->name2id_map);
    const char* p = names;
    for (int i = 0; i < num_seqs; ++i) {
        const char* name = p;
        khiter_t pos = kh_get_name2id_map_t(name2id_map, name);
        if (pos != kh_end(name2id_map)) HBN_ERR("sequence name '%s' already exists", name);
        int r = 0;
        pos = kh_put(name2id_map_t, name2id_map, name, &r);
        hbn_assert(r == 1);
        kh_value(name2id_map, pos) = i;
        map->id2name_map[i] = name;
        while (*p) ++p;
        ++p;
    }
    map->num_seqs = num_seqs;
}

int Name2IdMap_name2id(NameToIdMap* map, const char* name)
{
    khash_t(name2id_map_t)* name2id_map = (khash_t(name2id_map_t)*)(map->name2id_map);
    khiter_t pos = kh_get(name2id_map_t, name2id_map, name);
    hbn_assert(pos != kh_end(name2id_map));
    return kh_value(name2id_map, pos);
}

const char* Name2IdMap_id2name(NameToIdMap* map, const int id)
{
    hbn_assert(id < map->num_seqs, "id = %d, num_seqs = %d", id, map->num_seqs);
    return map->id2name_map[id];
} 