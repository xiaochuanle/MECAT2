#ifndef __SMALL_OBJECT_ALLOC_H
#define __SMALL_OBJECT_ALLOC_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char* data;
    u32 avail_data;
    u32 object_size;
} MemoryChunk;

typedef kvec_t(MemoryChunk*) vec_mem_chunk;

typedef struct {
    vec_mem_chunk   chunk_list;
    u32             first_avail_chunk;
    u32             object_size;
} SmallObjectAlloc;

SmallObjectAlloc*
SmallObjectAllocNew(const u32 object_size);

SmallObjectAlloc*
SmallObjectAllocFree(SmallObjectAlloc* alloc);

void
SmallObjectAllocClear(SmallObjectAlloc* alloc);

void*
SmallObjectAllocAlloc(SmallObjectAlloc* alloc, const u32 num_objects);

#ifdef __cplusplus
}
#endif

#endif // __SMALL_OBJECT_ALLOC_H