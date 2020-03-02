#include "small_object_alloc.h"

static const u32 kMemChunkSize = 8 * 1024 * 1024;

static MemoryChunk*
MemoryChunkNew(const u32 object_size)
{
    MemoryChunk* chunk = (MemoryChunk*)malloc(sizeof(MemoryChunk));
    chunk->data = (char*)malloc(kMemChunkSize);
    chunk->avail_data = 0;
    chunk->object_size = object_size;
    return chunk;
}

static MemoryChunk*
MemoryChunkFree(MemoryChunk* chunk)
{
    free(chunk->data);
    free(chunk);
    return NULL;
}

static void
MemoryChunkClear(MemoryChunk* chunk)
{
    chunk->avail_data = 0;
}

static void*
MemoryChunkAlloc(MemoryChunk* chunk, const u32 num_objects)
{
    u32 alloc_bytes = num_objects * chunk->object_size;
    if (chunk->avail_data + alloc_bytes > kMemChunkSize) return NULL;
    char* p = chunk->data + chunk->avail_data;
    chunk->avail_data += alloc_bytes;
    return (void*)(p);
}

static u32
round_up_16(const u32 s)
{
    const u32 mask = 15;
    u32 r = 16 - (s & mask);
    return s + r;
}

SmallObjectAlloc*
SmallObjectAllocNew(const u32 object_size)
{
    SmallObjectAlloc* alloc = (SmallObjectAlloc*)malloc(sizeof(SmallObjectAlloc));
    kv_init(alloc->chunk_list);
    alloc->first_avail_chunk = 0;
    alloc->object_size = round_up_16(object_size);
    MemoryChunk* chunk = MemoryChunkNew(alloc->object_size);
    kv_push(MemoryChunk*, alloc->chunk_list, chunk);
    return alloc;
}

SmallObjectAlloc*
SmallObjectAllocFree(SmallObjectAlloc* alloc)
{
    for (size_t i = 0; i < kv_size(alloc->chunk_list); ++i) {
        MemoryChunkFree(kv_A(alloc->chunk_list, i));
    }
    kv_destroy(alloc->chunk_list);
    free(alloc);
    return NULL;
}

void
SmallObjectAllocClear(SmallObjectAlloc* alloc)
{
    for (size_t i = 0; i <= alloc->first_avail_chunk; ++i) {
        MemoryChunkClear(kv_A(alloc->chunk_list, i));
    }
    alloc->first_avail_chunk = 0;
}

void*
SmallObjectAllocAlloc(SmallObjectAlloc* alloc, const u32 num_objects)
{
    void* p = MemoryChunkAlloc(kv_A(alloc->chunk_list, alloc->first_avail_chunk), num_objects);
    if (!p) {
        if (alloc->first_avail_chunk + 1 == kv_size(alloc->chunk_list)) {
            MemoryChunk* chunk = MemoryChunkNew(alloc->object_size);
            kv_push(MemoryChunk*, alloc->chunk_list, chunk);
        }
        ++alloc->first_avail_chunk;
        p = MemoryChunkAlloc(kv_A(alloc->chunk_list, alloc->first_avail_chunk), num_objects);
    }
    hbn_assert(p);
    return p;
}