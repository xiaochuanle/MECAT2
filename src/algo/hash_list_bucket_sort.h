#ifndef __HASH_LIST_BUCKET_SORT_H
#define __HASH_LIST_BUCKET_SORT_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define RADIX_SORT_QUIET

typedef u64 (*ValueExtractor)(void* list, const u64 i);
typedef void (*SetListValue)(void* src, const u64 src_idx, void* dst, const u64 dst_idx);

void
radix_sort(void* src, 
		   const u64 item_size,
		   const u64 item_count, 
		   const int num_threads,
		   ValueExtractor offset_extractor, 
		   ValueExtractor hash_extractor,
		   SetListValue set_list_value);

#ifdef __cplusplus
}
#endif

#endif // __HASH_LIST_BUCKET_SORT_H