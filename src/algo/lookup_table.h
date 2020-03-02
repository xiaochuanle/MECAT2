#ifndef __LOOKUP_TABLE_H
#define __LOOKUP_TABLE_H

#include "../corelib/seqdb.h"
#include "../ncbi_blast/setup/blast_sequence_blk.h"
#include "../ncbi_blast/setup/blast_query_info.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	u64*	offset_list;
	void*	kmer_stats;
} LookupTable;

#define OffsetBits	((u64)34)
#define HashBits	((u64)30)
#define CntBits		((u64)30)
#define OffsetMask	((U64_ONE << OffsetBits)-1)

#define KmerStats_Offset(u) ((u)&OffsetMask)
#define KmerStats_Cnt(u)	((u)>>OffsetBits)
#define PACK_OFFSET(h, o)	(((h) << OffsetBits) | (o))
#define OFFSET_HASH(u)		((u)>>OffsetBits)
#define OFFSET_OFFSET(u)	((u) & OffsetMask)

LookupTable*
build_lookup_table(const text_t* text, 
				   const int kmer_size, 
				   const int window_size,
				   const int num_threads,
				   const double repeat_frac);

LookupTable*
destroy_lookup_table(LookupTable* lktbl);

u64*
extract_kmer_list(const LookupTable* lktbl, const u64 hash, u64* n);

#ifdef __cplusplus
}
#endif

#endif // __LOOKUP_TABLE_H