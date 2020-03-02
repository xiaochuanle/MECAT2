#ifndef __HBN_LOOKUP_TABLE_H
#define __HBN_LOOKUP_TABLE_H

#include "../corelib/hbn_aux.h"
#include "../corelib/seqdb.h"
#include "../ncbi_blast/setup/blast_sequence_blk.h"
#include "../ncbi_blast/setup/blast_query_info.h"

#ifdef __cplusplus
extern "C" {
#endif

#define OffsetBits	((u64)34)
#define HashBits	((u64)30)
#define CntBits		((u64)30)
#define OffsetMask	((U64_ONE << OffsetBits)-1)

#define KmerStats_Offset(u) ((u)&OffsetMask)
#define KmerStats_Cnt(u)	((u)>>OffsetBits)

typedef struct {
    u64* offset_list;
    void* kmer_stats;
} LookupTable;

u64*
extract_kmer_list(const LookupTable* lktbl, const u64 hash, u64* n);

LookupTable*
destroy_lookup_table(LookupTable* lktbl);

LookupTable*
build_lookup_table(const text_t* db,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads);

LookupTable*
build_lookup_table_from_seq_chunk(BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads);

#ifdef __cplusplus
}
#endif

#endif // __HBN_LOOKUP_TABLE_H