#ifndef __WORD_FINDER_H
#define __WORD_FINDER_H

#include "../corelib/gapped_candidate.h"
#include "../corelib/seqdb.h"
#include "hbn_lookup_table.h"
#include "chain_dp.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BLOCK_DDF_KM_CNT    (80)

typedef struct {
    int qoff;
    idx soff;
} DDFKmerMatch;

typedef kvec_t(DDFKmerMatch) vec_ddfkm;

void ks_introsort_ddfkm_soff_lt(size_t n, DDFKmerMatch* a);

typedef struct {
    DDFKmerMatch ddfkm_array[BLOCK_DDF_KM_CNT];
    int ddfkm_count;
    int last_qoff;
} DDFKmerMatchBlock;

#define ddf_km_block_init(block) ((block).ddfkm_count = 0, (block).last_qoff = -1)

typedef struct {
    int block_index;
    int score;
} DDFKmerMatchBlockInfo;

void ks_introsort_kmblk_info_gt(size_t n, DDFKmerMatchBlockInfo* a);

typedef struct {
    DDFKmerMatchBlock*      _ddfkm_block_array;
    DDFKmerMatchBlock*      ddfkm_block_array;
    DDFKmerMatchBlockInfo*  ddfkm_block_info_array;
    int                     ddfkm_block_count;
} DDFKmerMatchBackbone;

void
DDFKmerMatchBackboneClear(DDFKmerMatchBackbone* backbone);

DDFKmerMatchBackbone*
DDFKmerMatchBackboneFree(DDFKmerMatchBackbone* backbone);

DDFKmerMatchBackbone*
DDFKmerMatchBackboneNew(const text_t* reference);

typedef struct {
    const text_t* reference;
    DDFKmerMatchBackbone* backbone;
    const LookupTable* lktbl;
    ChainWorkData* chain_data;
    int map_against_myself;
    int kmer_size;
    int window_size;
    int min_block_km;
    vec_int_pair seeding_subseqs;
    vec_u64 hash_list;
    vec_init_hit init_hit_list;
} WordFindData;

WordFindData*
WordFindDataNew(const text_t* reference,
    const LookupTable* lktbl,
    int kmer_size,
    int window_size,
    int min_block_km,
    int map_against_myself);

WordFindData*
WordFindDataFree(WordFindData* data);

void
ddfs_find_candidates(WordFindData* word_data,
    const u8* read,
    const int read_id,
    const int read_start_id,
    const int read_dir,
    const int read_size);

void set_kmer_block_size_info(const int desired_block_size);
int kmer_block_size_info_is_set();

#ifdef __cplusplus
}
#endif

#endif // __WORD_FINDER_H