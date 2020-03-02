#include "word_finder.h"

#include <math.h>

#include "../corelib/ksort.h"

static int block_size;
static int block_shift;
static idx block_mask;
static double s_ddfs_cutoff = .25;
static int block_size_info_is_set = 0;

void set_kmer_block_size_info(const int desired_block_size)
{
    int n = 1;
    int shift = 0;
    while (n < desired_block_size) {
        n <<= 1;
        shift++;
    }
    block_size = n;
    block_shift = shift;
    block_mask = n - 1;
    HBN_LOG("desired block size: %d, has been fixed to be (%d, %d, %d)", 
        desired_block_size, block_size, block_mask, block_shift);
    block_size_info_is_set = 1;
}

int kmer_block_size_info_is_set() 
{
    return block_size_info_is_set;
}

#define offset_2_blk_id(offset) ((offset)>>block_shift)
#define offset_2_blk_offset(offset) ((offset)&block_mask)
#define blk_id_2_offset(blk) ((blk)<<block_shift)
#define backbone_offset_2_offset(blk_id, blk_offset) (((blk_id)<<block_shift)|blk_offset)

#define ddfkm_soff_lt(a, b) ( \
    ((a).soff < (b).soff) \
    || \
    ((a).soff == (b).soff && (a).qoff < (b).qoff) \
)
KSORT_INIT(ddfkm_soff_lt, DDFKmerMatch, ddfkm_soff_lt);

#define kmblk_info_gt(a, b) ((a).score > (b).score)
KSORT_INIT(kmblk_info_gt, DDFKmerMatchBlockInfo, kmblk_info_gt);

void
DDFKmerMatchBackboneClear(DDFKmerMatchBackbone* backbone)
{
    for (int i = 0; i < backbone->ddfkm_block_count; ++i) {
        int block_index = backbone->ddfkm_block_info_array[i].block_index;
        ddf_km_block_init(backbone->ddfkm_block_array[block_index]);
    }
    backbone->ddfkm_block_count = 0;
}

DDFKmerMatchBackbone*
DDFKmerMatchBackboneFree(DDFKmerMatchBackbone* backbone)
{
    free(backbone->_ddfkm_block_array);
    free(backbone->ddfkm_block_info_array);
    free(backbone);
    return NULL;
}

DDFKmerMatchBackbone*
DDFKmerMatchBackboneNew(const text_t* reference)
{
    hbn_assert(kmer_block_size_info_is_set());
    size_t max_ref_size = seqdb_max_offset(reference);
    size_t num_block = offset_2_blk_id(max_ref_size) + 5;
    DDFKmerMatchBackbone* backbone = (DDFKmerMatchBackbone*)calloc(1, sizeof(DDFKmerMatchBackbone));
    backbone->_ddfkm_block_array = (DDFKmerMatchBlock*)calloc(num_block, sizeof(DDFKmerMatchBlock));
    for (size_t i = 0; i < num_block; ++i) ddf_km_block_init(backbone->_ddfkm_block_array[i]);
    backbone->ddfkm_block_array = backbone->_ddfkm_block_array + 1;
    backbone->ddfkm_block_info_array = (DDFKmerMatchBlockInfo*)calloc(num_block, sizeof(DDFKmerMatchBlockInfo));
    backbone->ddfkm_block_count = 0;
    return backbone;
}

static int
is_related(int qoff, 
    idx soff, 
    int seed_qoff, 
    idx seed_soff)
{
    int is_related = 0;
    if (qoff == seed_qoff && soff == seed_soff) {
        is_related = 1;
    } else if (qoff < seed_qoff && soff < seed_soff) {
        double dq = seed_qoff - qoff;
        double ds = seed_soff - soff;
        if (fabs(dq - ds) > 1500) {
            is_related = 0;
        } else {
            double s = fabs(dq/ds - 1.0);
            if (s < s_ddfs_cutoff) {
                is_related = 1;
            }
        }
    } else if (qoff > seed_qoff && soff > seed_soff) {
        double dq = seed_qoff - qoff;
        double ds = seed_soff - soff;
        if (fabs(dq - ds) > 1500) {
            is_related = 0;
        } else {
            double s = fabs(dq/ds - 1.0);
            if (s < s_ddfs_cutoff) {
                is_related = 1;
            }
        }
    }
    return is_related;
}

static int
extract_hash_values(const u8* read,
    const int read_size,
    const int kmer_size,
    const int window_size,
    vec_u64* hash_list)
{
    if (read_size < kmer_size) return 0;
    kv_clear(*hash_list);
    const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;

	if (!intersect) {
		for (size_t j = 0; j <= read_size - kmer_size; j += window_size) {
			u64 hash = 0;
			for (int k = 0; k < kmer_size; ++k) {
				size_t pos = j + k;
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
            kv_push(u64, *hash_list, hash);
		}
	} else {
		u64 hash = 0;
		for (int j = 0; j < kmer_size; ++j) {
			size_t pos = j;
			u8 c = read[pos];
            hbn_assert(c >= 0 && c < 4);
			hash = (hash << 2) | c;
		}
		kv_push(u64, *hash_list, hash);
		for (u64 j = window_size; j <= read_size - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
				size_t pos = j + k;
                hbn_assert(pos < read_size, "p = %d, read_size = %d, j = %d, k = %d, stride = %d", pos, read_size, j, k, stride);
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
			kv_push(u64, *hash_list, hash);
		}
	}
    return kv_size(*hash_list);
}

static void
add_one_block(DDFKmerMatchBackbone* backbone, const int block_index)
{
    backbone->ddfkm_block_info_array[backbone->ddfkm_block_count].block_index = block_index;
    backbone->ddfkm_block_count++;
}

static void
insert_one_ddfkm(DDFKmerMatch* km, DDFKmerMatchBackbone* backbone)
{
    size_t block_index = offset_2_blk_id(km->soff);
    DDFKmerMatchBlock* block = backbone->ddfkm_block_array + block_index;
    if (block->ddfkm_count >= BLOCK_DDF_KM_CNT) return;
    if (block->last_qoff == km->qoff) return;
    block->last_qoff = km->qoff;
    block->ddfkm_array[block->ddfkm_count] = *km;
    if (!block->ddfkm_count) add_one_block(backbone, block_index);
    ++block->ddfkm_count;
}

static void
collect_subseq_seeds(vec_u64* hash_list,
    const u8* read,
    const int read_from,
    const int read_to,
    const idx soff_max,
    const LookupTable* lktbl,
    const int kmer_size,
    const int window_size,
    DDFKmerMatchBackbone* backbone)
{
    const int SL = 150, SR = 200;
    //const int SL = 300, SR = 200;
    //const int SL = read_to - read_from, SR = 0;
    int n = read_to - read_from;
    int s = 0;
    while (s < n) {
        int e = s + SL;
        e = hbn_min(e, n);
        int n_kmer = extract_hash_values(read + read_from + s, e - s, kmer_size, window_size, hash_list);
        DDFKmerMatch ddfkm;
        for (int i = 0; i < n_kmer; ++i) {
            u64 n_km;
            u64* km_list = extract_kmer_list(lktbl, kv_A(*hash_list, i), &n_km);
            int qoff = read_from + s + i * window_size;
            ddfkm.qoff = qoff;
            for (u64 k = 0; k < n_km; ++k) {
                idx x = km_list[k];
                if (x >= soff_max) continue;
                //if (x <= 9332) HBN_LOG("find kmer match: [%d, %d]", qoff, x);
                ddfkm.soff = x;
                insert_one_ddfkm(&ddfkm, backbone);
            }
        }
        s = e + SR;                
    }
}

static void
collect_seeds(vec_u64* hash_list,
    const u8* read,
    const int read_id,
    const int read_start_id,
    const int read_size,
    const text_t* reference,
    const LookupTable* lktbl,
    const int kmer_size,
    const int window_size,
    const int map_against_myself,
    vec_int_pair* seeding_regions,
    DDFKmerMatchBackbone* backbone)
{
    idx soff_max = IDX_MAX;
    if (map_against_myself) {
        const int reference_start_id = reference->dbinfo.seq_start_id;
        int max_refid = reference_start_id + seqdb_num_seqs(reference);
        int g_read_id = read_id + read_start_id;
        if (g_read_id >= reference_start_id && g_read_id < max_refid) {
            soff_max = seqdb_seq_offset(reference, read_id);
        }
    }

    DDFKmerMatchBackboneClear(backbone);
    for (size_t s = 0; s < kv_size(*seeding_regions); ++s) {
        int from = kv_A(*seeding_regions, s).first;
        int to = kv_A(*seeding_regions, s).second;
        hbn_assert(to <= read_size);
        collect_subseq_seeds(hash_list,
            read,
            from,
            to,
            soff_max,
            lktbl,
            kmer_size,
            window_size,
            backbone);
    }

    for (int i = 0; i < backbone->ddfkm_block_count; ++i) {
        int block_index = backbone->ddfkm_block_info_array[i].block_index;
        DDFKmerMatchBlock* block = backbone->ddfkm_block_array + block_index;
        backbone->ddfkm_block_info_array[i].score = block->ddfkm_count + (block-1)->ddfkm_count;
    }
}

static void
calc_global_ddfs_info(const int seed_qoff,
    const idx seed_soff,
    const int read_size,
    const idx ref_off,
    const idx ref_end,
    int* ql,
    int* qr,
    idx* sl,
    idx* sr,
    int* min_block_id,
    int* max_block_id)
{
    const int E = 1000;
    idx ref_l = seed_soff - ref_off;
    idx ref_r = ref_end - seed_soff;
    int dq = hbn_min(seed_qoff, ref_l + E);
    int ds = hbn_min(seed_qoff + E, ref_l);
    *ql = seed_qoff - dq;
    *sl = seed_soff - ds;

    dq = hbn_min(read_size - seed_qoff, ref_r + E);
    ds = hbn_min(read_size - seed_qoff + E, ref_r);
    *qr = seed_qoff + dq;
    *sr = seed_soff + ds;

    *min_block_id = (*sl) / block_size;
    *max_block_id = (*sr) / block_size + 1;
}

#define lie_in_range(qx, sx, ql, qr, sl, sr) \
	( \
	  ((qx) >= (ql) && (qx) <= (qr)) \
	  && \
	  ((sx) >= (sl) && (sx) <= (sr)) \
	)

static int
scoring_one_km_block(DDFKmerMatchBlock* block, 
    int ql, int qr, idx sl, idx sr, 
    int seed_qoff, idx seed_soff)
{
    int score = 0;
    for (int i = 0; i < block->ddfkm_count; ++i) {
        int qoff = block->ddfkm_array[i].qoff;
        idx soff = block->ddfkm_array[i].soff;
        if (!lie_in_range(qoff, soff, ql, qr, sl, sr)) continue;
        if (is_related(qoff, soff, seed_qoff, seed_soff)) {
            ++score;
            block->ddfkm_array[i].qoff = -1;
        }
    }

    if (score == block->ddfkm_count) {
        block->ddfkm_count = 0;
    } else if (score > 0) {
        int k = 0;
        for (int i = 0; i < block->ddfkm_count; ++i) {
            if (block->ddfkm_array[i].qoff >= 0) block->ddfkm_array[k++] = block->ddfkm_array[i];
        }
        block->ddfkm_count = k;
    }
    return score;
}

static int
find_candidate_for_one_ddfkm_list(DDFKmerMatchBackbone* backbone,
    ChainWorkData* chain_data,
    int read_id,
    int read_dir,
    int read_size,
    int ref_id,
    size_t ref_off,
    size_t ref_size,
    size_t ref_end,
    DDFKmerMatch* ddfkm_array,
    int ddfkm_count,
    int kmer_size,
    vec_init_hit* init_hit_list)
{
    kv_clear(chain_data->seeds);
    ChainSeed seed;
    //if (ddfkm_array[0].soff < 9332) HBN_LOG("scoring");
    for (int i = 0; i < ddfkm_count; ++i) {
        seed.qoff = ddfkm_array[i].qoff;
        seed.soff = ddfkm_array[i].soff;
        seed.length = kmer_size;
        //if (ddfkm_array[0].soff < 9332) HBN_LOG("[%d, %d]", seed.qoff, seed.soff);
        kv_push(ChainSeed, chain_data->seeds, seed);
    }
    int best_km_index = -1;
    int best_km_score = -1;
    if (!find_best_kmer_match(chain_data, &best_km_index, &best_km_score)) return 0;
    int seed_qoff = kv_A(chain_data->seeds, best_km_index).qoff;
    idx seed_soff = kv_A(chain_data->seeds, best_km_index).soff;
    int ql, qr;
    idx sl, sr;
    int block_id_from, block_id_to;
    calc_global_ddfs_info(seed_qoff, seed_soff, read_size, ref_off, ref_end, 
        &ql, &qr, &sl, &sr, &block_id_from, &block_id_to);
    int score = 0;
    for (int i = block_id_from; i < block_id_to; ++i) {
        DDFKmerMatchBlock* block = backbone->ddfkm_block_array + i;
        score += scoring_one_km_block(block, ql, qr, sl, sr, seed_qoff, seed_soff);
    }

    HbnInitHit init_hit;
    init_hit.qid = read_id;
    init_hit.qdir = read_dir;
    init_hit.sid = ref_id;
    init_hit.sdir = FWD;
    init_hit.qoff = seed_qoff;
    init_hit.qsize = read_size;
    hbn_assert(seed_soff >= ref_off);
    init_hit.soff = seed_soff - ref_off;
    init_hit.ssize = ref_size;
    init_hit.score = score;
    kv_push(HbnInitHit, *init_hit_list, init_hit);
    return 1;
}

static int
find_candidate_for_one_block(const text_t* reference,
    DDFKmerMatchBackbone* backbone,
    ChainWorkData* chain_data,
    int block_index,
    const int ddf_score_cutoff,
    int read_id,
    int read_dir,
    int read_size,
    int kmer_size,
    vec_init_hit* init_hist_list)
{
    DDFKmerMatchBlock* block = backbone->ddfkm_block_array + block_index;
    if (block->ddfkm_count <= ddf_score_cutoff) return 0;

    DDFKmerMatch ddfkm_array[BLOCK_DDF_KM_CNT*2];
    int ddfkm_count = 0;
    if (block->ddfkm_count < 20) {
        --block;
        if (block->ddfkm_count) {
            ks_introsort_ddfkm_soff_lt(block->ddfkm_count, block->ddfkm_array);
            int i = 0;
            while (i < block->ddfkm_count) {
                int j = i + 1;
                while (j < block->ddfkm_count && block->ddfkm_array[i].soff == block->ddfkm_array[j].soff) ++j;
                if (j - i <= 3) {
                    for (int k = i; k < j; ++k) ddfkm_array[ddfkm_count++] = block->ddfkm_array[k];
                }
                i = j;
            }
        }
        ++block;
    }
    if (block->ddfkm_count) {
        ks_introsort_ddfkm_soff_lt(block->ddfkm_count, block->ddfkm_array);
        int i = 0;
        while (i < block->ddfkm_count) {
            int j = i + 1;
            while (j < block->ddfkm_count && block->ddfkm_array[i].soff == block->ddfkm_array[j].soff) ++j;
            if (j - i <= 3) {
                for (int k = i; k < j; ++k) ddfkm_array[ddfkm_count++] = block->ddfkm_array[k];
            }
            i = j;
        }
    }
    //ks_introsort_ddfkm_soff_lt(ddfkm_count, ddfkm_array);

    int i = 0;
    int added_can = 0;
    while (i < ddfkm_count) {
        int ref_id = seqdb_offset_to_seq_id(reference, ddfkm_array[i].soff);
        size_t ref_off = seqdb_seq_offset(reference, ref_id);
        size_t ref_size = seqdb_seq_size(reference, ref_id);
        size_t ref_end = ref_off + ref_size;
        int j = i + 1;
        while (j < ddfkm_count && ddfkm_array[j].soff < ref_end) ++j;
        if (j - i > ddf_score_cutoff) {
            added_can += find_candidate_for_one_ddfkm_list(backbone, 
                            chain_data,
                            read_id,
                            read_dir,
                            read_size,
                            ref_id,
                            ref_off,
                            ref_size,
                            ref_end,
                            ddfkm_array + i,
                            j - i,
                            kmer_size,
                            init_hist_list);
        }
        i = j;
    }
    return added_can;
}

void
ddfs_find_candidates(WordFindData* word_data,
    const u8* read,
    const int read_id,
    const int read_start_id,
    const int read_dir,
    const int read_size)
{
    hbn_assert(block_size_info_is_set);
    collect_seeds(&word_data->hash_list, 
        read, 
        read_id, 
        read_start_id, 
        read_size, 
        word_data->reference, 
        word_data->lktbl, 
        word_data->kmer_size, 
        1, 
        word_data->map_against_myself, 
        &word_data->seeding_subseqs,
        word_data->backbone);

    ks_introsort_kmblk_info_gt(word_data->backbone->ddfkm_block_count, word_data->backbone->ddfkm_block_info_array);
    int added_can = 0;
    for (int i = 0; i < word_data->backbone->ddfkm_block_count; ++i) {
        int r = (word_data->backbone->ddfkm_block_info_array[i].score >= word_data->min_block_km * 2)
                ||
                (read_size < block_size * 2 && word_data->backbone->ddfkm_block_info_array[i].score > word_data->min_block_km);
        if (!r) continue;
        added_can += find_candidate_for_one_block(word_data->reference,
            word_data->backbone, 
            word_data->chain_data,
            word_data->backbone->ddfkm_block_info_array[i].block_index,
            word_data->min_block_km,
            read_id,
            read_dir,
            read_size,
            word_data->kmer_size,
            &word_data->init_hit_list);      
    }
}

WordFindData*
WordFindDataNew(const text_t* reference,
    const LookupTable* lktbl,
    int kmer_size,
    int window_size,
    int min_block_km,
    int map_against_myself)
{
    WordFindData* data = (WordFindData*)calloc(1, sizeof(WordFindData));
    data->reference = reference;
    data->backbone = DDFKmerMatchBackboneNew(reference);
    data->lktbl = lktbl;
    data->chain_data = ChainWorkDataNew(min_block_km, min_block_km * kmer_size * 0.8);
    data->map_against_myself = map_against_myself;
    data->kmer_size = kmer_size;
    data->window_size = window_size;
    data->min_block_km = min_block_km;
    kv_init(data->seeding_subseqs);
    kv_init(data->hash_list);
    kv_init(data->init_hit_list);

    return data;
}

WordFindData*
WordFindDataFree(WordFindData* data)
{
    DDFKmerMatchBackboneFree(data->backbone);
    ChainWorkDataFree(data->chain_data);
    kv_destroy(data->seeding_subseqs);
    kv_destroy(data->hash_list);
    kv_destroy(data->init_hit_list);
    free(data);
    return NULL;
}