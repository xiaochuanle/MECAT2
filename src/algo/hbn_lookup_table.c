#include "hbn_lookup_table.h"

#include "../corelib/cstr_util.h"
#include "../corelib/ksort.h"
#include "../corelib/khash.h"
#include "hash_list_bucket_sort.h"

KHASH_MAP_INIT_INT64(KmerHashToOffsetMap, u64);

typedef struct {
    u64 hash;
    i64 offset;
} KmerHashAndOffset;

static u64 
calc_num_kmers(const text_t* db,
    const int kmer_size,
    const int window_size)
{
    u64 num_kmers = 0;
    const int num_subjects = seqdb_num_seqs(db);
    for (int i = 0; i < num_subjects; ++i) {
        u64 subject_size = seqdb_seq_size(db, i);
        if (subject_size < kmer_size) continue;
        u64 n = (subject_size - kmer_size) / window_size + 1;
        num_kmers += n;
    }
    return num_kmers;
}

static KmerHashAndOffset*
get_khao_array(const text_t* db,
    const int kmer_size,
    const int window_size,
    u64* khao_count)
{
    hbn_timing_begin(__FUNCTION__);

    hbn_assert(db->unpacked_seq != NULL);
    u64 num_kmers = calc_num_kmers(db, kmer_size, window_size);
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)calloc(num_kmers, sizeof(KmerHashAndOffset));
    const int kIntersect = kmer_size > window_size;
    const int kStride = kmer_size - window_size;
    const u64 kIntersectMask = kIntersect ? ((U64_ONE << (kStride<<1)) - 1) : 0;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1);
    const int num_subjects = seqdb_num_seqs(db);
    size_t cnt = 0;

    HBN_LOG("kmer size = %d, window_size = %d", kmer_size, window_size);

    for (int i = 0; i < num_subjects; ++i) {
        const u64 subject_size = seqdb_seq_size(db, i);
        if (subject_size < kmer_size) continue;
        const u64 start = seqdb_seq_offset(db, i);

        if (!kIntersect) {
            for (u64 j = 0; j <= subject_size - kmer_size; j += window_size) {
                u64 hash = 0;
                for (int k = 0; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    const u8 c = db->unpacked_seq[pos];
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        } else {
            u64 hash = 0;
            for (int j = 0; j < kmer_size; ++j) {
                const u64 pos = start + j;
                const u8 c = db->unpacked_seq[pos];
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            KmerHashAndOffset khao = { hash, start };
            khao_array[cnt++] = khao;
            for (u64 j = window_size; j <= subject_size - kmer_size; j += window_size) {
                hash &= kIntersectMask;
                for (int k = kStride; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    const u8 c = db->unpacked_seq[pos];
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        }
    }
    hbn_assert(cnt == num_kmers);
    *khao_count = num_kmers;
    hbn_timing_end(__FUNCTION__);
    return khao_array;
}

static u64 
calc_num_kmers_from_seq_chunk(BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info,
    const int kmer_size,
    const int window_size)
{
    u64 num_kmers = 0;
    for (int i = seq_info->first_context; i <= seq_info->last_context; ++i) {
        u64 subject_size = seq_info->contexts[i].query_length;
        if (subject_size < kmer_size) continue;
        u64 n = (subject_size - kmer_size) / window_size + 1;
        num_kmers += n;
    }
    return num_kmers;
}

static KmerHashAndOffset*
get_khao_array_from_seq_chunk(BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info,
    const int kmer_size,
    const int window_size,
    u64* khao_count)
{
    HBN_LOG("kmer_size = %d, window_size = %d", kmer_size, window_size);
    u64 num_kmers = calc_num_kmers_from_seq_chunk(seq_blk, seq_info, kmer_size, window_size);
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)calloc(num_kmers, sizeof(KmerHashAndOffset));
    const int kIntersect = kmer_size > window_size;
    const int kStride = kmer_size - window_size;
    const u64 kIntersectMask = kIntersect ? ((U64_ONE << (kStride<<1)) - 1) : 0;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1);
    size_t cnt = 0;

    for (int i = seq_info->first_context; i <= seq_info->last_context; ++i) {
        const u64 subject_size = seq_info->contexts[i].query_length;
        if (subject_size < kmer_size) continue;
        const u64 start = seq_info->contexts[i].query_offset;
        const u8* seq = seq_blk->sequence;

        if (!kIntersect) {
            for (u64 j = 0; j <= subject_size - kmer_size; j += window_size) {
                u64 hash = 0;
                for (int k = 0; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    const u8 c = seq[pos];
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        } else {
            u64 hash = 0;
            for (int j = 0; j < kmer_size; ++j) {
                const u64 pos = start + j;
                const u8 c = seq[pos];
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            KmerHashAndOffset khao = { hash, start };
            khao_array[cnt++] = khao;
            for (u64 j = window_size; j <= subject_size - kmer_size; j += window_size) {
                hash &= kIntersectMask;
                for (int k = kStride; k < kmer_size; ++k) {
                    const u64 pos = start + j + k;
                    const u8 c = seq[pos];
                    hash = (hash << 2) | c;
                }
                hbn_assert(hash < kMaxHashValue);
                KmerHashAndOffset khao = { hash, start + j };
                khao_array[cnt++] = khao;
            }
        }
    }
    hbn_assert(cnt == num_kmers);
    *khao_count = num_kmers;
    return khao_array;
}

static u64
remove_repetitive_kmers(KmerHashAndOffset* khao_array,
    const u64 khao_count,
    const int max_kmer_occ)
{
    u64 distinct_kmers = 0;
    u64 removed_distinct_kmers = 0;
    u64 removed_kemrs = 0;

    u64 i = 0;
    while (i < khao_count) {
        ++distinct_kmers;
        u64 j = i + 1;
        while (j < khao_count && khao_array[i].hash == khao_array[j].hash) ++j;
        u64 n = j - i;
        if (n > max_kmer_occ) {
            ++removed_distinct_kmers;
            removed_kemrs += n;
            for (u64 k = i; k < j; ++k) khao_array[k].offset = U64_MAX;
        }
        i = j;
    }

    {
        char buf1[64], buf2[64], buf3[64];
        u64_to_string_comma(khao_count, buf1);
        u64_to_string_comma(removed_kemrs, buf2);
        double perc = 100.0 * removed_kemrs / khao_count;
        double_to_string(perc, buf3);
        HBN_LOG("Total kmers: %s, %s (%s%%) are filtered out.", buf1, buf2, buf3);

        u64_to_string_comma(distinct_kmers, buf1);
        u64_to_string_comma(removed_distinct_kmers, buf2);
        perc = 100.0 * removed_distinct_kmers / distinct_kmers;
        double_to_string(perc, buf3);
        HBN_LOG("Distinct kmers: %s, %s (%s%%) are fltered out.", buf1, buf2, buf3);
    }

    i = 0;
    for (u64 k = 0; k < khao_count; ++k) {
        if (khao_array[k].offset != U64_MAX) khao_array[i++] = khao_array[k];
    }
    return i;
}

static void
build_lktbl_from_khao_array(KmerHashAndOffset* khao_array,
    u64 khao_count,
    const int max_kmer_occ,
    khash_t(KmerHashToOffsetMap)** hash_2_offset_map_pp,
    u64** offset_array_pp)
{
    khao_count = remove_repetitive_kmers(khao_array, khao_count, max_kmer_occ);
    khash_t(KmerHashToOffsetMap)* hash_2_offset_map = kh_init(KmerHashToOffsetMap);
    u64 i = 0;
    while (i < khao_count) {
        u64 j = i + 1;
        while (j < khao_count && khao_array[i].hash == khao_array[j].hash) ++j;
        u64 n = j - i;
        hbn_assert(n <= max_kmer_occ);
        khiter_t iter = kh_get(KmerHashToOffsetMap, hash_2_offset_map, khao_array[i].hash);
        hbn_assert(iter == kh_end(hash_2_offset_map));
        int r = 0;
        iter = kh_put(KmerHashToOffsetMap, hash_2_offset_map, khao_array[i].hash, &r);
        hbn_assert(r == 1);
        //u64 u = (i << OffsetBits) | n;
        u64 u = (n << OffsetBits) | i;
        kh_value(hash_2_offset_map, iter) = u;
        i = j;
    }

    u64* offset_array = (u64*)calloc(khao_count, sizeof(u64));
    for (i = 0; i < khao_count; ++i) offset_array[i] = khao_array[i].offset;
    
    *hash_2_offset_map_pp = hash_2_offset_map;
    *offset_array_pp = offset_array;
}

static u64 
hash_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].hash;
}

static u64 
offset_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].offset;
}

static void
set_khao_array_item_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
    KmerHashAndOffset* src_khao_array = (KmerHashAndOffset*)(src);
    KmerHashAndOffset* dst_khao_array = (KmerHashAndOffset*)(dst);
    dst_khao_array[dst_idx] = src_khao_array[src_idx];
}

LookupTable*
build_lookup_table(const text_t* db,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads)
{
    u64 khao_count = 0;
    KmerHashAndOffset* khao_array = get_khao_array(db, kmer_size, window_size, &khao_count);
    radix_sort(khao_array, 
        sizeof(KmerHashAndOffset), 
        khao_count, 
        num_threads,
        offset_extractor,
        hash_extractor,
        set_khao_array_item_value);
    for (u64 i = 0; i < khao_count - 1; ++i) hbn_assert(khao_array[i].hash <= khao_array[i+1].hash);
    khash_t(KmerHashToOffsetMap)* hash_2_offset_map = NULL;
    u64* offset_array = NULL;
    build_lktbl_from_khao_array(khao_array, khao_count, max_kmer_occ, &hash_2_offset_map, &offset_array);
    free(khao_array);

    LookupTable* lktbl = (LookupTable*)calloc(1, sizeof(LookupTable));
    lktbl->offset_list = offset_array;
    lktbl->kmer_stats = hash_2_offset_map;
    return lktbl;
}

LookupTable*
build_lookup_table_from_seq_chunk(BLAST_SequenceBlk* seq_blk,
    BlastQueryInfo* seq_info,
    const int kmer_size,
    const int window_size,
    const int max_kmer_occ,
    const int num_threads)
{
    u64 khao_count = 0;
    KmerHashAndOffset* khao_array = get_khao_array_from_seq_chunk(seq_blk, seq_info, kmer_size, window_size, &khao_count);
    radix_sort(khao_array, 
        sizeof(KmerHashAndOffset), 
        khao_count, 
        num_threads,
        offset_extractor,
        hash_extractor,
        set_khao_array_item_value);
    for (u64 i = 0; i < khao_count - 1; ++i) hbn_assert(khao_array[i].hash <= khao_array[i+1].hash);
    khash_t(KmerHashToOffsetMap)* hash_2_offset_map = NULL;
    u64* offset_array = NULL;
    build_lktbl_from_khao_array(khao_array, khao_count, max_kmer_occ, &hash_2_offset_map, &offset_array);
    free(khao_array);

    LookupTable* lktbl = (LookupTable*)calloc(1, sizeof(LookupTable));
    lktbl->offset_list = offset_array;
    lktbl->kmer_stats = hash_2_offset_map;
    return lktbl;
}

LookupTable*
destroy_lookup_table(LookupTable* lktbl)
{
	if (lktbl->offset_list) free(lktbl->offset_list);
	khash_t(KmerHashToOffsetMap)* hash_2_offset_map = (khash_t(KmerHashToOffsetMap)*)(lktbl->kmer_stats);
	kh_destroy(KmerHashToOffsetMap, hash_2_offset_map);
	free(lktbl);
	return 0;
}

u64*
extract_kmer_list(const LookupTable* lktbl, const u64 hash, u64* n)
{
	u64* list = 0;
	*n = 0;

	khash_t(KmerHashToOffsetMap)* hash_2_offset_map = (khash_t(KmerHashToOffsetMap)*)(lktbl->kmer_stats);
	khiter_t pos = kh_get(KmerHashToOffsetMap, hash_2_offset_map, hash);
	if (pos == kh_end(hash_2_offset_map)) return NULL;
	u64 u = kh_value(hash_2_offset_map, pos);
	u64 cnt = KmerStats_Cnt(u);
	u64 start = KmerStats_Offset(u);
	if (cnt) {
		list = lktbl->offset_list + start;
		*n = cnt;
	}
    return list;
}
	