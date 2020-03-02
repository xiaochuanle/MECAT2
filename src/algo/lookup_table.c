#include "lookup_table.h"

#include "../corelib/ksort.h"
#include "hash_list_bucket_sort.h"
#include "../corelib/khash.h"

KHASH_MAP_INIT_INT64(KmerHash2OffsetMap, u64);

static u64
calc_num_kmers(const text_t* text,
	const int kmer_size,
	const int window_size)
{
	u64 nkm = 0;
	const int num_seqs = seqdb_num_seqs(text);
	for (int i = 0; i < num_seqs; ++i) {
		u64 size = seqdb_seq_size(text, i);
		if (size < kmer_size) continue;
		u64 n = (size - kmer_size) / window_size + 1;
		nkm += n;
	}	
	return nkm;
}

static u64*
get_offset_list(const text_t* text,
				const int kmer_size,
				const int window_size,
				size_t* num_kmers)
{
	hbn_timing_begin(__func__);
	hbn_assert(text->unpacked_seq != NULL);
	u64 nkm = calc_num_kmers(text, kmer_size, window_size);
	u64* offset_list = (u64*)calloc(nkm, sizeof(u64));
	const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;
	
	const int num_seqs = seqdb_num_seqs(text);
	size_t cnt = 0;
	for (int i = 0; i < num_seqs; ++i) {
		u64 size = seqdb_seq_size(text, i);
		if (size < kmer_size) continue;
		u64 start = seqdb_seq_offset(text, i);

		if (!intersect) {
			for (size_t j = 0; j <= size - kmer_size; j += window_size) {
				u64 hash = 0;
				for (int k = 0; k < kmer_size; ++k) {
					size_t pos = start + j + k;
					//u8 c = _get_pac(text->packed_seq, pos);
					u8 c = text->unpacked_seq[pos];
					hash = (hash << 2) | c;
				}
				u64 pos = start + j;
				offset_list[cnt++] = PACK_OFFSET(hash, pos);
			}
		} else {
			u64 hash = 0;
			for (int j = 0; j < kmer_size; ++j) {
				size_t pos = start + j;
				//u8 c = _get_pac(text->packed_seq, pos);
				u8 c = text->unpacked_seq[pos];
				hash = (hash << 2) | c;
			}
			u64 pos = start;
			offset_list[cnt++] = PACK_OFFSET(hash, pos);
			for (u64 j = window_size; j <= size - kmer_size; j += window_size) {
				hash &= intersect_mask;
				for (int k = stride; k < kmer_size; ++k) {
					size_t pos = start + j + k;
					//u8 c = _get_pac(text->packed_seq, pos);
					u8 c = text->unpacked_seq[pos];
					hash = (hash << 2) | c;
				}
				u64 pos = start + j;
				offset_list[cnt++] = PACK_OFFSET(hash, pos);
			}
		}
	}
	hbn_assert(cnt == nkm);
	*num_kmers = cnt;
	hbn_timing_end(__func__);
	return offset_list;
}

static khash_t(KmerHash2OffsetMap)*
get_kmer_counts(u64* offset_list,
				const size_t num_kmers,
				const int kmer_size,
				const double repeat_frac)
{
	hbn_timing_begin(__func__);
	const u64 kMaxHash = U64_ONE << (kmer_size << 1);
	khash_t(KmerHash2OffsetMap)* hash_2_offset_map = kh_init(KmerHash2OffsetMap);
	int r = 0;
	size_t i = 0;
	while (i < num_kmers) {
		u64 hash_i = OFFSET_HASH(offset_list[i]);
		size_t j = i + 1;
		while (j < num_kmers) {
			u64 hash_j = OFFSET_HASH(offset_list[j]);
			if (hash_i != hash_j) break;
			++j;
		}
		size_t n = j - i;
		khiter_t pos = kh_get(KmerHash2OffsetMap, hash_2_offset_map, hash_i);
		hbn_assert(pos == kh_end(hash_2_offset_map));
		pos = kh_put(KmerHash2OffsetMap, hash_2_offset_map, hash_i, &r);
		hbn_assert(r == 1);
		kh_value(hash_2_offset_map, pos) = n;
		i = j;
	}
	hbn_timing_end(__func__);
	return hash_2_offset_map;
}

void
build_kmer_starts(u64* offset_list, const u64 num_kmers, const int kmer_size, khash_t(KmerHash2OffsetMap)* hash_2_offset_map)
{
	hbn_timing_begin(__func__);
	
	u64 i = 0;
	while (i < num_kmers) {
		u64 hash = OFFSET_HASH(offset_list[i]);
		u64 j = i + 1;
		while (j < num_kmers) {
			u64 h = OFFSET_HASH(offset_list[j]);
			if (h != hash) break;
			++j;
		}
		khiter_t pos = kh_get(KmerHash2OffsetMap, hash_2_offset_map, hash);
		u64 n = 0;
		if (pos != kh_end(hash_2_offset_map)) n = kh_value(hash_2_offset_map, pos);
		if (n) {
			n = n << OffsetBits;
			n |= i;
			kh_value(hash_2_offset_map, pos) = n;
		}
		i = j;
	}
	
	hbn_timing_end(__func__);
}

typedef struct {
	u32 hash;
	int cnt;
} KmerInfo;

#define kmer_info_cnt_gt(a, b) ((a).cnt > (b).cnt)
KSORT_INIT(kmer_info_cnt_gt, KmerInfo, kmer_info_cnt_gt);

static void
remove_repetitive_kmers(u64* offset_list, const u64 num_kmers, double frac, khash_t(KmerHash2OffsetMap)* hash_2_offset_map)
{
	u64 distinct_kmers = 0;
	u64 i = 0;
	while (i < num_kmers) {
		u64 hash = OFFSET_HASH(offset_list[i]);
		u64 j = i + 1;
		while (j < num_kmers) {
			u64 h = OFFSET_HASH(offset_list[j]);
			if (h != hash) break;
			++j;
		}
		++distinct_kmers;
		i = j;
	}
	//HBN_LOG("number of distinct kmers: %zu", distinct_kmers);
	KmerInfo* vki = (KmerInfo*)calloc(distinct_kmers, sizeof(KmerInfo));
	i = 0;
	u64 k = 0;
	while (i < num_kmers) {
		u64 hash = OFFSET_HASH(offset_list[i]);
		u64 j = i + 1;
		while (j < num_kmers) {
			u64 h = OFFSET_HASH(offset_list[j]);
			if (h != hash) break;
			++j;
		}
		vki[k].cnt = j - i;
		vki[k].hash = hash;
		++k;
		i = j;
	}	
	hbn_assert(k == distinct_kmers);
	ks_introsort_kmer_info_cnt_gt(distinct_kmers, vki);

	u64 removed_distinct_kmers = distinct_kmers * frac;
	HBN_LOG("distinct_kmers = %zu, removed_distinct_kmers = %zu", distinct_kmers, removed_distinct_kmers);
	while (removed_distinct_kmers && vki[removed_distinct_kmers].cnt < 200) --removed_distinct_kmers;
	while (removed_distinct_kmers < distinct_kmers && vki[removed_distinct_kmers].cnt > 500) ++removed_distinct_kmers;
	u64 removed_kmers = 0;
	for (i = 0; i < removed_distinct_kmers; ++i) {
		removed_kmers += vki[i].cnt;
		khiter_t pos = kh_get(KmerHash2OffsetMap, hash_2_offset_map, vki[i].hash);
		hbn_assert(pos != kh_end(hash_2_offset_map));
		kh_value(hash_2_offset_map, pos) = 0;
	}

	HBN_LOG("kmer occ cutoff: %d", vki[removed_distinct_kmers].cnt);
	double p = 100.0 * removed_kmers / num_kmers;
	HBN_LOG("number of kmers: %zu, %zu (%.2f%%) are filtered out.", num_kmers, removed_kmers, p);
	p = 100.0 * removed_distinct_kmers / distinct_kmers;
	HBN_LOG("number of distinct kmers: %zu, %zu (%.2f%%) are filtered out.", 
		distinct_kmers, removed_distinct_kmers, p);
	free(vki);
}

void
clear_hash_in_offset_list(u64* offset_list, const u64 num_kmers, const u64 max_offset)
{
	hbn_timing_begin(__func__);
	
	for (u64 i = 0; i != num_kmers; ++i) {
		offset_list[i] = offset_list[i] & OffsetMask;
		hbn_assert(offset_list[i] < max_offset);
	}
	
	hbn_timing_end(__func__);
}

u64 hash_extractor(void* list, const u64 i)
{
	u64* ul = list;
	u64 u = ul[i];
	return OFFSET_HASH(u);
}

u64 offset_extractor(void* list, const u64 i)
{
	u64* ul = list;
	u64 u = ul[i];
	return OFFSET_OFFSET(u);
}

void set_offset_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
	u64* su = src;
	u64* du = dst;
	du[dst_idx] = su[src_idx];
}

LookupTable*
build_lookup_table(const text_t* text, 
				   const int kmer_size, 
				   const int window_size, 
				   const int num_threads,
				   const double repeat_frac)
{
	size_t num_kmers = 0;
	u64* offset_list = get_offset_list(text, kmer_size, window_size, &num_kmers);
	radix_sort(offset_list, sizeof(u64), num_kmers, num_threads, offset_extractor, hash_extractor, set_offset_value);
	khash_t(KmerHash2OffsetMap)* hash_2_offset_map = get_kmer_counts(offset_list, num_kmers, kmer_size, repeat_frac);
	build_kmer_starts(offset_list, num_kmers, kmer_size, hash_2_offset_map);
	remove_repetitive_kmers(offset_list, num_kmers, repeat_frac, hash_2_offset_map);
	clear_hash_in_offset_list(offset_list, num_kmers, seqdb_max_offset(text));

	LookupTable* lktbl = (LookupTable*)malloc(sizeof(LookupTable));
	lktbl->offset_list = offset_list;
	lktbl->kmer_stats = hash_2_offset_map;
	return lktbl;
}

LookupTable*
destroy_lookup_table(LookupTable* lktbl)
{
	if (lktbl->offset_list) free(lktbl->offset_list);
	khash_t(KmerHash2OffsetMap)* hash_2_offset_map = (khash_t(KmerHash2OffsetMap)*)(lktbl->kmer_stats);
	kh_destroy(KmerHash2OffsetMap, hash_2_offset_map);
	free(lktbl);
	return 0;
}

u64*
extract_kmer_list(const LookupTable* lktbl, const u64 hash, u64* n)
{
	u64* list = 0;
	*n = 0;

	khash_t(KmerHash2OffsetMap)* hash_2_offset_map = (khash_t(KmerHash2OffsetMap)*)(lktbl->kmer_stats);
	khiter_t pos = kh_get(KmerHash2OffsetMap, hash_2_offset_map, hash);
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
