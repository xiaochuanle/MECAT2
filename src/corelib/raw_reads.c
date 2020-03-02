#include "raw_reads.h"

static int
rr_set_raw_reads_to_be_loaded(const char* seqdb_dir,
    const char* seqdb_title,
    RawReads* raw_reads,
    const HbnConsensusInitHit* cns_hit_array,
    const size_t cns_hit_count)
{
    hbn_assert(raw_reads->id_maps == NULL);
    const int num_reads = seqdb_load_num_reads(seqdb_dir, seqdb_title);
    if (num_reads == 0) return 0;
    raw_reads->id_maps = (int*)calloc(num_reads, sizeof(int));
    for (size_t i = 0; i < cns_hit_count; ++i) {
        int qid = cns_hit_array[i].qid;
        int sid = cns_hit_array[i].sid;
        hbn_assert(qid < num_reads);
        hbn_assert(sid < num_reads);
        raw_reads->id_maps[qid] = 1;
        raw_reads->id_maps[sid] = 1;
    }
    int num_loaded_reads = 0;
    for (int i = 0; i < num_reads; ++i) num_loaded_reads += raw_reads->id_maps[i];
    return num_loaded_reads;
}

static size_t 
rr_calc_pac_res(const char* seqdb_dir, const char* seqdb_title, RawReads* raw_reads)
{
    const int num_reads = seqdb_load_num_reads(seqdb_dir, seqdb_title);
    CSeqInfo* seq_info_array = load_seq_infos(seqdb_dir, seqdb_title, 0, num_reads);
    size_t pac_res = 0;
    for (int i = 0; i < num_reads; ++i) {
        if (!raw_reads->id_maps[i]) continue;
        CSeqInfo si = seq_info_array[i];
        size_t s = ((si.seq_size + 3) >> 2) << 2;
        pac_res += s;
    }
    free(seq_info_array);
    return pac_res;
}

static size_t 
rr_calc_hdr_size(const char* seqdb_dir, const char* seqdb_title, RawReads* raw_reads)
{
    const int num_reads = seqdb_load_num_reads(seqdb_dir, seqdb_title);
    CSeqInfo* seq_info_array = load_seq_infos(seqdb_dir, seqdb_title, 0, num_reads);
    size_t hdr_size = 0;
    for (int i = 0; i < num_reads; ++i) {
        if (!raw_reads->id_maps[i]) continue;
        hdr_size += (seq_info_array[i].hdr_size + 1);
    }
    free(seq_info_array);
    return hdr_size;
}

static void
rr_load_one_volume(const char* seqdb_dir, 
    const char* seqdb_title,
    RawReads* raw_reads,
    int* next_lid,
    size_t* pac_bytes_offset,
    size_t* hdr_offset,
    const int vid)
{
    CSeqDB* vol = seqdb_load(seqdb_dir, seqdb_title, vid);
    CSeqDBInfo volinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, vid);
    size_t pbs = *pac_bytes_offset;
    size_t hbs = *hdr_offset;
    int gid = volinfo.seq_start_id;
    int lid = *next_lid;
    for (int i = 0; i < volinfo.num_seqs; ++i, ++gid) {
        if (!raw_reads->id_maps[gid]) continue;
        raw_reads->id_maps[gid] = (lid << 1) | 1;
        size_t s = seqdb_seq_size(vol, i);
        raw_reads->raw_reads->seq_info_list[lid].seq_size = s;
        raw_reads->raw_reads->seq_info_list[lid].seq_offset = pbs << 2;
        raw_reads->raw_reads->seq_info_list[lid].hdr_offset = hbs;
        raw_reads->raw_reads->seq_info_list[lid].hdr_size = vol->seq_info_list[i].hdr_size;

        size_t offset = seqdb_seq_offset(vol, i);
        hbn_assert((offset & 3) == 0);
        offset >>= 2;
        s = (s + 2) >> 2;
        memcpy(raw_reads->raw_reads->packed_seq + pbs, vol->packed_seq + offset, s);
        pbs += s;

        s = vol->seq_info_list[i].hdr_size;
        const char* hdr = seqdb_seq_name(vol, i);
        memcpy(raw_reads->raw_reads->seq_header_list + hbs, hdr, s + 1);
        hbs += (s + 1);

        ++lid;
    }
    CSeqDBFree(vol);
    *next_lid = lid;
    *pac_bytes_offset = pbs;
    *hdr_offset = hbs;
}

static void
rr_mem_alloc(RawReads* raw_reads,
    const int num_raw_reads_to_be_loaded,
    const size_t pac_residues,
    const size_t hdr_size)
{
    raw_reads->raw_reads = CSeqDBNew();
    raw_reads->raw_reads->seq_info_list = (CSeqInfo*)calloc(num_raw_reads_to_be_loaded, sizeof(CSeqInfo));
    raw_reads->raw_reads->seq_header_list = (char*)calloc(hdr_size, sizeof(char));
    hbn_assert((pac_residues & 3) == 0);
    raw_reads->raw_reads->packed_seq = (u8*)calloc(pac_residues>>2, sizeof(u8));

    CSeqDBInfo* dbinfo = &raw_reads->raw_reads->dbinfo;
    dbinfo->db_size = pac_residues;
    dbinfo->hdr_offset_from = 0;
    dbinfo->hdr_offset_to = hdr_size;
    dbinfo->num_seqs = num_raw_reads_to_be_loaded;
    dbinfo->seq_offset_from = 0;
    dbinfo->seq_offset_to = pac_residues;
    dbinfo->seq_start_id = 0;
}

static void
rr_load_raw_reads(const char* seqdb_dir,
    const char* seqdb_title,
    RawReads* raw_reads,
    const HbnConsensusInitHit* cns_hit_array,
    const size_t cns_hit_count)
{
    if (raw_reads->id_maps) {
        free(raw_reads->id_maps);
        raw_reads->id_maps = NULL;
    }
    if (raw_reads->raw_reads) {
        raw_reads->raw_reads = CSeqDBFree(raw_reads->raw_reads);
    }
    raw_reads->max_global_id = seqdb_load_num_reads(seqdb_dir, seqdb_title);
    raw_reads->max_local_id = rr_set_raw_reads_to_be_loaded(seqdb_dir, seqdb_title, raw_reads, cns_hit_array, cns_hit_count);
    const size_t pac_res = rr_calc_pac_res(seqdb_dir, seqdb_title, raw_reads);
    const size_t hdr_size = rr_calc_hdr_size(seqdb_dir, seqdb_title, raw_reads);
    rr_mem_alloc(raw_reads, raw_reads->max_local_id, pac_res, hdr_size);

    int local_id = 0;
    size_t pac_bytes_idx = 0;
    size_t hdr_idx = 0;
    const int num_vols = seqdb_load_num_volumes(seqdb_dir, seqdb_title);
    for (int i = 0; i < num_vols; ++i) {
        rr_load_one_volume(seqdb_dir, seqdb_title, raw_reads, &local_id, &pac_bytes_idx, &hdr_idx, i);
    }
    hbn_assert(local_id == raw_reads->max_local_id);
    hbn_assert(pac_bytes_idx * 4 == pac_res);
    hbn_assert(hdr_idx == hdr_size);

    HBN_LOG("load %s raw reads (%s residues)", i64_to_string_with_comma(raw_reads->max_local_id), i64_to_string_with_comma(pac_res));
}

void
raw_reads_load(const char* seqdb_dir,
    const char* seqdb_title,
    const HbnConsensusInitHit* cns_hit_array,
    const size_t cns_hit_count,
    const int small_memory,
    RawReads* raw_reads)
{
    if (!small_memory) {
        if (raw_reads->raw_reads) {
            hbn_assert(raw_reads->id_maps);
            return;
        }
        hbn_assert(raw_reads->id_maps == NULL);
        const int num_reads = seqdb_load_num_reads(seqdb_dir, seqdb_title);
        raw_reads->id_maps = (int*)calloc(num_reads, sizeof(int));
        for (int i = 0; i < num_reads; ++i) {
            raw_reads->id_maps[i] = (i << 1) | 1;
            hbn_assert(raw_reads->id_maps[i] & 1);
            hbn_assert((raw_reads->id_maps[i]>>1) == i);
        }
        raw_reads->raw_reads = seqdb_load(seqdb_dir, seqdb_title, -1);
        raw_reads->max_global_id = num_reads;
        raw_reads->max_local_id = num_reads;
        return;
    }
    rr_load_raw_reads(seqdb_dir, seqdb_title, raw_reads, cns_hit_array, cns_hit_count);
}

const char*
raw_reads_seq_name(const RawReads* raw_reads, const int gid)
{
    hbn_assert(gid < raw_reads->max_global_id);
    hbn_assert(raw_reads->id_maps[gid] & 1);
    int lid = raw_reads->id_maps[gid] >> 1;
    hbn_assert(lid < raw_reads->max_local_id);
    hbn_assert(lid < raw_reads->raw_reads->dbinfo.num_seqs);
    return raw_reads->raw_reads->seq_header_list + raw_reads->raw_reads->seq_info_list[lid].hdr_offset;
}

int 
raw_reads_seq_size(const RawReads* raw_reads, const int gid)
{
    hbn_assert(gid < raw_reads->max_global_id);
    hbn_assert(raw_reads->id_maps[gid] & 1);
    int lid = raw_reads->id_maps[gid] >> 1;
    hbn_assert(lid < raw_reads->max_local_id);
    hbn_assert(lid < raw_reads->raw_reads->dbinfo.num_seqs);
    return raw_reads->raw_reads->seq_info_list[lid].seq_size;    
}

void
raw_reads_extract_sequence(const RawReads* raw_reads, const int gid, const int strand, vec_u8* seq)
{
    hbn_assert(gid < raw_reads->max_global_id);
    hbn_assert(raw_reads->id_maps[gid] & 1);
    int lid = raw_reads->id_maps[gid] >> 1;
    hbn_assert(lid < raw_reads->max_local_id);
    hbn_assert(lid < raw_reads->raw_reads->dbinfo.num_seqs);
    seqdb_extract_sequence(raw_reads->raw_reads, lid, strand, seq);    
}