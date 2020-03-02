#include "seqdb.h"

void
make_ambig_subseq_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "ambig");
}

CAmbigSubseq*
load_ambig_subseqs(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    const size_t bytes_from = sizeof(CAmbigSubseq) * from;
    const size_t bytes_to = sizeof(CAmbigSubseq) * to;
    char path[HBN_MAX_PATH_LEN];
    make_ambig_subseq_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(bytes_to <= total_bytes);
    size_t n = to - from;
    CAmbigSubseq* ambig_subseq_list = (CAmbigSubseq*)malloc(bytes_to - bytes_from);
    hbn_dfopen(in, path, "rb");
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(ambig_subseq_list, sizeof(CAmbigSubseq), n, in);
    hbn_fclose(in);
    return ambig_subseq_list;
}

void
make_header_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "hdr");
}

char* 
load_seq_headers(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    char path[HBN_MAX_PATH_LEN];
    make_header_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(to <= total_bytes, "to = %zu, total_bytes = %zu", to, total_bytes);
    const size_t n = to - from;
    char* hdr = (char*)malloc(n);
    hbn_dfopen(in, path, "rb");
    fseek(in, from, SEEK_SET);
    hbn_fread(hdr, 1, n, in);
    hbn_fclose(in);
    return hdr;
}

void
make_seq_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "seq_info");
}

void
make_seqdb_summary_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "summary");
}

CSeqInfo*
load_seq_infos(const char* data_dir, const char* db_name, const size_t from, const size_t to)
{
    hbn_assert(from >= 0);
    hbn_assert(from <= to);
    const size_t bytes_from = sizeof(CSeqInfo) * from;
    const size_t bytes_to = sizeof(CSeqInfo) * to;
    char path[HBN_MAX_PATH_LEN];
    make_seq_info_path(data_dir, db_name, path);
    const size_t total_bytes = hbn_file_size(path);
    hbn_assert(bytes_to <= total_bytes);
    size_t n = to - from;
    CSeqInfo* seq_info_list = (CSeqInfo*)malloc(sizeof(CSeqInfo) * n);
    hbn_dfopen(in, path, "rb");
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(seq_info_list, sizeof(CSeqInfo), n, in);
    hbn_fclose(in);
    return seq_info_list;
}

u8*
seqdb_load_pac(const char* seqdb_dir, const char* seqdb_title, const size_t res_from, size_t res_to)
{
    hbn_assert((res_from&3) == 0);
    hbn_assert(res_from <= res_to);
    const size_t bytes_from = res_from >> 2;
    const size_t bytes_to = (res_to + 3) >> 2;
    const size_t bytes_cnt = bytes_to - bytes_from;
    char path[HBN_MAX_PATH_LEN];
    make_packed_seq_path(seqdb_dir, seqdb_title, path);
    hbn_dfopen(in, path, "rb");
    u8* pac = (u8*)calloc(bytes_cnt, sizeof(u8));
    fseek(in, bytes_from, SEEK_SET);
    hbn_fread(pac, sizeof(u8), bytes_cnt, in);
    hbn_fclose(in);
    return pac;
}

u8*
seqdb_load_unpac(const char* seqdb_dir, const char* seqdb_title, const size_t res_from, size_t res_to)
{
    u8* pac = seqdb_load_pac(seqdb_dir, seqdb_title, res_from, res_to);
    size_t res_cnt = res_to - res_from;
    u8* unpac = (u8*)calloc(res_cnt, sizeof(u8));
    for (size_t i = 0; i < res_cnt; ++i) unpac[i] = _get_pac(pac, i);
    free(pac);
    return unpac;
}

void
make_ascii_volume_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "volume_info.txt");
}

void
make_bin_volume_info_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "volume_info.bin");
}

CSeqDBInfo
seqdb_load_volume_info(const char* data_dir, const char* db_name, int volume_index)
{
    char path[HBN_MAX_PATH_LEN];
    make_bin_volume_info_path(data_dir, db_name, path);
    const size_t file_size = hbn_file_size(path);
    hbn_assert(file_size % sizeof(CSeqDBInfo) == 0);
    const int nvol = file_size / sizeof(CSeqDBInfo);
    hbn_assert(volume_index >= 0, ": %d", volume_index);
    hbn_assert(volume_index < nvol);
    CSeqDBInfo* vol_info_list = (CSeqDBInfo*)malloc(file_size);
    hbn_dfopen(in, path, "rb");
    hbn_fread(vol_info_list, sizeof(CSeqDBInfo), nvol, in);
    hbn_fclose(in);
    CSeqDBInfo vol_info = vol_info_list[volume_index];
    free(vol_info_list);
    return vol_info;
}

void
make_packed_seq_path(const char* data_dir, const char* db_name, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (db_name) {
        strcat(path, db_name);
        strcat(path, ".");
    }
    strcat(path, "pac");
}

size_t seqdb_seq_offset(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs);
    return seqdb->seq_info_list[seq_id].seq_offset;
}

size_t seqdb_seq_size(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs, 
        "seq_id = %d, num_seqs = %d", seq_id, seqdb->dbinfo.num_seqs);
    return seqdb->seq_info_list[seq_id].seq_size;    
}

const char* seqdb_seq_name(const CSeqDB* seqdb, const int seq_id)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs);
    return seqdb->seq_header_list + seqdb->seq_info_list[seq_id].hdr_offset;
}

int seqdb_num_seqs(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.num_seqs;
}

size_t seqdb_size(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.db_size;
}

int seqdb_offset_to_seq_id(const CSeqDB* seqdb, const size_t offset)
{
    int left = 0, mid = 0, right = seqdb->dbinfo.num_seqs;
    while (left < right) {
        mid = (left + right) >> 1;
        if (offset >= seqdb->seq_info_list[mid].seq_offset) {
            if (mid == seqdb->dbinfo.num_seqs - 1) break;
            if (offset < seqdb->seq_info_list[mid+1].seq_offset) break;
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return mid;    
}

size_t seqdb_max_offset(const CSeqDB* seqdb)
{
    return seqdb->dbinfo.seq_offset_to - seqdb->dbinfo.seq_offset_from;
}

int seqdb_load_num_volumes(const char* seqdb_dir, const char* seqdb_title)
{
    char path[HBN_MAX_PATH_LEN];
    make_bin_volume_info_path(seqdb_dir, seqdb_title, path);
    size_t file_size = hbn_file_size(path);
    hbn_assert(file_size > 0);
    hbn_assert(file_size % sizeof(CSeqDBInfo) == 0);
    const int num_vols = file_size / sizeof(CSeqDBInfo);
    hbn_assert(num_vols > 1);
    return num_vols - 1;
}

int seqdb_load_num_reads(const char* seqdb_dir, const char* seqdb_title)
{
    CSeqDBInfo dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, 0);
    return dbinfo.num_seqs;
}

void
seqdb_extract_subsequence(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq)
{
    hbn_assert(seq_id < seqdb->dbinfo.num_seqs);
    size_t start = seqdb_seq_offset(seqdb, seq_id);
    size_t size = seqdb_seq_size(seqdb, seq_id);
    hbn_assert(from <= to && to <= size);
    kv_clear(*seq);
    if (strand == FWD) {
        size_t sfrom = start + from;
        size_t sto = start + to;
        for (size_t i = sfrom; i < sto; ++i) {
            u8 c = _get_pac(seqdb->packed_seq, i);
            kv_push(u8, *seq, c);
        }        
    } else {
        size_t sfrom = start + (size - to);
        size_t sto = start + (size - from);
        size_t i = sto;
        while (i > sfrom) {
            --i;
            u8 c = _get_pac(seqdb->packed_seq, i);
            c = 3 - c;
            kv_push(u8, *seq, c);
        }
    }
}

void
seqdb_extract_sequence(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq)
{
    const size_t seq_size = seqdb_seq_size(seqdb, seq_id);
    seqdb_extract_subsequence(seqdb, seq_id, 0, seq_size, strand, seq);
}

CSeqDB*
seqdb_load(const char* seqdb_dir, const char* seqdb_title, int vol_id)
{
    CSeqDB* vol = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    ++vol_id;
    vol->dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, vol_id);
    vol->packed_seq = seqdb_load_pac(seqdb_dir, seqdb_title, vol->dbinfo.seq_offset_from, vol->dbinfo.seq_offset_to);
    vol->seq_header_list = load_seq_headers(seqdb_dir, seqdb_title, vol->dbinfo.hdr_offset_from, vol->dbinfo.hdr_offset_to);
    vol->seq_info_list = load_seq_infos(seqdb_dir, seqdb_title, vol->dbinfo.seq_start_id, vol->dbinfo.seq_start_id + vol->dbinfo.num_seqs);
    vol->ambig_subseq_list = load_ambig_subseqs(seqdb_dir, seqdb_title, vol->dbinfo.ambig_offset_from, vol->dbinfo.ambig_offset_to);

    const size_t hdr_offset_from = vol->dbinfo.hdr_offset_from;
    const size_t seq_offset_from = vol->dbinfo.seq_offset_from;
    for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
        hbn_assert(vol->seq_info_list[i].hdr_offset >= hdr_offset_from);
        vol->seq_info_list[i].hdr_offset -= hdr_offset_from;
        hbn_assert(vol->seq_info_list[i].seq_offset >= seq_offset_from);
        vol->seq_info_list[i].seq_offset -= seq_offset_from;
    }

    return vol;
}

CSeqDB*
seqdb_load_unpacked(const char* seqdb_dir, const char* seqdb_title, int vol_id)
{
    CSeqDB* vol = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    ++vol_id;
    vol->dbinfo = seqdb_load_volume_info(seqdb_dir, seqdb_title, vol_id);
    vol->unpacked_seq = seqdb_load_unpac(seqdb_dir, seqdb_title, vol->dbinfo.seq_offset_from, vol->dbinfo.seq_offset_to);
    vol->seq_header_list = load_seq_headers(seqdb_dir, seqdb_title, vol->dbinfo.hdr_offset_from, vol->dbinfo.hdr_offset_to);
    vol->seq_info_list = load_seq_infos(seqdb_dir, seqdb_title, vol->dbinfo.seq_start_id, vol->dbinfo.seq_start_id + vol->dbinfo.num_seqs);
    vol->ambig_subseq_list = load_ambig_subseqs(seqdb_dir, seqdb_title, vol->dbinfo.ambig_offset_from, vol->dbinfo.ambig_offset_to);

    const size_t hdr_offset_from = vol->dbinfo.hdr_offset_from;
    const size_t seq_offset_from = vol->dbinfo.seq_offset_from;
    for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
        hbn_assert(vol->seq_info_list[i].hdr_offset >= hdr_offset_from);
        vol->seq_info_list[i].hdr_offset -= hdr_offset_from;
        hbn_assert(vol->seq_info_list[i].seq_offset >= seq_offset_from);
        vol->seq_info_list[i].seq_offset -= seq_offset_from;        
    }
    return vol;
}

CSeqDB*
CSeqDBFree(CSeqDB* vol)
{
    free(vol->packed_seq);
    free(vol->unpacked_seq);
    free(vol->seq_header_list);
    free(vol->seq_info_list);
    free(vol->ambig_subseq_list);
    free(vol);
    return NULL;
}

CSeqDB*
CSeqDBNew()
{
    CSeqDB* db = (CSeqDB*)calloc(1, sizeof(CSeqDB));
    return db;
}