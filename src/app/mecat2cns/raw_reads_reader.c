#include "raw_reads_reader.h"

#include "../../ncbi_blast/c_ncbi_blast_aux.h"

RawReadsReader*
RawReadsReaderNew(const char* db_dir, const char* db_title, const BOOL use_batch_mode)
{
    RawReadsReader* reader = (RawReadsReader*)calloc(1, sizeof(RawReadsReader));
    reader->db_dir = db_dir;
    reader->db_title = db_title;
    reader->dbinfo = seqdb_load_volume_info(db_dir, db_title, 0);
    hbn_assert(reader->dbinfo.hdr_offset_from == 0);
    hbn_assert(reader->dbinfo.seq_offset_from == 0);
    hbn_assert(reader->dbinfo.seq_start_id == 0);
    reader->seq_names = load_seq_headers(db_dir, db_title, 0, reader->dbinfo.hdr_offset_to);
    reader->seqinfo_array = load_seq_infos(db_dir, db_title, 0, reader->dbinfo.num_seqs);
    reader->raw_reads_offset_array = NULL;
    reader->use_batch_mode = use_batch_mode;
    if (!use_batch_mode) {
        reader->packed_seq = seqdb_load_pac(db_dir, db_title, 0, reader->dbinfo.seq_offset_to);
        reader->pac_stream = NULL;
    } else {
        reader->raw_reads_offset_array = (size_t*)malloc(sizeof(size_t) * reader->dbinfo.num_seqs);
        char path[HBN_MAX_PATH_LEN];
        make_packed_seq_path(db_dir, db_title, path);
        hbn_fopen(reader->pac_stream, path, "rb");
    }
    return reader;
}

RawReadsReader*
RawReadsReaderFree(RawReadsReader* reader)
{
    if (!reader) return NULL;
    if (reader->seq_names) sfree(reader->seq_names);
    if (reader->seqinfo_array) sfree(reader->seqinfo_array);
    if (reader->raw_reads_offset_array) sfree(reader->raw_reads_offset_array);
    if (reader->packed_seq) sfree(reader->packed_seq);
    if (reader->pac_stream) hbn_fclose(reader->pac_stream);
    sfree(reader);
    return NULL;
}

void
RawReadsReaderLoadRawReads(const HbnConsensusInitHit* cns_hit_array, 
    const size_t cns_hit_count,
    RawReadsReader* reader)
{
    if (!reader->use_batch_mode) return;
    if (reader->packed_seq) sfree(reader->packed_seq);
    hbn_assert(reader->raw_reads_offset_array);

    memset(reader->raw_reads_offset_array, 0, sizeof(size_t) * reader->dbinfo.num_seqs);
    size_t num_res = 0;
    int min_template_id = cns_hit_array[0].sid;
    int max_template_id = cns_hit_array[cns_hit_count-1].sid + 1;
    for (int i = min_template_id; i < max_template_id; ++i) {
        hbn_assert(reader->raw_reads_offset_array[i] == 0);
        reader->raw_reads_offset_array[i] = 1;
        size_t s = reader->seqinfo_array[i].seq_size;
        s = (s + 3) / 4;
        s *= 4;
        num_res += s;
    }
    for (size_t i = 0; i < cns_hit_count; ++i) {
        const HbnConsensusInitHit* hit = cns_hit_array + i;
        if (reader->raw_reads_offset_array[hit->qid]) continue;
        reader->raw_reads_offset_array[hit->qid] = 1;
        size_t s = reader->seqinfo_array[hit->qid].seq_size;
        s = (s + 3) / 4;
        s *= 4;
        num_res += s;
    }
    hbn_assert((num_res % 4) == 0);

    size_t num_bytes = num_res / 4;
    reader->packed_seq = (u8*)calloc(num_bytes, 1);
    size_t bytes_idx = 0;
    int loaded_seqs = 0;
    size_t loaded_res = 0;
    for (int i = 0; i < reader->dbinfo.num_seqs; ++i) {
        if (!reader->raw_reads_offset_array[i]) continue;
        size_t s = reader->seqinfo_array[i].seq_offset;
        hbn_assert((s % 4) == 0);
        s /= 4;
        size_t n = reader->seqinfo_array[i].seq_size;
        ++loaded_seqs;
        loaded_res += n;
        n = (n + 3) / 4;
        reader->raw_reads_offset_array[i] = bytes_idx * 4 + 1;
        fseek(reader->pac_stream, s, SEEK_SET);
        u8* p = reader->packed_seq + bytes_idx;
        hbn_fread(p, 1, n, reader->pac_stream);
        bytes_idx += n;
    }
    hbn_assert(bytes_idx == num_bytes);

    HBN_LOG("load %d sequences, %zu residues", loaded_seqs, loaded_res);
}

void
RawReadsReaderExtractRead(RawReadsReader* reader, int id, int strand, vec_u8* seqv)
{
    hbn_assert(id < reader->dbinfo.num_seqs);
    hbn_assert(strand == FWD || strand == REV);
    size_t res_from, res_to, res_cnt;
    if (!reader->use_batch_mode) {
        hbn_assert(reader->raw_reads_offset_array == NULL);
        res_from = reader->seqinfo_array[id].seq_offset;
        res_cnt = reader->seqinfo_array[id].seq_size;
        res_to = res_from + res_cnt;
    } else {
        hbn_assert(reader->raw_reads_offset_array[id] > 0, "id = %d", id);
        res_from = reader->raw_reads_offset_array[id] - 1;
        res_cnt = reader->seqinfo_array[id].seq_size;
        res_to = res_from + res_cnt;
    }
    hbn_assert((res_from % 4) == 0);

    kv_resize(u8, *seqv, res_cnt);
    int pos = 0;
    if (strand == FWD) {
        for (size_t i = res_from; i < res_to; ++i) {
            u8 c = _get_pac(reader->packed_seq, i);
            kv_A(*seqv, pos) = c;
            ++pos;
        }
    } else {
        size_t i = res_to;
        while (i > res_from) {
            --i;
            u8 c = _get_pac(reader->packed_seq, i);
            kv_A(*seqv, pos) = 3 - c;
            ++pos;
        }
    }
}