#ifndef __SEQDB_H
#define __SEQDB_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define INIT_DB_DIR HBN_PACKAGE_NAME "_db"
#define INIT_QUERY_DB_TITLE "query"
#define INIT_SUBJECT_DB_TITLE "subject"

typedef struct {
    size_t offset;
    char ambig_residue;
    int count;
} CAmbigSubseq;

typedef struct {
    size_t seq_offset;
    size_t seq_size;
    size_t hdr_offset;
    size_t hdr_size;
    size_t ambig_offset;
    size_t ambig_size;
} CSeqInfo;

typedef struct {
    int seq_start_id;
    int num_seqs;
    size_t db_size;
    size_t seq_offset_from;
    size_t seq_offset_to;
    size_t hdr_offset_from;
    size_t hdr_offset_to;
    size_t ambig_offset_from;
    size_t ambig_offset_to;
} CSeqDBInfo;

#define dump_seqdb_info(output_func, stream, dbinfo) \
    output_func(stream, \
        "%d\t%d\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\t%zu\n", \
        (dbinfo).seq_start_id, \
        (dbinfo).num_seqs, \
        (dbinfo).db_size, \
        (dbinfo).seq_offset_from, \
        (dbinfo).seq_offset_to, \
        (dbinfo).hdr_offset_from, \
        (dbinfo).hdr_offset_to, \
        (dbinfo).ambig_offset_from, \
        (dbinfo).ambig_offset_to \
    )

typedef struct {
    CSeqDBInfo dbinfo;
    CSeqInfo* seq_info_list;
    char* seq_header_list;
    CAmbigSubseq* ambig_subseq_list;
    u8* packed_seq;
    u8* unpacked_seq;
    //char* raw_seq;
} CSeqDB;

typedef CSeqDB text_t;

CSeqDB*
CSeqDBFree(CSeqDB* vol);

CSeqDB*
CSeqDBNew();

CAmbigSubseq*
load_ambig_subseqs(const char* data_dir, const char* db_name, const size_t from, const size_t to);

CSeqDB*
seqdb_load(const char* seqdb_dir, const char* seqdb_title, int vol_id);

CSeqDB*
seqdb_load_unpacked(const char* seqdb_dir, const char* seqdb_title, int vol_id);

void
make_ambig_subseq_path(const char* data_dir, const char* db_name, char path[]);

void
make_ascii_volume_info_path(const char* data_dir, const char* db_name, char path[]);

void
make_bin_volume_info_path(const char* data_dir, const char* db_name, char path[]);

void
make_seq_info_path(const char* data_dir, const char* db_name, char path[]);

void
make_header_path(const char* data_dir, const char* db_name, char path[]);

char* 
load_seq_headers(const char* data_dir, const char* db_name, const size_t from, const size_t to);

void
make_packed_seq_path(const char* data_dir, const char* db_name, char path[]);

CSeqInfo*
load_seq_infos(const char* data_dir, const char* db_name, const size_t from, const size_t to);

size_t seqdb_seq_offset(const CSeqDB* seqdb, const int seq_id);

size_t seqdb_seq_size(const CSeqDB* seqdb, const int seq_id);

const char* seqdb_seq_name(const CSeqDB* seqdb, const int seq_id);

int seqdb_num_seqs(const CSeqDB* seqdb);

size_t seqdb_size(const CSeqDB* seqdb);

int seqdb_offset_to_seq_id(const CSeqDB* seqdb, const size_t offset);

size_t seqdb_max_offset(const CSeqDB* seqdb);

CSeqDBInfo
seqdb_load_volume_info(const char* data_dir, const char* db_name, int volume_index);

int seqdb_load_num_volumes(const char* seqdb_dir, const char* seqdb_title);

int seqdb_load_num_reads(const char* seqdb_dir, const char* seqdb_title);

u8*
seqdb_load_pac(const char* seqdb_dir, const char* seqdb_title, const size_t res_from, size_t res_to);

void
seqdb_extract_subsequence(const CSeqDB* seqdb,
    const int seq_id,
    const size_t from,
    const size_t to,
    const int strand,
    vec_u8* seq);

void
seqdb_extract_sequence(const CSeqDB* seqdb,
    const int seq_id,
    const int strand,
    vec_u8* seq);

#ifdef __cplusplus
}
#endif

#endif // __SEQDB_H