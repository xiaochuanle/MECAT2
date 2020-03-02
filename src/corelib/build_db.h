#ifndef __BUILD_DB_H
#define __BUILD_DB_H

#include "seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

int 
hbndb_is_built(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq);

void
hbndb_make_built(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq);

void
build_db(const char* input,
    const char* seqdb_dir,
    const char* seqdb_title,
    const int min_seq_size,
    const int max_file_seqs,
    const size_t max_file_res,
    const int rename_seq);

#ifdef __cplusplus
}
#endif

#endif // __BUILD_DB_H