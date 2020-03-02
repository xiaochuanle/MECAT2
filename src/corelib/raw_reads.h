#ifndef __RAW_READS_H
#define __RAW_READS_H

#include "gapped_candidate.h"
#include "seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    CSeqDB* raw_reads;
    int* id_maps;
    int max_global_id;
    int max_local_id;
} RawReads;

void
raw_reads_extract_sequence(const RawReads* raw_reads, const int gid, const int strand, vec_u8* seq);

int 
raw_reads_seq_size(const RawReads* raw_reads, const int gid);

const char*
raw_reads_seq_name(const RawReads* raw_reads, const int gid);

void
raw_reads_load(const char* seqdb_dir,
    const char* seqdb_title,
    const HbnConsensusInitHit* cns_hit_array,
    const size_t cns_hit_count,
    const int small_memory,
    RawReads* raw_reads);

#ifdef __cplusplus
}
#endif

#endif // __RAW_READS_H