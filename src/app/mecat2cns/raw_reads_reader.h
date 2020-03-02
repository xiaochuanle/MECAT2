#ifndef __RAW_READS_READER_H
#define __RAW_READS_READER_H

#include "../../corelib/seqdb.h"
#include "../../corelib/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* db_dir;
    const char* db_title;
    CSeqDBInfo dbinfo;
    const char* seq_names;
    CSeqInfo* seqinfo_array;
    FILE* pac_stream;
    u8* packed_seq;
    BOOL use_batch_mode;
    size_t* raw_reads_offset_array;
} RawReadsReader;

void
RawReadsReaderExtractRead(RawReadsReader* reader, int id, int strand, vec_u8* seqv);

void
RawReadsReaderLoadRawReads(const HbnConsensusInitHit* cns_hit_array,
    const size_t cns_hit_count,
    RawReadsReader* reader);

RawReadsReader*
RawReadsReaderFree(RawReadsReader* reader);

RawReadsReader*
RawReadsReaderNew(const char* db_dir, const char* db_title, const BOOL use_batch_mode);

#ifdef __cplusplus
}
#endif

#endif // __RAW_READS_READER_H