#ifndef __FASTA_H
#define __FASTA_H

#include "line_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* filename;
    HbnLineReader* line_reader;
    kstring_t name;
    kstring_t comment;
    kstring_t sequence;
    kstring_t plus;
    kstring_t qual;
    BOOL skip_error_formated_sequences;
    void* cached_names;
} HbnFastaReader;

HbnFastaReader*
HbnFastaReaderNew(const char* filename);

HbnFastaReader*
HbnFastaReaderFree(HbnFastaReader* reader);

kstring_t* HbnFastaReaderName(HbnFastaReader* reader);

kstring_t* HbnFastaReaderComment(HbnFastaReader* reader);

kstring_t* HbnFastaReaderSequence(HbnFastaReader* reader);

kstring_t* HbnFastaReaderPlusLine(HbnFastaReader* reader);

kstring_t* HbnFastaReaderQualityScoresLine(HbnFastaReader* reader);

void HbnFastaReaderSkipErrorFormatedSequences(HbnFastaReader* reader);

size_t HbnFastaReaderLineNumber(const HbnFastaReader* reader);

BOOL HbnFastaReaderReadOneSeq(HbnFastaReader* reader);

#ifdef __cplusplus
}
#endif

#endif // __FASTA_H