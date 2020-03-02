#ifndef __LINE_READER_H
#define __LINE_READER_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    gzFile      stream;
    BOOL        eof;
    BOOL        ungetline;
    size_t      last_read_size;
    size_t      buffer_size;
    char*       buffer;
    const char* pos;
    const char* end;
    kstring_t   line;
    size_t      input_pos;
    size_t      line_number;
} HbnBufferedLineReader;

typedef HbnBufferedLineReader HbnLineReader;

kstring_t* HbnLineReaderLine(HbnLineReader* reader);

size_t HbnLineReaderLineNumber(const HbnLineReader* reader);

size_t HbnLineReaderPosition(const HbnLineReader* reader);

void HbnLineReaderReadOneLine(HbnLineReader* reader);

void HbnLineReaderUngetline(HbnLineReader* reader);

char HbnLineReaderPeekChar(HbnLineReader* reader);

BOOL HbnLineReaderAtEof(const HbnLineReader* reader);

HbnLineReader*
HbnLineReaderFree(HbnLineReader* reader);

HbnLineReader*
HbnLineReaderNew(const char* filename);

#ifdef __cplusplus
}
#endif

#endif // __LINE_READER_H