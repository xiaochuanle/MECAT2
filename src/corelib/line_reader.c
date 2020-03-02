#include "line_reader.h"

/// Result codes for I/O operations.
/// @note
///     Exceptions (if any) thrown by IReader/IWriter interfaces should be
///     treated as unrecoverable errors (eRW_Error).
/// @sa
///   IReader, IWriter, IReaderWriter
typedef enum {
    eRW_NotImplemented = -1,  ///< Action / information is not available
    eRW_Success = 0,          ///< Everything is okay, I/O completed
    eRW_Timeout,              ///< Timeout expired, try again later
    eRW_Error,                ///< Unrecoverable error, no retry possible
    eRW_Eof                   ///< End of data, should be considered permanent
} ERW_Result;

HbnBufferedLineReader*
HbnBufferedLineReaderNew(const char* filename)
{
    HbnBufferedLineReader* reader = (HbnBufferedLineReader*)calloc(1, sizeof(HbnBufferedLineReader));
    hbn_gzopen(reader->stream, filename, "r");
    reader->eof = FALSE;
    reader->ungetline = FALSE;
    reader->buffer_size = 32 * 1024;
    reader->buffer = (char*)malloc(reader->buffer_size);
    reader->pos = reader->buffer;
    reader->end = reader->pos;
    reader->input_pos = 0;
    ks_init(reader->line);
    reader->line_number = 0;

    return reader;
}

HbnBufferedLineReader*
HbnBufferedLineReaderFree(HbnBufferedLineReader* reader)
{
    hbn_gzclose(reader->stream);
    free(reader->buffer);
    ks_destroy(reader->line);
    free(reader);
    return NULL;
}

BOOL HbnBufferedLineReaderAtEof(const HbnBufferedLineReader* reader)
{
    return gzeof(reader->stream) && (reader->pos >= reader->end) && (!reader->ungetline);
}

char HbnBufferedLineReaderPeekChar(HbnBufferedLineReader* reader)
{
    /* If at EOF - undefined behavior */
    if (HbnBufferedLineReaderAtEof(reader)) {
        return *reader->pos;
    }
    /* If line was ungot - return its first symbol */
    if (reader->ungetline) {
        /* if line is empty - return 0 */
        if (ks_empty(reader->line)) {
            return 0;
        }
        return ks_front(reader->line);
    }
    /* If line is empty - return 0 */
    if ( (*reader->pos == '\n') || (*reader->pos == '\r') ) {
        return 0;
    }
    return *reader->pos;
}

void HbnBufferedLineReaderUngetline(HbnBufferedLineReader* reader)
{
    hbn_assert((!reader->ungetline) && ks_s(reader->line));
    /* If after Ungetline() or after constructor - noop */
    if (reader->ungetline || ks_s(reader->line) == NULL) return;
    --reader->line_number;
    reader->ungetline = TRUE;
}

static ERW_Result
HbnBufferedLineReaderLoadData(gzFile stream, char* buffer, int count, int* bytes_read)
{
    ERW_Result result = eRW_Success;
    int n = gzread(stream, buffer, count);
    if (n < count) {
        hbn_assert(gzeof(stream));
        if (gzeof(stream)) {
            result = eRW_Eof;
        } else {
            int errno;
            const char* error_string = gzerror(stream, &errno);
            if (errno) {
                HBN_ERR("%s", error_string);
            }
        }
    }
    *bytes_read = n;
    return result;
}

static BOOL HbnBufferedLineReaderReadBuffer(HbnBufferedLineReader* reader)
{
    if (HbnLineReaderAtEof(reader)) return FALSE;

    reader->input_pos += (reader->end - reader->buffer);
    reader->pos = reader->end = reader->buffer;
    while (1) {
        int size;
        ERW_Result result = HbnBufferedLineReaderLoadData(reader->stream, 
                                reader->buffer, 
                                reader->buffer_size, 
                                &size);
        switch (result) {
        case eRW_NotImplemented:
        case eRW_Error:
            HBN_ERR("Read error");
            /* NOT REACHED */
            break;
        case eRW_Timeout:
            // keep spinning around
            break;
        case eRW_Eof:
            reader->eof = TRUE;
            /* fall through */
        case eRW_Success:
            reader->end = reader->pos + size;
            return (result == eRW_Success || size > 0);
        default:
            hbn_assert(0);
            break;
        }
    } // while
    return FALSE;
}

static void HbnBufferedLineReaderLoadLong(HbnBufferedLineReader* reader)
{
    const char* start = reader->pos;
    const char* end = reader->end;
    kputsn(start, end - start, &reader->line);
    while (HbnBufferedLineReaderReadBuffer(reader)) {
        start = reader->pos;
        end = reader->end;
        for (const char* p = start; p < end; ++p) {
            char c = *p;
            if (c == '\r' || c == '\n') {
                kputsn(start, p - start, &reader->line);
                reader->last_read_size = ks_size(reader->line) + 1;
                if (++p == end) {
                    if (HbnBufferedLineReaderReadBuffer(reader)) {
                        p = reader->pos;
                        end = reader->end;
                        if (p < end && c == '\r' && *p == '\n') {
                            ++p;
                            reader->pos = p;
                            ++reader->last_read_size;
                        }
                    }
                } else {
                    if (c == '\r' && *p == '\n') {
                        if (++p == end) {
                            HbnBufferedLineReaderReadBuffer(reader);
                            p = reader->pos;
                            ++reader->last_read_size;
                        }
                    }
                    reader->pos = p;
                }
                return;
            }
        }
        kputsn(start, end - start, &reader->line);
    }
    reader->last_read_size = ks_size(reader->line);
    return;
}

void HbnBufferedLineReaderReadOneLine(HbnBufferedLineReader* reader)
{
    /* If at EOF - noop */
    if (HbnBufferedLineReaderAtEof(reader)) {
        ks_destroy(reader->line);
        return;
    }
    ++reader->line_number;
    if (reader->ungetline) {
        hbn_assert(ks_s(reader->line));
        reader->ungetline = FALSE;
        return;
    }
    // check if we are at the buffer end
    ks_clear(reader->line);
    const char* start = reader->pos;
    const char* end = reader->end;
    for (const char* p = start; p < end; ++p) {
        if (*p == '\n') {
            kputsn(start, p - start, &reader->line);
            reader->last_read_size = p + 1 - start;
            reader->pos = ++p;
            if (p == end) {
                HbnBufferedLineReaderReadBuffer(reader);
            }
            return;
        } else if (*p == '\r') {
            kputsn(start, p - start, &reader->line);
            reader->last_read_size = p + 1 - start;
            reader->pos = ++p;
            if (p == end) {
                if (HbnBufferedLineReaderReadBuffer(reader)) {
                    p = reader->pos;
                    if (*p == '\n') {
                        reader->pos = p + 1;
                        ++reader->last_read_size;
                    }
                }
                return;
            }
            if (*p != '\n') return;
            ++reader->last_read_size;
            reader->pos = ++p;
            if (p == end) {
                HbnBufferedLineReaderReadBuffer(reader);
            }
            return;
        }
    }
    HbnBufferedLineReaderLoadLong(reader);
    return;
}

size_t HbnBufferedLineReaderPosition(const HbnBufferedLineReader* reader)
{
    size_t offset = reader->pos - reader->buffer;
    if (reader->ungetline) {
        offset -= reader->last_read_size;
    }
    return reader->input_pos + offset;
}

size_t HbnBufferedLineReaderLineNumber(const HbnBufferedLineReader* reader)
{
    /* Right after constructor (m_LineNumber is 0 and UngetLine() was not run)
       - 0 */
    /* If at EOF - returns the number of the last string */
    /* After UngetLine() - number of the previous string */
    /* Not at EOF, not after UngetLine() - number of the current string */
    return reader->line_number;
}

kstring_t* HbnBufferedLineReaderLine(HbnBufferedLineReader* reader)
{
    return &reader->line;
}

/// abstract type

kstring_t* HbnLineReaderLine(HbnLineReader* reader)
{
    return HbnBufferedLineReaderLine(reader);
}

size_t HbnLineReaderLineNumber(const HbnLineReader* reader)
{
    return HbnBufferedLineReaderLineNumber(reader);
}

size_t HbnLineReaderPosition(const HbnLineReader* reader)
{
    return HbnBufferedLineReaderPosition(reader);
}

void HbnLineReaderReadOneLine(HbnLineReader* reader)
{
    HbnBufferedLineReaderReadOneLine(reader);
}

void HbnLineReaderUngetline(HbnLineReader* reader)
{
    HbnBufferedLineReaderUngetline(reader);
}

char HbnLineReaderPeekChar(HbnLineReader* reader)
{
    return HbnBufferedLineReaderPeekChar(reader);
}

BOOL HbnLineReaderAtEof(const HbnLineReader* reader)
{
    return HbnBufferedLineReaderAtEof(reader);
}

HbnLineReader*
HbnLineReaderFree(HbnLineReader* reader)
{
    return HbnBufferedLineReaderFree(reader);
}

HbnLineReader*
HbnLineReaderNew(const char* filename)
{
    return HbnBufferedLineReaderNew(filename);
}