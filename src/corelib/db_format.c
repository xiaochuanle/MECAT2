#include "db_format.h"

#include "line_reader.h"

EDbFormat
hbn_guess_db_format(const char* path)
{
    EDbFormat fmt = eDbFormatEmptyFile;
    HbnLineReader* reader = HbnLineReaderNew(path);
    while (!HbnLineReaderAtEof(reader)) {
        HbnLineReaderReadOneLine(reader);
        kstring_t* line = &reader->line;
        truncate_end_spaces(line);
        // ignore lines containing only whitespace
        if (ks_empty(*line)) continue; 
        char c = ks_front(*line);
        // no content, just a comment or blank line
        if (c == '!' || c == '#' || c == ';') continue;
        if (c == '>') {
            fmt = eDbFormatFasta;
        } else if (c == '@') {
            fmt = eDbFormatFastq;
        } else {
            fmt = eDbFormatUnknown;
        }
        // we just test the first valid line (not empty line, not comment line)
        break;
    }
    HbnLineReaderFree(reader);
    return fmt;
}