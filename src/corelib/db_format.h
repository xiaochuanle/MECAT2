#ifndef __DB_FORMAT_H
#define __DB_FORMAT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    eDbFormatUnknown,
    eDbFormatFasta,
    eDbFormatFastq,
    eDbFormatEmptyFile
} EDbFormat;

EDbFormat
hbn_guess_db_format(const char* path);

#ifdef __cplusplus
}
#endif

#endif // __DB_FORMAT_H