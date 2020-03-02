#ifndef __SEQUENCE_H
#define __SEQUENCE_H

#include "line_reader.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    kstring_t sequence;
    kstring_t name;
    HbnLineReader* reader;
} HbnSequence;

#define hbn_seq_init(seq) (ks_init((seq).sequence), ks_init((seq).name), (seq).reader = NULL)

#define hbn_seq_dinit(name) HbnSequence name; hbn_seq_init(name)

#define hbn_seq_at_eof(seq) (hbn_line_reader_at_eof((seq).reader))

void hbn_sequence_read_one_seq(HbnSequence* seq);

void
hbn_sequence_reset(HbnSequence* seq);

HbnSequence*
hbn_sequence_free(HbnSequence* seq);

HbnSequence*
hbn_sequence_new(const char* path);

#ifdef __cplusplus
}
#endif

#endif // __SEQUENCE_H