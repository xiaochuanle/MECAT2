#ifndef __HBN_FORMATH_H
#define __HBN_FORMATH_H

#include "hbn_hit.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    eOutFmtCan = 0,
    eOutFmtCanBin,
    eOutFmtM4,
    eOutFmtM4Bin,
    eOutFmtM4x,
    eOutFmtM4y,
    eOutFmtSam,
    eOutFmtError
} EOutputFormat;

#define outfmt_is_can(outfmt) ((outfmt)>=eOutFmtCan && (outfmt)<=eOutFmtCanBin)
#define outfmt_is_m4(outfmt) ((outfmt)>=eOutFmtM4 && (outfmt)<=eOutFmtM4y)
#define outfmt_is_sam(outfmt) ((outfmt)==eOutFmtSam)

const char*
outfmt_2_string(EOutputFormat outfmt);

EOutputFormat
string_2_outfmt(const char* str);

void
dump_hbn_align_results(const text_t* queries,
    const text_t* database,
    HbnAlignResults* results,
    EOutputFormat outfmt);

#ifdef __cplusplus
}
#endif

#endif // __HBN_FORMATH_H