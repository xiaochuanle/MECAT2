#ifndef __HBN_DEFS_H
#define __HBN_DEFS_H

#include <stdint.h>
#include <inttypes.h>

#include "kvec.h"

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef int8_t      i8;
typedef uint8_t     u8;
typedef int16_t     i16;
typedef uint16_t    u16;
typedef int32_t     i32;
typedef uint32_t    u32;
typedef int64_t     i64;
typedef uint64_t    u64;

typedef i8          Int1;
typedef u8          Uint1;
typedef i16         Int2;
typedef u16         Uint2;
typedef i32         Int4;
typedef u32         Uint4;
typedef i64         Int8;
typedef u64         Uint8;

#ifndef Boolean
typedef u8          Boolean;
#endif

#ifndef BOOL
typedef u8          BOOL;
#endif

#ifndef TRUE
#define TRUE        1
#define FALSE       0
#endif

typedef i64         idx;
typedef u64         uidx;
#define PRIdx       PRId64
#define PRUidx      PRIu64

#define I8_MAX      INT8_MAX
#define U8_MAX      UINT8_MAX
#define I16_MIN     INT16_MIN
#define I16_MAX     INT16_MAX
#define U16_MAX     UINT16_MAX
#define I32_MIN     INT32_MIN
#define I32_MAX     INT32_MAX
#define I64_MIN     INT64_MIN
#define I64_MAX     INT64_MAX
#define U64_MAX     UINT64_MAX
#define IDX_MAX     INT64_MAX
#define UIDX_MAX    UINT64_MAX

#define U8_ONE      ((u8)1)
#define U32_ONE		((u32)1)
#define U64_ONE		((u64)1)
#define UIDX_ONE    U64_ONE

#define FWD         (0)
#define REV         (1)
#define F_R         (2)
#define CMP_STRAND(__s)   (1-(__s))

#define MATCH_REWARD        (8)
#define MISMATCH_PENALTY    (-1)
#define AMB_PENALTY         (-1)
#define GAP_OPEN            (0)
#define GAP_EXTEND          (1)
#define KSW_ZDROP           (30)
#define KSW_BAND_WIDTH      (100)

#define GAP_CHAR            ('-')
#define GAP_CODE            (4)
#define DECODE_RESIDUE(__r) ("ACGT-"[(u64)(__r)])

typedef kvec_t(char)    vec_char;
typedef kvec_t(idx)     vec_idx;
typedef kvec_t(uidx)    vec_uidx;
typedef kvec_t(int)     vec_int;
typedef kvec_t(u8)      vec_u8;
typedef kvec_t(u64)     vec_u64;
typedef kvec_t(size_t)  vec_size_t;
typedef kvec_t(void*)   vec_void_ptr;

typedef struct {
    int first;
    int second;
} IntPair;

typedef kvec_t(IntPair) vec_int_pair;

#define HBN_MAX_PATH_LEN    2000

#ifdef __cplusplus
}
#endif

#endif // __HBN_DEFS_H