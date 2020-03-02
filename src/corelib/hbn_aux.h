#ifndef __HBN_AUX_H
#define __HBN_AUX_H

#include <stdarg.h>
#include <stdio.h>
#include <syslog.h>
#include <time.h>
#include <zlib.h>
#include <sys/time.h>

#include "kstring.h"
#include "hbn_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

#define __HBN_MACRO_ARG_CNT_AUX(\
	 _0,  _1,  _2,  _3,  _4,  _5,  _6,  _7,  _8,  _9, \
	_10, _11, _12, _13, _14, _15, _16, _17, _18, _19, \
	_20, _21, _22, _23, _24, _25, _26, _27, _28, _29, \
	_30, _31, _32, _33, _34, _35, _36, _37, _38, _39, \
	_40, _41, _42, _43, _44, _45, _46, _47, _48, _49, \
	_50, _51, _52, _53, _54, _55, _56, _57, _58, _59, \
	_60, _61, _62, _63, _64, N, ...) N

#define HBN_MACRO_ARG_CNT(...) __HBN_MACRO_ARG_CNT_AUX(0, ##__VA_ARGS__,\
	64, 63, 62, 61, 60, \
	59, 58, 57, 56, 55, 54, 53, 52, 51, 50, \
	49, 48, 47, 46, 45, 44, 43, 42, 41, 40, \
	39, 38, 37, 36, 35, 34, 33, 32, 31, 30, \
	29, 28, 27, 26, 25, 24, 23, 22, 21, 20, \
	19, 18, 17, 16, 15, 14, 13, 12, 11, 10, \
	 9,  8,  7,  6,  5,  4,  3,  2,  1,  0)

#define HBN_LOG_ARGS_DEFAULT    __FILE__, __FUNCTION__, __LINE__
#define HBN_LOG_ARGS_GENERIC    file, func, line
#define HBN_LOG_PARAMS_GENERIC  const char* file, const char* func, const int line

void
what_is_time_now(char now[]);

void
hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, ...);

#define HBN_LOG(fmt, args...) \
    hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_INFO, fmt, ##args)

#define HBN_WARN(fmt, args...) \
    hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_WARNING, fmt, ##args)

#define HBN_ERR(fmt, args...) \
    do { hbn_dump_message(HBN_LOG_ARGS_DEFAULT, LOG_ERR, fmt, ##args); abort(); } while(0)

#define HBN_ERR_GENERIC(fmt, args...) \
    do { hbn_dump_message(HBN_LOG_ARGS_GENERIC, LOG_ERR, fmt, ##args); abort(); } while(0)

/// gzFile

int safe_gzread(HBN_LOG_PARAMS_GENERIC, gzFile stream, void* buf, unsigned int len);
int err_gzread(gzFile stream, void* buf, unsigned int len);
gzFile safe_gzopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode);
int safe_gzclose(HBN_LOG_PARAMS_GENERIC, gzFile stream);

#define hbn_gzopen(stream, path, mode)      stream = safe_gzopen(HBN_LOG_ARGS_DEFAULT, path, mode)
#define hbn_dgzopen(stream, path, mode)     gzFile stream; hbn_gzopen(stream, path, mode)
#define hbn_gzclose(stream)                 safe_gzclose(HBN_LOG_ARGS_DEFAULT, stream)
#define hbn_gzread(stream, buf, len)        safe_gzread(HBN_LOG_ARGS_DEFAULT, stream, buf, len)

/// file

#define HBN_SCANF(input_func, stream, nread, fmt, ...) \
do { \
    int __nread = input_func(stream, fmt, __VA_ARGS__); \
    if (nread != __nread) { \
        HBN_ERR("scanf error: should read %d items, but only %d have been read", nread, __nread); \
    } \
} while(0)

FILE* safe_fopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode);
size_t safe_fwrite(HBN_LOG_PARAMS_GENERIC, const void* buf, size_t size, size_t nmemb, FILE* stream);
size_t safe_fread(HBN_LOG_PARAMS_GENERIC, void* buf, size_t size, size_t nmemb, FILE* stream);
int safe_fclose(HBN_LOG_PARAMS_GENERIC, FILE* stream);
off_t hbn_get_file_size(HBN_LOG_PARAMS_GENERIC, const char* path);

#define hbn_fopen(stream, path, mode)           stream = safe_fopen(HBN_LOG_ARGS_DEFAULT, path, mode)
#define hbn_dfopen(stream, path, mode)          FILE* stream; hbn_fopen(stream, path, mode)
#define hbn_fwrite(buf, size, nmemb, stream)    safe_fwrite(HBN_LOG_ARGS_DEFAULT, buf, size, nmemb, stream)
#define hbn_fread(buf, size, nmemb, stream)     safe_fread(HBN_LOG_ARGS_DEFAULT, buf, size, nmemb, stream)
#define hbn_fclose(stream)                      safe_fclose(HBN_LOG_ARGS_DEFAULT, stream)
#define hbn_file_size(path)                     hbn_get_file_size(HBN_LOG_ARGS_DEFAULT, path)

/// timing

double hbn_time_diff(const struct timeval* begin, const struct timeval* end);

#define hbn_timing_begin(title) \
    struct timeval hbn_begin; \
    gettimeofday(&hbn_begin, NULL); \
    HBN_LOG("'%s' BEGINS", title)

#define hbn_timing_end(title) \
    struct timeval hbn_end; \
    gettimeofday(&hbn_end, NULL); \
    double hbn_time_dur = hbn_time_diff(&hbn_begin, &hbn_end); \
    HBN_LOG("'%s' takes %.2lf secs.", title, hbn_time_dur)

/// assert

void
hbn_exception(const char* expr, HBN_LOG_PARAMS_GENERIC, const char* fmt, ...);

#define __hbn_assert(expr, ...) \
    do { \
        if (!(expr)) { \
            hbn_exception(#expr, __VA_ARGS__, NULL); \
            abort(); \
        } \
    } while(0)

#define hbn_assert(expr, args...) __hbn_assert(expr, HBN_LOG_ARGS_DEFAULT, ##args)

/// others

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

#define hbn_system(cmd) do { if(system(cmd) != 0) HBN_ERR("fail to running '%s'", cmd); } while(0)

#define hbn_min(a, b) ((a) < (b) ? (a) : (b))
#define hbn_max(a, b) ((a) > (b) ? (a) : (b))

extern u8 nst_nt4_table[256];
extern u8 nst_nt16_table[256];

/// sort

void ks_introsort_i32(size_t n, i32* a);
void ks_introsort_idx(size_t n, idx* a);
void ks_introsort_u32(size_t n, u32* a);
void ks_introsort_u64(size_t n, u64* a);
void ks_introsort_uidx(size_t n, uidx* a);
void ks_introsort_size_t(size_t n, size_t* a);
void ks_introsort_int_pair(size_t n, IntPair* a);

void i64_to_string_with_comma_r(i64 n, char str_n[]);

char* i64_to_string_with_comma(i64 n);

#define HBN_DIGIT_WIDTH 8

void u64_to_fixed_width_string_r(u64 n, char n_str[], const int width);

char* u64_to_fixed_width_string(u64 n, const int width);

void
print_digit_with_comma(FILE* out, size_t num);

void
print_fixed_width_string(FILE* out, const char* s, const int width);

int string_is_valid_number(const char * arg);

int hbn_get_cpu_count();

#ifdef __cplusplus
}
#endif

#endif // __HBN_AUX_H