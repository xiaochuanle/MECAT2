#ifndef __CSTR_UTIL_H
#define __CSTR_UTIL_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

BOOL is_blank_string(const char* str);

int truncate_both_end_spaces(char* cstr);

const char* u64_to_string_datasize(u64 n, char buf[]);

const char* u64_to_string_comma(u64 n, char buf[]);

const char* double_to_string(double n, char buf[]);

u64 datasize_to_u64(const char* str);

int string_to_int(const char* str);

#ifdef __cplusplus
}
#endif

#endif // __CSTR_UTIL_H