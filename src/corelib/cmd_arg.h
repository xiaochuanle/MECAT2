#ifndef __CMD_ARG_H
#define __CMD_ARG_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    eCmdArgParseSuccess,
    eCmdArgParseExitNormally,
    eCmdArgParseError
} ECmdArgParseStatus;

ECmdArgParseStatus
validate_cmd_arg_cnt(const char* parser,
    const char* opt_name,
    const int argc,
    const int arg_idx,
    const int arg_cnt_expected);

ECmdArgParseStatus
parse_int_arg(const char* parser, const char* opt_name, char* cmd_arg, int64_t* ret);

ECmdArgParseStatus
parse_real_arg(const char* parser, const char* opt_name, char* cmd_arg, double* ret);

#ifdef __cplusplus
}
#endif

#endif // __CMD_ARG_H