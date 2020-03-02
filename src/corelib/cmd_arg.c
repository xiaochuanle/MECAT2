#include "cmd_arg.h"

#include <stdio.h>

#include "hbn_aux.h"

ECmdArgParseStatus
validate_cmd_arg_cnt(const char* parser,
    const char* opt_name,
    const int argc,
    const int arg_idx,
    const int arg_cnt_expected)
{
    if (arg_idx + arg_cnt_expected > argc) {
        fprintf(stderr, "[%s] ERROR: parameters to option '%s' are not enough (%d expected)\n",
            parser, opt_name, arg_cnt_expected);
        return eCmdArgParseError;
    }
    return eCmdArgParseSuccess;
}

ECmdArgParseStatus
parse_int_arg(const char* parser, const char* opt_name, char* cmd_arg, int64_t* ret)
{
    if (!string_is_valid_number(cmd_arg)) {
        if (parser) {
            fprintf(stderr, "[%s] ", parser);
        }
        fprintf(stderr, "ERROR: argument %s", cmd_arg);
        if (opt_name) {
            fprintf(stderr, " to option '%s'", opt_name);
        }
        fprintf(stderr, " seems not a plausible integer.\n");
        return eCmdArgParseError;
    }
    
    char* endp = NULL;
    int64_t i = strtoll(cmd_arg, &endp, 0);
    hbn_assert(endp > cmd_arg);
    *ret = i;
    return eCmdArgParseSuccess;
}

ECmdArgParseStatus
parse_real_arg(const char* parser, const char* opt_name, char* cmd_arg, double* ret)
{
    if (!string_is_valid_number(cmd_arg)) {
        if (parser) {
            fprintf(stderr, "[%s] ", parser);
        }
        fprintf(stderr, "ERROR: argument %s", cmd_arg);
        if (opt_name) {
            fprintf(stderr, " to option '%s'", opt_name);
        }
        fprintf(stderr, " seems not a plausible real number.\n");
        return eCmdArgParseError;
    }
    
    char* endp = NULL;
    double e = strtod(cmd_arg, &endp);
    hbn_assert(endp > cmd_arg);
    *ret = e;
    return eCmdArgParseSuccess; 
}