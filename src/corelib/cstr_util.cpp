#include "cstr_util.h"

#include "../ncbi_blast/str_util/ncbistr.hpp"

using namespace std;

extern "C"
BOOL is_blank_string(const char* str)
{
    return NStr::IsBlank(str);
}

extern "C"
int truncate_both_end_spaces(char* cstr)
{
    string str(cstr);
    NStr::TruncateSpaces(str);
    memcpy(cstr, str.c_str(), str.size());
    cstr[str.size()] = '\0';
    return str.size();
}

extern "C"
const char* u64_to_string_datasize(u64 n, char buf[])
{
    string str = NStr::UInt8ToString_DataSize(n);
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = '\0';
    return buf;
}

extern "C"
const char* u64_to_string_comma(u64 n, char buf[])
{
    string str = NStr::UInt8ToString(n, NStr::fWithCommas, 10);
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = '\0';
    return buf;
}

extern "C"
const char* double_to_string(double n, char buf[])
{
    string str = NStr::DoubleToString(n);
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = '\0';
    return buf;
}

extern "C"
u64 datasize_to_u64(const char* str)
{
    return NStr::StringToUInt8_DataSize(str);
}

extern "C"
int string_to_int(const char* str)
{
    return NStr::StringToInt(str);
}