#ifndef __NUMERIC_STR_INTERCONV_HPP
#define __NUMERIC_STR_INTERCONV_HPP

#include "../ncbi_blast_aux.hpp"
#include "tempstr.hpp"
#include "../../corelib/hbn_aux.h"

BEGIN_HBNSTR_SCOPE

    using namespace ncbi;

    enum class EHbnStrConvErrorMode:int {
        eHbnStrConvErrorAbort = (1 << 0),
        eHbnStrConvErrorSetErrno = (1 << 2)
    };
    EHbnStrConvErrorMode GetStrConvErrorMode();
    void SetStrConvErrorMode(EHbnStrConvErrorMode mode);

    /// Number to string conversion flags.
    ///
    /// NOTE: 
    ///   If specified base in the *ToString() methods is not default 10,
    ///   that some flags like fWithSign and fWithCommas will be ignored.
    enum ENumToStringFlags {
        fWithSign                = (1 <<  6),  ///< Prefix the output value with a sign
        fWithCommas              = (1 <<  7),  ///< Use commas as thousands separator
        fDoubleFixed             = (1 <<  8),  ///< Use n.nnnn format for double
        fDoubleScientific        = (1 <<  9),  ///< Use scientific format for double
        fDoublePosix             = (1 << 10),  ///< Use C locale
        fDoubleGeneral           = fDoubleFixed | fDoubleScientific,
        // Additional flags to convert "software" qualifiers (see UInt8ToString_DataSize)
        fDS_Binary               = (1 << 11),
        fDS_NoDecimalPoint       = (1 << 12),
        fDS_PutSpaceBeforeSuffix = (1 << 13),
        fDS_ShortSuffix          = (1 << 14),
        fDS_PutBSuffixToo        = (1 << 15)
    };
    typedef int TNumToStringFlags;    ///< Bitwise OR of "ENumToStringFlags"

    /// String to number conversion flags.
    enum EStringToNumFlags {
        fMandatorySign           = (1 << 17),  ///< See 'ENumToStringFlags::fWithSign'
        fAllowCommas             = (1 << 18),  ///< See 'ENumToStringFlags::fWithCommas'
        fAllowLeadingSpaces      = (1 << 19),  ///< Can have leading spaces
        fAllowLeadingSymbols     = (1 << 20) | fAllowLeadingSpaces,
                                               ///< Can have leading non-nums
        fAllowTrailingSpaces     = (1 << 21),  ///< Can have trailing spaces
        fAllowTrailingSymbols    = (1 << 22) | fAllowTrailingSpaces,
                                               ///< Can have trailing non-nums
        fDecimalPosix            = (1 << 23),  ///< For decimal point, use C locale
        fDecimalPosixOrLocal     = (1 << 24),  ///< For decimal point, try both C and current locale
        fDecimalPosixFinite      = (1 << 25),  ///< Keep result finite and normalized:
                                               ///< if DBL_MAX < result < INF,     result becomes DBL_MAX
                                               ///< if       0 < result < DBL_MIN, result becomes DBL_MIN
        // Additional flags to convert "software" qualifiers (see StringToUInt8_DataSize)
        fDS_ForceBinary          = (1 << 26),
        fDS_ProhibitFractions    = (1 << 27),
        fDS_ProhibitSpaceBeforeSuffix = (1 << 28)
    };
    typedef int TStringToNumFlags;   ///< Bitwise OR of "EStringToNumFlags"

    /// Convert string to int.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    int StringToInt(const CTempString str,
                           TStringToNumFlags flags = 0,
                           int               base  = 10);

    /// Convert string to unsigned int.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    unsigned int StringToUInt(const CTempString str,
                                     TStringToNumFlags flags = 0,
                                     int               base  = 10);  

    /// Convert string to long.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    long StringToLong(const CTempString str,
                             TStringToNumFlags flags = 0,
                             int               base  = 10); 

    /// Convert string to unsigned long.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    unsigned long StringToULong(const CTempString str,
                                       TStringToNumFlags flags = 0,
                                       int               base  = 10); 

    /// Convert string to double-precision value (analog of strtod function)
    ///
    /// @param str
    ///   String to be converted.
    /// @param endptr
    ///   Pointer to character that stops scan.
    /// @return
    ///   Double-precision value.
    ///   This function always uses dot as decimal separator.
    ///   - on overflow, it returns HUGE_VAL and sets errno to ERANGE;
    ///   - on underflow, it returns 0 and sets errno to ERANGE;
    ///   - if conversion was impossible, it returns 0 and sets errno.
    ///   Also, when input string equals (case-insensitive) to
    ///   - "NAN", the function returns NaN;
    ///   - "INF" or "INFINITY", the function returns HUGE_VAL;
    ///   - "-INF" or "-INFINITY", the function returns -HUGE_VAL;
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, set errno to non-zero and return zero.
    ///   - Denormal or infinite results are considered successful conversion.
    ///   - To enforce finite and normalized result, use fDecimalPosixFinite flag.
    ///   - This function is meant to be more "low-level" than other
    ///     StringToXxx functions - for example, it allows trailing characters
    ///     (and doesn't include a flags parameter for tweaking such behavior).
    ///     This could result in strings like "nanosecond" being converted to
    ///     NaN, "-inf=input_file" being converted to -INF, or other unexpected
    ///     behavior. Therefore, please consider using StringToDouble unless
    ///     you specifically need this functionality.
    double StringToDoublePosix(const char* str, char** endptr=0,
                                      TStringToNumFlags flags=0);

    /// Convert string to double.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    ///   Do not support fAllowCommas flag.
    /// @return
    ///   - If invalid flags are passed, throw an exception.
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    /// @note
    ///   - Denormal or infinite results are considered successful conversion.
    ///   - To enforce finite and normalized result, use fDecimalPosixFinite flag.
    double StringToDouble(const CTempStringEx str,
                                 TStringToNumFlags   flags = 0);

    /// This version accepts zero-terminated string
    /// @deprecated
    ///   It is unsafe to use this method directly, please use StringToDouble()
    ///   instead.
    NCBI_DEPRECATED
    double StringToDoubleEx(const char* str, size_t size,
                                   TStringToNumFlags flags = 0);  

    /// Convert string to Int8.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    Int8 StringToInt8(const CTempString str,
                             TStringToNumFlags flags = 0,
                             int               base  = 10);  

    /// Convert string to Uint8.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Radix base. Allowed values are 0, 2..36. Zero means to use the
    ///   first characters to determine the base - a leading "0x" or "0X"
    ///   means base 16; otherwise a leading 0 means base 8; otherwise base 10.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    Uint8 StringToUInt8(const CTempString str,
                               TStringToNumFlags flags = 0,
                               int               base  = 10);  

    /// Convert string that can contain "software" qualifiers to Uint8. 
    ///
    /// String can contain "software" qualifiers: G(giga-), MB(mega-),
    /// KiB (kibi-) etc.
    /// Example: 100MB, 1024KiB, 5.7G.
    /// Meaning of qualifiers depends on flags and by default is 1000-based
    /// (i.e. K=1000, M=10^6 etc.) except in cases when qualifiers with "iB"
    /// are used, i.e. KiB=1024, MiB=1024^2 etc. When flags parameter contains
    /// fDS_ForceBinary then qualifiers without "iB" (i.e. "K" or "MB") will
    /// also be 1024-based.
    /// String can contain a decimal fraction (except when fDS_ProhibitFractions
    /// flag is used), in this case the resultant Uint8 number will be rounded
    /// to fit into integer value.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @return
    ///   - If invalid flags are passed, throw an exception.
    ///   - If conversion succeeds, return the converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    Uint8 StringToUInt8_DataSize(const CTempString str,
                                        TStringToNumFlags flags = 0); 

    /// Convert string to number of bytes. 
    ///
    /// String can contain "software" qualifiers: MB(megabyte), KB (kilobyte).
    /// Example: 100MB, 1024KB
    /// Note the qualifiers are power-of-2 based, aka kibi-, mebi- etc, so that
    /// 1KB = 1024B (not 1000B), 1MB = 1024KB = 1048576B, etc.
    ///
    /// @param str
    ///   String to be converted.
    /// @param flags
    ///   How to convert string to value.
    /// @param base
    ///   Numeric base of the number (before the qualifier). Allowed values
    ///   are 0, 2..20. Zero means to use the first characters to determine
    ///   the base - a leading "0x" or "0X" means base 16; otherwise a
    ///   leading 0 means base 8; otherwise base 10.
    ///   The base is limited to 20 to prevent 'K' from being interpreted as
    ///   a digit in the number.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted value.
    ///   - Otherwise, if fConvErr_NoThrow is not set, throw an exception.
    ///   - Otherwise, set errno to non-zero and return zero.
    /// @deprecated  Use StringToUInt8_DataSize(str, flags) instead.
    NCBI_DEPRECATED
    Uint8 StringToUInt8_DataSize(const CTempString str,
                                        TStringToNumFlags flags,
                                        int               base);

    /// Convert int to string.
    ///
    /// @param value
    ///   Integer value to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string IntToString(int value, TNumToStringFlags flags = 0,
                              int base = 10);

    string IntToString(unsigned int value, TNumToStringFlags flags = 0,
                              int base = 10);

    /// Convert int to string.
    ///
    /// @param out_str
    ///   Output string variable.
    /// @param value
    ///   Integer value to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void IntToString(string& out_str, int value, 
                            TNumToStringFlags flags = 0,
                            int               base  = 10);

    void IntToString(string& out_str, unsigned int value, 
                            TNumToStringFlags flags = 0,
                            int               base  = 10);

    /// Convert UInt to string.
    ///
    /// @param value
    ///   Integer value (unsigned long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string UIntToString(unsigned int      value,
                               TNumToStringFlags flags = 0,
                               int               base  = 10);

    string UIntToString(int               value,
                               TNumToStringFlags flags = 0,
                               int               base  = 10);

    /// Convert UInt to string.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Integer value (unsigned long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void UIntToString(string& out_str, unsigned int value,
                             TNumToStringFlags flags = 0,
                             int               base  = 10);

    void UIntToString(string& out_str, int value,
                             TNumToStringFlags flags = 0,
                             int               base  = 10);

    /// Convert Int to string.
    ///
    /// @param value
    ///   Integer value (long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string LongToString(long value, TNumToStringFlags flags = 0,
                               int base = 10);

    /// Convert Int to string.
    ///
    /// @param out_str
    ///   Output string variable.
    /// @param value
    ///   Integer value (long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void LongToString(string& out_str, long value, 
                             TNumToStringFlags flags = 0,
                             int               base  = 10);

    /// Convert unsigned long to string.
    ///
    /// @param value
    ///   Integer value (unsigned long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string ULongToString(unsigned long     value,
                                TNumToStringFlags flags = 0,
                                int               base  = 10);

    /// Convert unsigned long to string.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Integer value (unsigned long) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void ULongToString(string& out_str, unsigned long value,
                              TNumToStringFlags flags = 0,
                              int               base  = 10); 

    /// Convert Int8 to string.
    ///
    /// @param value
    ///   Integer value (Int8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string Int8ToString(Int8 value,
                               TNumToStringFlags flags = 0,
                               int               base  = 10); 

    /// Convert Int8 to string.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Integer value (Int8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void Int8ToString(string& out_str, Int8 value,
                             TNumToStringFlags flags = 0,
                             int               base  = 10);   

    /// Convert UInt8 to string.
    ///
    /// @param value
    ///   Integer value (UInt8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string UInt8ToString(Uint8 value,
                                TNumToStringFlags flags = 0,
                                int               base  = 10);

    /// Convert UInt8 to string.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Integer value (UInt8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void UInt8ToString(string& out_str, Uint8 value,
                              TNumToStringFlags flags = 0,
                              int               base  = 10);   

    /// Convert UInt8 to string using "software" qualifiers.
    /// 
    /// Result of conversion will be limited to max_digits digits so that e.g.
    /// 1024 will be converted to 1.02KB. Conversion will be made using
    /// rounding so that 1025 will be converted to 1.03KB. By default function
    /// uses 1000-based qualifiers (as in examples above) but with fDS_Binary
    /// flag it will use 1024-based qualifiers, e.g. 1100 will be converted to
    /// 1.07KiB. With fDS_ShortSuffix flag function will omit "B" in 1000-based
    /// and "iB" in 1024-based qualifiers. When the result of conversion doesn't
    /// need any qualifiers then the result of this function will be equivalent
    /// to result of UInt8ToString() above except if fDS_PutBSuffixToo flag
    /// is passed. In the latter case "B" will be added to the number.
    /// 
    /// Function will always try to use a maximum possible qualifier and
    /// a number with decimal point except if fDS_NoDecimalPoint flag is passed.
    /// In that case function will return only whole number and try to use a
    /// minimum possible qualifier (which makes difference only if
    /// max_digits > 3).
    ///
    /// @param value
    ///   Integer value (UInt8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param max_digits
    ///   Maximum number of digits to use (cannot be less than 3)
    /// @return
    ///   - If invalid flags are passed, throw an exception.
    ///   - If conversion succeeds, return the converted value.
    string UInt8ToString_DataSize(Uint8 value,
                                         TNumToStringFlags flags = 0,
                                         unsigned int max_digits = 3); 

    /// Convert UInt8 to string using "software" qualifiers.
    /// 
    /// See notes and details of how function works in the comments to 
    /// UInt8ToString_DataSize() above.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Integer value (UInt8) to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param max_digits
    ///   Maximum number of digits to use (cannot be less than 3)
    void UInt8ToString_DataSize(string& out_str,
                                       Uint8 value,
                                       TNumToStringFlags flags = 0,
                                       unsigned int max_digits = 3); 

    /// Convert double to string.
    ///
    /// @param value
    ///   Double value to be converted.
    /// @param precision
    ///   Precision value for conversion. If precision is more that maximum
    ///   for current platform, then it will be truncated to this maximum.
    ///   If it is negative, that double will be converted to number in
    ///   scientific notation.
    /// @param flags
    ///   How to convert value to string.
    ///   If double format flags are not specified, that next output format
    ///   will be used by default:
    ///     - fDoubleFixed,   if 'precision' >= 0.
    ///     - fDoubleGeneral, if 'precision' < 0.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string DoubleToString(double value, int precision = -1,
                                 TNumToStringFlags flags = 0);

    /// Convert double to string.
    ///
    /// @param out_str
    ///   Output string variable
    /// @param value
    ///   Double value to be converted.
    /// @param precision
    ///   Precision value for conversion. If precision is more that maximum
    ///   for current platform, then it will be truncated to this maximum.
    ///   If it is negative, that double will be converted to number in
    ///   scientific notation.
    /// @param flags
    ///   How to convert value to string.
    ///   If double format flags are not specified, that next output format
    ///   will be used by default:
    ///     - fDoubleFixed,   if 'precision' >= 0.
    ///     - fDoubleGeneral, if 'precision' < 0.
    /// @note
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value in 'out_str'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    void DoubleToString(string& out_str, double value,
                               int precision = -1,
                               TNumToStringFlags flags = 0);

    /// Convert double to string with specified precision and place the result
    /// in the specified buffer.
    ///
    /// @param value
    ///   Double value to be converted.
    /// @param precision
    ///   Precision value for conversion. If precision is more that maximum
    ///   for current platform, then it will be truncated to this maximum.
    /// @param buf
    ///   Put result of the conversion into this buffer.
    /// @param buf_size
    ///   Size of buffer, "buf".
    /// @param flags
    ///   How to convert value to string.
    ///   Default output format is fDoubleFixed.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     number of bytes stored in "buf", not counting the
    ///     terminating '\0'.
    ///   - Otherwise, set errno to non-zero, value of 'out_str' is undefined.
    SIZE_TYPE DoubleToString(double value, unsigned int precision,
                                    char* buf, SIZE_TYPE buf_size,
                                    TNumToStringFlags flags = 0);

    /// Convert double to string with specified precision and put the result
    /// into a character buffer, in scientific format.
    ///
    /// NOTE:
    ///   The output character buffer is NOT zero-terminated.
    ///   The decimal separator is dot, always.
    ///   This function DOES NOT check 'value' for being finite or not-a-number;
    ///   if it is, the result is unpredictable.
    ///   This function is less precise for a small fraction of values
    ///   (the difference is in the last significant digit) than its
    ///   'DoubleToString' siblings, but it is much faster.
    ///
    /// @param value
    ///   Double value to be converted.
    /// @param precision
    ///   Maximum number of significant digits to preserve. If precision is greater than
    ///   maximum for the current platform, it will be truncated to this maximum.
    /// @param buf
    ///   Put result of the conversion into this buffer.
    /// @param buf_size
    ///   Size of buffer, "buf".
    /// @return
    ///   The number of bytes written into "buf".
    SIZE_TYPE DoubleToStringPosix(double value,unsigned int precision,
                                         char* buf, SIZE_TYPE buf_size);


    /// Convert double to string with specified precision.
    /// 
    /// The result consists of three parts: significant digits, exponent and sign.
    /// For example, input value -12345.67 will produce
    /// buf = "1234567" , *dec = 4, and *sign = -1.
    /// NOTE:
    ///   The output character buffer is NOT zero-terminated.
    ///   The buffer is NOT padded with zeros.
    ///   This function DOES NOT check 'value' for being finite or not-a-number;
    ///   if it is, the result is unpredictable.
    ///   This function is less precise for a small fraction of values
    ///   (the difference is in the last significant digit) than its
    ///   'DoubleToString' siblings, but it is much faster.
    ///
    /// @param value
    ///   Double value to be converted.
    /// @param precision
    ///   Maximum number of significant digits to preserve. If precision is greater than
    ///   maximum for the current platform, it will be truncated to this maximum.
    /// @param buf
    ///   Put result of the conversion into this buffer.
    /// @param buf_size
    ///   Size of buffer, "buf".
    /// @param dec
    ///   Exponent
    /// @param sign
    ///   Sign of the value
    /// @return
    ///   The number of bytes written into "buf".
    SIZE_TYPE DoubleToString_Ecvt(double value,unsigned int precision,
                                         char* buf, SIZE_TYPE buf_size,
                                         int* dec, int* sign);

    /// Convert size_t to string.
    ///
    /// @param value
    ///   Value to be converted.
    /// @param flags
    ///   How to convert value to string.
    /// @param base
    ///   Radix base. Default is 10. Allowed values are 2..36.
    ///   Bases 8 and 16 do not add leading '0' and '0x' accordingly.
    ///   If necessary you should add it yourself.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return the
    ///     converted string value.
    ///   - Otherwise, set errno to non-zero and return empty string.
    string SizetToString(size_t value,
                                TNumToStringFlags flags = 0,
                                int               base  = 10);

    /// Convert bool to string.
    ///
    /// @param value
    ///   Boolean value to be converted.
    /// @return
    ///   One of: 'true, 'false'
    /// @note
    ///   Don't change errno.
    const string BoolToString(bool value);

    /// Convert string to bool.
    ///
    /// @param str
    ///   Boolean string value to be converted.  Can recognize
    ///   case-insensitive version as one of:  'true, 't', 'yes', 'y'
    ///   for TRUE; and  'false', 'f', 'no', 'n' for FALSE.
    /// @return
    ///   - If conversion succeeds, set errno to zero and return
    ///     TRUE or FALSE.
    ///   - Otherwise, set errno to non-zero and throw an exception.
    bool StringToBool(const CTempString str);

END_HBNSTR_SCOPE

inline
std::string NStr::IntToString(int value,
                         TNumToStringFlags flags, int base)
{
    string ret;
    IntToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::IntToString(unsigned int value,
                         TNumToStringFlags flags, int base)
{
    string ret;
    IntToString(ret, (int)value, flags, base);
    return ret;
}

inline
void NStr::IntToString(string& out_str, unsigned int value, 
                       TNumToStringFlags flags, int base)
{
    IntToString(out_str, (int)value, flags, base);
}

inline
std::string NStr::UIntToString(unsigned int value,
                          TNumToStringFlags flags, int base)
{
    string ret;
    ULongToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::UIntToString(int value,
                          TNumToStringFlags flags, int base)
{
    string ret;
    UIntToString(ret, (unsigned int)value, flags, base);
    return ret;
}

inline
void NStr::UIntToString(string& out_str, unsigned int value,
                        TNumToStringFlags flags, int base)
{
    ULongToString(out_str, value, flags, base);
}

inline
void NStr::UIntToString(string& out_str, int value,
                        TNumToStringFlags flags, int base)
{
    UIntToString(out_str, (unsigned int)value, flags, base);
}

inline
std::string NStr::LongToString(long value,
                          TNumToStringFlags flags, int base)
{
    string ret;
    LongToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::ULongToString(unsigned long value,
                           TNumToStringFlags flags, int base)
{
    string ret;
    ULongToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::Int8ToString(Int8 value,
                          TNumToStringFlags flags, int base)
{
    string ret;
    NStr::Int8ToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::UInt8ToString(Uint8 value,
                           TNumToStringFlags flags, int base)
{
    string ret;
    NStr::UInt8ToString(ret, value, flags, base);
    return ret;
}

inline
std::string NStr::UInt8ToString_DataSize(Uint8 value,
                                    TNumToStringFlags flags /* = 0 */,
                                    unsigned int max_digits /* = 3 */)
{
    string ret;
    NStr::UInt8ToString_DataSize(ret, value, flags, max_digits);
    return ret;
}

inline
std::string NStr::DoubleToString(double value, int precision,
                            TNumToStringFlags flags)
{
    string str;
    DoubleToString(str, value, precision, flags);
    return str;
}

#endif // __NUMERIC_STR_INTERCONV_HPP