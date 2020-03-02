#include "numeric_str_interconv.hpp"

#include <cerrno>
#include <clocale>
#include <cmath>
#include <cfloat>
#include <climits>
#include <vector>
#include <sstream>
#include <iostream>

#include "str_cmp.hpp"
#include "str_util.hpp"

BEGIN_HBNSTR_SCOPE

using namespace ncbi;
using namespace std;

EHbnStrConvErrorMode gStrConvErrorMode = EHbnStrConvErrorMode::eHbnStrConvErrorAbort;

EHbnStrConvErrorMode GetStrConvErrorMode() 
{ 
    return gStrConvErrorMode; 
}

void SetStrConvErrorMode(EHbnStrConvErrorMode mode) 
{ 
    gStrConvErrorMode = mode; 
} 

// Digits (up to base 36)
static const char kDigit[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";


static inline
SIZE_TYPE s_DiffPtr(const char* end, const char* start)
{
    return end ? (SIZE_TYPE)(end - start) : (SIZE_TYPE) 0;
}

static string MakeConvErrorMessage(const CTempString str, const char* to_type, const CTempString msg)
{
    string s;
    s.reserve(str.length() + msg.length() + 50);
    s += "Cannot convert string '";
    s += str;
    s += "' to ";
    s += to_type;
    if ( !msg.empty() ) {
        s += ", ";
        s += msg;
    }
    return s;
}

#define S2N_CONVERT_GUARD(flags) 
#define S2N_CONVERT_GUARD_EX(flags)

#define S2N_CONVERT_ERROR(to_type, msg, errcode, pos) do { \
    string s2n_err_ = MakeConvErrorMessage(str, #to_type, msg); \
    cerr << "[" << __FILE__ << ", " << __FUNCTION__ << ", " << __LINE__ << "] " << s2n_err_; \
    cerr << " (Error takes place at position " << pos << ")"; \
    if (GetStrConvErrorMode() == EHbnStrConvErrorMode::eHbnStrConvErrorSetErrno) { \
        cerr << endl; \
        errno = errcode; \
        return (to_type)(0); \
    } else { \
        hbn_assert(GetStrConvErrorMode() == EHbnStrConvErrorMode::eHbnStrConvErrorAbort); \
        cerr << "\n\tWe abort the program and exit with 'EXIT_FAILURE'" << endl; \
        abort(); \
    } \
} while(0)

#define S2N_CONVERT_ERROR_INVAL(to_type)                              \
    S2N_CONVERT_ERROR(to_type, kEmptyStr, EINVAL, pos)

#define S2N_CONVERT_ERROR_RADIX(to_type, msg)                         \
    S2N_CONVERT_ERROR(to_type, msg, EINVAL, pos)

#define S2N_CONVERT_ERROR_OVERFLOW(to_type)                           \
    S2N_CONVERT_ERROR(to_type, "overflow", ERANGE, pos)

#define CHECK_ENDPTR(to_type)                                         \
    if ( str[pos] ) {                                                 \
        S2N_CONVERT_ERROR(to_type, kEmptyStr, EINVAL, pos);           \
    }

#define CHECK_ENDPTR_SIZE(to_type)                                    \
    if ( pos < size ) {                                               \
        S2N_CONVERT_ERROR(to_type, kEmptyStr, EINVAL, pos);           \
    }

#define CHECK_COMMAS                                                  \
    /* Check on possible commas */                                    \
    if (flags & NStr::fAllowCommas) {                                 \
        if (ch == ',') {                                              \
            if ((numpos == pos)  ||                                   \
                ((comma >= 0)  &&  (comma != 3)) ) {                  \
                /* Not first comma, sitting on incorrect place */     \
                break;                                                \
            }                                                         \
            /* Skip it */                                             \
            comma = 0;                                                \
            pos++;                                                    \
            continue;                                                 \
        } else {                                                      \
            if (comma >= 0) {                                         \
                /* Count symbols between commas */                    \
                comma++;                                              \
            }                                                         \
        }                                                             \
    }

/// @internal
// Check that symbol 'ch' is good symbol for number with radix 'base'.
static inline
bool s_IsGoodCharForRadix(char ch, int base, int* value = 0)
{
    if ( base <= 10 ) {
        // shortcut for most frequent case
        int delta = ch-'0';
        if ( unsigned(delta) < unsigned(base) ) {
            if ( value ) {
                *value = delta;
            }
            return true;
        }
        return false;
    }
    if (!isalnum((unsigned char) ch)) {
        return false;
    }
    // Corresponding numeric value of *endptr
    int delta;
    if (isdigit((unsigned char) ch)) {
        delta = ch - '0';
    } else {
        ch = (char) tolower((unsigned char) ch);
        delta = ch - 'a' + 10;
    }
    if ( value ) {
        *value = delta;
    }
    return delta < base;
 }


// Skip all allowed chars (all except used for digit composition). 
// Update 'ptr' to current position in the string.
enum ESkipMode {
    eSkipAll,           // all symbols
    eSkipAllAllowed,    // all symbols, except digit/+/-/.
    eSkipSpacesOnly     // spaces only 
};

static inline
bool s_IsDecimalPoint(unsigned char ch, NStr::TStringToNumFlags  flags)
{
    if ( ch != '.' && ch != ',') {
        return false;
    }
    if (flags & NStr::fDecimalPosix) {
        return ch == '.';
    }
    else if (flags & NStr::fDecimalPosixOrLocal) {
        return ch == '.' || ch == ',';
    }
    struct lconv* conv = localeconv();
    return ch == *(conv->decimal_point);
}

static inline
void s_SkipAllowedSymbols(const CTempString       str,
                          SIZE_TYPE&              pos,
                          ESkipMode               skip_mode,
                          NStr::TStringToNumFlags flags)
{
    if (skip_mode == eSkipAll) {
        pos = str.length();
        return;
    }

    for ( SIZE_TYPE len = str.length();  pos < len;  ++pos ) {
        unsigned char ch = str[pos];
        if ( isdigit(ch)  ||  ch == '+' ||  ch == '-'  ||  s_IsDecimalPoint(ch,flags) ) {
            break;
        }
        if ( (skip_mode == eSkipSpacesOnly)  &&  !isspace(ch) ) {
            break;
        }
    }
}


// Check radix base. If it is zero, determine base using first chars
// of the string. Update 'base' value.
// Update 'ptr' to current position in the string.
static inline
bool s_CheckRadix(const CTempString str, SIZE_TYPE& pos, int& base)
{
    if ( base == 10  ||  base == 8 ) {
        // shortcut for most frequent case
        return true;
    }
    // Check base
    if ( base < 0  ||  base == 1  ||  base > 36 ) {
        return false;
    }
    // Try to determine base using first chars of the string
    unsigned char ch   = str[pos];
    unsigned char next = str[pos+1];
    if ( base == 0 ) {
        if ( ch != '0' ) {
            base = 10;
        } else if (next == 'x' || next == 'X') {
            base = 16;
        } else {
            base = 8;
        }
    }
    // Remove leading '0x' for hex numbers
    if ( base == 16 ) {
        if (ch == '0'  &&  (next == 'x' || next == 'X')) {
            pos += 2;
        }
    }
    return true;
}

int StringToInt(const CTempString str, TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD_EX(flags);
    Int8 value = StringToInt8(str, flags, base);
    if ( value < kMin_Int  ||  value > kMax_Int ) {
        S2N_CONVERT_ERROR(int, "overflow", ERANGE, 0);
    }
    return (int) value;
}

unsigned int
StringToUInt(const CTempString str, TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD_EX(flags);
    Uint8 value = StringToUInt8(str, flags, base);
    if ( value > kMax_UInt ) {
        S2N_CONVERT_ERROR(unsigned int, "overflow", ERANGE, 0);
    }
    return (unsigned int) value;
}

long StringToLong(const CTempString str, TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD_EX(flags);
    Int8 value = StringToInt8(str, flags, base);
    if ( value < kMin_Long  ||  value > kMax_Long ) {
        S2N_CONVERT_ERROR(long, "overflow", ERANGE, 0);
    }
    return (long) value;
}

unsigned long
StringToULong(const CTempString str, TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD_EX(flags);
    Uint8 value = StringToUInt8(str, flags, base);
    if ( value > kMax_ULong ) {
        S2N_CONVERT_ERROR(unsigned long, "overflow", ERANGE, 0);
    }
    return (unsigned long) value;
}

double StringToDoublePosix(const char* ptr, char** endptr, TStringToNumFlags flags)
{
    S2N_CONVERT_GUARD(NStr::fConvErr_NoThrow);

    const char* start = ptr;
    char c = *ptr++;

    // skip leading blanks
    while ( isspace((unsigned char)c) ) {
        c = *ptr++;
    }

    int sign = 0;
    if ( c == '-' ) {
        sign = -1;
        c = *ptr++;
    }
    else if ( c == '+' ) {
        sign = +1;
        c = *ptr++;
    }
    
    if (c == 0) {
        if (endptr) {
            *endptr = (char*)start;
        }
        //err_guard.Set(EINVAL);
        errno = EINVAL;
        return 0.;
    }

    // short-cut - single digit
    if ( !*ptr && c >= '0' && c <= '9' ) {
        if (endptr) {
            *endptr = (char*)ptr;
        }
        double result = c-'0';
        // some compilers fail to negate zero
        return sign < 0 ? (c == '0' ? -0. : -result) : result;
    }

    bool         dot = false, expn = false, anydigits = false;
    int          digits = 0, dot_position = 0;
    unsigned int first=0, second=0, first_mul=1;
    long double  second_mul = NCBI_CONST_LONGDOUBLE(1.),
                 third      = NCBI_CONST_LONGDOUBLE(0.);

    // up to exponent
    for ( ; ; c = *ptr++ ) {
        if (c >= '0' && c <= '9') {
            // digits: accumulate
            c = (char)(c - '0');
            anydigits = true;
            ++digits;
            if (first == 0) {
                first = c;
                if ( first == 0 ) {
                    // omit leading zeros
                    --digits;
                    if (dot) {
                        --dot_position;
                    }
                }
            } else if (digits <= 9) {
                // first 9 digits come to 'first'
                first = first*10 + c;
            } else if (digits <= 18) {
                // next 9 digits come to 'second'
                first_mul *= 10;
                second = second*10 + c;
            } else {
                // other digits come to 'third'
                second_mul *= NCBI_CONST_LONGDOUBLE(10.);
                third = third * NCBI_CONST_LONGDOUBLE(10.) + c;
            }
        }
        else if (c == '.') {
            // dot
            // if second dot, stop
            if (dot) {
                --ptr;
                break;
            }
            dot_position = digits;
            dot = true;
        }
        else if (c == 'e' || c == 'E') {
            // if exponent, stop
            if (!anydigits) {
                --ptr;
                break;
            }
            expn = true;
            break;
        }
        else {
            --ptr;
            if (!anydigits) {
                if ( !dot && (c == 'n' || c == 'N') &&
                     NStr::strncasecmp(ptr,"nan",3)==0) {
                    if (endptr) {
                        *endptr = (char*)(ptr+3);
                    }
                    return HUGE_VAL/HUGE_VAL; /* NCBI_FAKE_WARNING */
                }
                if ( (c == 'i' || c == 'I') ) {
                    if ( NStr::strncasecmp(ptr,"inf",3)==0) {
                        ptr += 3;
                        if ( NStr::strncasecmp(ptr,"inity",5)==0) {
                            ptr += 5;
                        }
                        if (endptr) {
                            *endptr = (char*)ptr;
                        }
                        return sign < 0 ? -HUGE_VAL : HUGE_VAL;
                    }
                }
            }
            break;
        }
    }
    // if no digits, stop now - error
    if (!anydigits) {
        if (endptr) {
            *endptr = (char*)start;
        }
        //err_guard.Set(EINVAL);
        errno = EINVAL;
        return 0.;
    }
    int exponent = dot ? dot_position - digits : 0;

    // read exponent
    if (expn && *ptr) {
        int expvalue = 0;
        bool expsign = false, expnegate= false;
        int expdigits= 0;
        for( ; ; ++ptr) {
            c = *ptr;
            // sign: should be no digits at this point
            if (c == '-' || c == '+') {
                // if there was sign or digits, stop
                if (expsign || expdigits) {
                    break;
                }
                expsign = true;
                expnegate = c == '-';
            }
            // digits: accumulate
            else if (c >= '0' && c <= '9') {
                ++expdigits;
                int newexpvalue = expvalue*10 + (c-'0');
                if (newexpvalue > expvalue) {
                    expvalue = newexpvalue;
                }
            }
            else {
                break;
            }
        }
        // if no digits, rollback
        if (!expdigits) {
            // rollback sign
            if (expsign) {
                --ptr;
            }
            // rollback exponent
            if (expn) {
                --ptr;
            }
        }
        else {
            exponent = expnegate ? exponent - expvalue : exponent + expvalue;
        }
    }
    long double ret;
    if ( first_mul > 1 ) {
        _ASSERT(first);
        ret = ((long double)first * first_mul + second)* second_mul + third;
    }
    else {
        _ASSERT(first_mul == 1);
        _ASSERT(second == 0);
        _ASSERT(second_mul == 1);
        _ASSERT(third == 0);
        ret = first;
    }
    // calculate exponent
    if ( first && exponent ) {
        // multiply by power of 10 only non-zero mantissa
        if (exponent > 2*DBL_MAX_10_EXP) {
            ret = (flags & fDecimalPosixFinite) ? DBL_MAX :  HUGE_VAL;
            //err_guard.Set(ERANGE);
            errno = ERANGE;
        } else if (exponent < 2*DBL_MIN_10_EXP) {
            ret = (flags & fDecimalPosixFinite) ? DBL_MIN : 0.;
            //err_guard.Set(ERANGE);
            errno = ERANGE;
        } else {
            if ( exponent > 0 ) {
                static const double mul1[16] = {
                    1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7,
                    1e8, 1e9, 1e10, 1e11, 1e12, 1e13, 1e14, 1e15
                };
                ret *= mul1[exponent&15];
                if ( exponent >>= 4 ) {
                    static const long double mul2[16] = {
                        NCBI_CONST_LONGDOUBLE(1e0),
                        NCBI_CONST_LONGDOUBLE(1e16),
                        NCBI_CONST_LONGDOUBLE(1e32),
                        NCBI_CONST_LONGDOUBLE(1e48),
                        NCBI_CONST_LONGDOUBLE(1e64),
                        NCBI_CONST_LONGDOUBLE(1e80),
                        NCBI_CONST_LONGDOUBLE(1e96),
                        NCBI_CONST_LONGDOUBLE(1e112),
                        NCBI_CONST_LONGDOUBLE(1e128),
                        NCBI_CONST_LONGDOUBLE(1e144),
                        NCBI_CONST_LONGDOUBLE(1e160),
                        NCBI_CONST_LONGDOUBLE(1e176),
                        NCBI_CONST_LONGDOUBLE(1e192),
                        NCBI_CONST_LONGDOUBLE(1e208),
                        NCBI_CONST_LONGDOUBLE(1e224),
                        NCBI_CONST_LONGDOUBLE(1e240)
                    };
                    ret *= mul2[exponent&15];
                    for ( exponent >>= 4; exponent; --exponent ) {
                        ret *= NCBI_CONST_LONGDOUBLE(1e256);
                    }
                }
                if (!isfinite(double(ret))) {
                    if (flags & fDecimalPosixFinite) {
                        ret = DBL_MAX;
                    }
                    //err_guard.Set(ERANGE);
                    errno = ERANGE;
                }
            }
            else {
                exponent = -exponent;
                static const long double mul1[16] = {
                    NCBI_CONST_LONGDOUBLE(1e-0),
                    NCBI_CONST_LONGDOUBLE(1e-1),
                    NCBI_CONST_LONGDOUBLE(1e-2),
                    NCBI_CONST_LONGDOUBLE(1e-3),
                    NCBI_CONST_LONGDOUBLE(1e-4),
                    NCBI_CONST_LONGDOUBLE(1e-5),
                    NCBI_CONST_LONGDOUBLE(1e-6),
                    NCBI_CONST_LONGDOUBLE(1e-7),
                    NCBI_CONST_LONGDOUBLE(1e-8),
                    NCBI_CONST_LONGDOUBLE(1e-9),
                    NCBI_CONST_LONGDOUBLE(1e-10),
                    NCBI_CONST_LONGDOUBLE(1e-11),
                    NCBI_CONST_LONGDOUBLE(1e-12),
                    NCBI_CONST_LONGDOUBLE(1e-13),
                    NCBI_CONST_LONGDOUBLE(1e-14),
                    NCBI_CONST_LONGDOUBLE(1e-15)
                };
                ret *= mul1[exponent&15];
                if ( exponent >>= 4 ) {
                    static const long double mul2[16] = {
                        NCBI_CONST_LONGDOUBLE(1e-0),
                        NCBI_CONST_LONGDOUBLE(1e-16),
                        NCBI_CONST_LONGDOUBLE(1e-32),
                        NCBI_CONST_LONGDOUBLE(1e-48),
                        NCBI_CONST_LONGDOUBLE(1e-64),
                        NCBI_CONST_LONGDOUBLE(1e-80),
                        NCBI_CONST_LONGDOUBLE(1e-96),
                        NCBI_CONST_LONGDOUBLE(1e-112),
                        NCBI_CONST_LONGDOUBLE(1e-128),
                        NCBI_CONST_LONGDOUBLE(1e-144),
                        NCBI_CONST_LONGDOUBLE(1e-160),
                        NCBI_CONST_LONGDOUBLE(1e-176),
                        NCBI_CONST_LONGDOUBLE(1e-192),
                        NCBI_CONST_LONGDOUBLE(1e-208),
                        NCBI_CONST_LONGDOUBLE(1e-224),
                        NCBI_CONST_LONGDOUBLE(1e-240)
                    };
                    ret *= mul2[exponent&15];
                    for ( exponent >>= 4; exponent; --exponent ) {
                        ret *= NCBI_CONST_LONGDOUBLE(1e-256);
                    }
                }
                if ( ret < DBL_MIN ) {
                    if (flags & fDecimalPosixFinite) {
                        ret = DBL_MIN;
                    }
                    //err_guard.Set(ERANGE);
                    errno = ERANGE;
                }
            }
        }
    }
    if ( sign < 0 ) {
        ret = -ret;
    }
    // done
    if (endptr) {
        *endptr = (char*)ptr;
    }
    return (double)ret;
}

/// @internal
static double s_StringToDouble(const char* str, size_t size,
                               NStr::TStringToNumFlags flags)
{
    _ASSERT(str[size] == '\0');
    if ((flags & NStr::fDecimalPosix) && (flags & NStr::fDecimalPosixOrLocal)) {
        //NCBI_THROW2(CStringException, eBadArgs,
        //            "NStr::StringToDouble():  mutually exclusive flags specified",0);
        errno = EINVAL;
        cerr << "[" << __FILE__ << ", " << __FUNCTION__ << ", " << __LINE__ << "] "
             << "NStr::StringToDouble():  mutually exclusive flags specified" << endl;
        return (double)(0.0);
    }
    S2N_CONVERT_GUARD_EX(flags);

    // Current position in the string
    SIZE_TYPE pos  = 0;

    // Skip allowed leading symbols
    if (flags & NStr::fAllowLeadingSymbols) {
        bool spaces = ((flags & NStr::fAllowLeadingSymbols) == 
                       NStr::fAllowLeadingSpaces);
        s_SkipAllowedSymbols(CTempString(str, size), pos,
                             spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
    }
    // Check mandatory sign
    if (flags & NStr::fMandatorySign) {
        switch (str[pos]) {
        case '-':
        case '+':
            break;
        default:
            S2N_CONVERT_ERROR_INVAL(double);
        }
    }
    // For consistency make additional check on incorrect leading symbols.
    // Because strtod() may just skip such symbols.
    if (!(flags & NStr::fAllowLeadingSymbols)) {
        char c = str[pos];
        if ( !isdigit((unsigned char)c)  &&  !s_IsDecimalPoint(c,flags)  &&  c != '-'  &&  c != '+') {
            S2N_CONVERT_ERROR_INVAL(double);
        }
    }

    // Conversion
    int& errno_ref = errno;
    errno_ref = 0;

    char* endptr = 0;
    const char* begptr = str + pos;

    double n;
    if (flags & NStr::fDecimalPosix) {
        n = NStr::StringToDoublePosix(begptr, &endptr, flags);
    } else {
        n = strtod(begptr, &endptr);
    }
    if (flags & NStr::fDecimalPosixOrLocal) {
        char* endptr2 = 0;
        double n2 = NStr::StringToDoublePosix(begptr, &endptr2, flags);
        if (!endptr || (endptr2 && endptr2 > endptr)) {
            n = n2;
            endptr = endptr2;
        }
    }
    if ( !endptr  ||  endptr == begptr ) {
        S2N_CONVERT_ERROR(double, kEmptyStr, EINVAL, s_DiffPtr(endptr, begptr) + pos);
    }
    // some libs set ERANGE, others do not
    // here, we do not consider ERANGE as error
    if ( errno_ref && errno_ref != ERANGE ) {
        S2N_CONVERT_ERROR(double, kEmptyStr, errno_ref, s_DiffPtr(endptr, begptr) + pos);
    }
    // special cases
    if ((flags & NStr::fDecimalPosixFinite) && n != 0. && !std::isnan(n))
    {
        bool is_negative = n < 0.;
        if (is_negative) {
            n = -n;
        }
        if ( n < DBL_MIN) {
            n = DBL_MIN;
        } else if (!isfinite(n)) {
            n = DBL_MAX;
        }
        if (is_negative) {
            n = -n;
        }
    }

    pos += s_DiffPtr(endptr, begptr);

    // Skip allowed trailing symbols
    if (flags & NStr::fAllowTrailingSymbols) {
        bool spaces = ((flags & NStr::fAllowTrailingSymbols) ==
                       NStr::fAllowTrailingSpaces);
        s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
    }
    CHECK_ENDPTR(double);
    return n;
}


double StringToDoubleEx(const char* str, size_t size,
                              TStringToNumFlags flags)
{
    return s_StringToDouble(str, size, flags);
}

double StringToDouble(const CTempStringEx str, TStringToNumFlags flags)
{
    size_t size = str.size();
    if ( str.HasZeroAtEnd() ) {
        // string has zero at the end already
        return s_StringToDouble(str.data(), size, flags);
    }
    char buf[256]; // small temporary buffer on stack for appending zero char
    if ( size < sizeof(buf) ) {
        memcpy(buf, str.data(), size);
        buf[size] = '\0';
        return s_StringToDouble(buf, size, flags);
    }
    else {
        // use std::string() to allocate memory for appending zero char
        return s_StringToDouble(string(str).c_str(), size, flags);
    }
}

Int8 StringToInt8(const CTempString str, TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD(flags);

    // Current position in the string
    SIZE_TYPE pos = 0;

    // Skip allowed leading symbols
    if (flags & fAllowLeadingSymbols) {
        bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
        s_SkipAllowedSymbols(str, pos,
                             spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
    }
    // Determine sign
    bool sign = false;
    switch (str[pos]) {
    case '-':
        sign = true;
        /*FALLTHRU*/
    case '+':
        pos++;
        break;
    default:
        if (flags & fMandatorySign) {
            S2N_CONVERT_ERROR_INVAL(Int8);
        }
        break;
    }
    SIZE_TYPE pos0 = pos;
    // Check radix base
    if ( !s_CheckRadix(str, pos, base) ) {
        S2N_CONVERT_ERROR_RADIX(Int8, "bad numeric base '" + 
                                NStr::IntToString(base)+ "'");
    }

    // Begin conversion
    Int8 n = 0;
    Int8 limdiv = base==10? kMax_I8 / 10: kMax_I8 / base;
    Int8 limoff = (base==10? kMax_I8 % 10: kMax_I8 % base) + (sign ? 1 : 0);

    // Number of symbols between two commas. '-1' means -- no comma yet.
    int       comma  = -1;  
    SIZE_TYPE numpos = pos;

    while (char ch = str[pos]) {
        int  delta;   // corresponding numeric value of 'ch'

        // Check on possible commas
        CHECK_COMMAS;
        // Sanity check
        if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
            break;
        }
        // Overflow check
        if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
            S2N_CONVERT_ERROR_OVERFLOW(Int8);
        }
        n *= base;
        n += delta;
        pos++;
    }

    // Last checks
    if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
        S2N_CONVERT_ERROR_INVAL(Int8);
    }
    // Skip allowed trailing symbols
    if (flags & fAllowTrailingSymbols) {
        bool spaces = ((flags & fAllowTrailingSymbols) ==
                       fAllowTrailingSpaces);
        s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
    }
    // Assign sign before the end pointer check
    n = sign ? -n : n;
    CHECK_ENDPTR(Int8);

    return n;
}

Uint8 StringToUInt8(const CTempString str,
                          TStringToNumFlags flags, int base)
{
    S2N_CONVERT_GUARD(flags);

    const TStringToNumFlags slow_flags =
        fMandatorySign|fAllowCommas|fAllowLeadingSymbols|fAllowTrailingSymbols;

    if ( base == 10  &&  (flags & slow_flags) == 0 ) {
        // fast conversion

        // Current position in the string
        CTempString::const_iterator ptr = str.begin(), end = str.end();

        // Determine sign
        if ( ptr != end && *ptr == '+' ) {
            ++ptr;
        }
        if ( ptr == end ) {
            S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
        }

        // Begin conversion
        Uint8 n = 0;

        const Uint8 limdiv = kMax_UI8/10;
        const int   limoff = int(kMax_UI8 % 10);

        do {
            char ch = *ptr;
            int  delta = ch - '0';
            if ( unsigned(delta) >= 10 ) {
                S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, ptr-str.begin());
            }
            // Overflow check
            if ( n >= limdiv && (n > limdiv || delta > limoff) ) {
                S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, ptr-str.begin());
            }
            n = n*10+delta;
        } while ( ++ptr != end );

        return n;
    }

    // Current position in the string
    SIZE_TYPE pos = 0, size = str.size();

    // Skip allowed leading symbols
    if (flags & fAllowLeadingSymbols) {
        bool spaces = ((flags & fAllowLeadingSymbols) == fAllowLeadingSpaces);
        s_SkipAllowedSymbols(str, pos,
                             spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
    }
    // Determine sign
    if (str[pos] == '+') {
        pos++;
    } else {
        if (flags & fMandatorySign) {
            S2N_CONVERT_ERROR_INVAL(Uint8);
        }
    }
    SIZE_TYPE pos0 = pos;

    // Begin conversion
    Uint8 n = 0;
    // Check radix base
    if ( !s_CheckRadix(str, pos, base) ) {
        S2N_CONVERT_ERROR_RADIX(Uint8, "bad numeric base '" +
                                NStr::IntToString(base) + "'");
    }

    Uint8 limdiv = kMax_UI8 / base;
    int   limoff = int(kMax_UI8 % base);

    // Number of symbols between two commas. '-1' means -- no comma yet.
    int       comma  = -1;  
    SIZE_TYPE numpos = pos;

    while (char ch = str[pos]) {
        int delta;  // corresponding numeric value of 'ch'

        // Check on possible commas
        CHECK_COMMAS;
        // Sanity check
        if ( !s_IsGoodCharForRadix(ch, base, &delta) ) {
            break;
        }
        // Overflow check
        if ( n >= limdiv  &&  (n > limdiv  ||  delta > limoff) ) {
            S2N_CONVERT_ERROR_OVERFLOW(Uint8);
        }
        n *= base;
        n += delta;
        pos++;
    }

    // Last checks
    if ( pos == pos0  || ((comma >= 0)  &&  (comma != 3)) ) {
        S2N_CONVERT_ERROR_INVAL(Uint8);
    }
    // Skip allowed trailing symbols
    if (flags & fAllowTrailingSymbols) {
        bool spaces = ((flags & fAllowTrailingSymbols) ==
                       fAllowTrailingSpaces);
        s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
    }
    CHECK_ENDPTR_SIZE(Uint8);
    return n;
}

/// @internal
static Uint8 s_DataSizeConvertQual(const CTempString       str,
                                   SIZE_TYPE&              pos, 
                                   Uint8                   value,
                                   NStr::TStringToNumFlags flags)
{
    S2N_CONVERT_GUARD(flags);

    unsigned char ch = str[pos];
    if ( !ch ) {
        return value;
    }

    ch = (unsigned char)toupper(ch);
    Uint8 v   = value;
    bool  err = false;

    switch(ch) {
    case 'K':
        pos++;
        if ((kMax_UI8 / 1024) < v) {
            err = true;
        }
        v *= 1024;
        break;
    case 'M':
        pos++;
        if ((kMax_UI8 / 1024 / 1024) < v) {
            err = true;
        }
        v *= 1024 * 1024;
        break;
    case 'G':
        pos++;
        if ((kMax_UI8 / 1024 / 1024 / 1024) < v) {
            err = true;
        }
        v *= 1024 * 1024 * 1024;
        break;
    default:
        // error -- the "qual" points to the last unprocessed symbol
        S2N_CONVERT_ERROR_INVAL(Uint8);
    }
    if ( err ) {
#ifndef DataSize
typedef Uint8 DataSize;
#endif
        S2N_CONVERT_ERROR_OVERFLOW(DataSize);
    }

    ch = str[pos];
    if ( ch  &&  toupper(ch) == 'B' ) {
        pos++;
    }
    return v;
}

Uint8 StringToUInt8_DataSize(const CTempString str, 
                                   TStringToNumFlags flags, 
                                   int               base)
{
    // We have a limited base range here
    if ( base < 2  ||  base > 16 ) {
        //NCBI_THROW2(CStringException, eConvert,  
        //            "Bad numeric base '" + NStr::IntToString(base)+ "'", 0);
        errno = EINVAL;
        cerr << "[" << __FILE__ << ", " << __FUNCTION__ << ", " << __LINE__ << "] "
             << "Bad numeric base '" << NStr::IntToString(base) << "'" << endl;
        return 0;    
    }
    S2N_CONVERT_GUARD_EX(flags);

    // Current position in the string
    SIZE_TYPE pos = 0;

    // Find end of number representation
    {{
        // Skip allowed leading symbols
        if (flags & fAllowLeadingSymbols) {
            bool spaces = ((flags & fAllowLeadingSymbols) ==
                           fAllowLeadingSpaces);
            s_SkipAllowedSymbols(str, pos,
                           spaces ? eSkipSpacesOnly : eSkipAllAllowed, flags);
        }
        // Determine sign
        if (str[pos] == '+') {
            pos++;
            // strip fMandatorySign flag
            flags &= ~fMandatorySign;
        } else {
            if (flags & fMandatorySign) {
                S2N_CONVERT_ERROR_INVAL(Uint8);
            }
        }
        // Check radix base
        if ( !s_CheckRadix(str, pos, base) ) {
            S2N_CONVERT_ERROR_RADIX(Uint8, "bad numeric base '" +
                                    NStr::IntToString(base) + "'");
        }
    }}

    SIZE_TYPE numpos = pos;
    char ch = str[pos];
    while (ch) {
        if ( !s_IsGoodCharForRadix(ch, base)  &&
             ((ch != ',')  ||  !(flags & fAllowCommas)) ) {
            break;
        }
        ch = str[++pos];
    }
    // If string is empty, just use whole remaining string for conversion
    // (for correct error reporting)
    if (pos-numpos == 0) {
        pos = str.length();
    }

    // Convert to number
    Uint8 n = StringToUInt8(CTempString(str.data()+numpos, pos-numpos),
                            flags, base);
    if ( !n && errno ) {
        // If exceptions are enabled that it has been already thrown.
        // The errno is also set, so just return a zero.
        return 0;
    }
    // Check trailer (KB, MB, ...)
    if ( ch ) {
        n = s_DataSizeConvertQual(str, pos, n, flags);
    }
    // Skip allowed trailing symbols
    if (flags & fAllowTrailingSymbols) {
        bool spaces = ((flags & fAllowTrailingSymbols) ==
                       fAllowTrailingSpaces);
        s_SkipAllowedSymbols(str, pos, spaces ? eSkipSpacesOnly : eSkipAll, flags);
    }
    CHECK_ENDPTR(Uint8);
    return n;
}

Uint8 StringToUInt8_DataSize(const CTempString str,
                                   TStringToNumFlags flags /* = 0 */)
{
    TStringToNumFlags allowed_flags = /* fConvErr_NoThrow + */
                                      fMandatorySign +
                                      fAllowCommas +
                                      fAllowLeadingSymbols +
                                      fAllowTrailingSymbols +
                                      fDS_ForceBinary +
                                      fDS_ProhibitFractions +
                                      fDS_ProhibitSpaceBeforeSuffix;

    if ((flags & allowed_flags) != flags) {
        //NCBI_THROW2(CStringException, eConvert, "Wrong set of flags", 0);
        errno = EINVAL;
        cerr << "[" << __FILE__ << ", " << __FUNCTION__ << ", " << __LINE__ << "] "
             << "Wrong set of flags" << endl;
        return 0;
    }
    S2N_CONVERT_GUARD(flags);

    const char* str_ptr = str.data();
    const char* str_end = str_ptr + str.size();
    if (flags & fAllowLeadingSymbols) {
        bool allow_all = (flags & fAllowLeadingSymbols) != fAllowLeadingSpaces;
        for (; str_ptr < str_end; ++str_ptr) {
            char c = *str_ptr;
            if (isdigit(c))
                break;
            if (isspace(c))
                continue;
            if ((c == '+'  ||  c == '-')  &&  (flags & fMandatorySign)
                &&  str_ptr + 1 < str_end  &&  isdigit(*(str_ptr + 1)))
            {
                break;
            }
            if (!allow_all)
                break;
        }
    }

    if (str_ptr < str_end  &&  *str_ptr == '+') {
        ++str_ptr;
    }
    else if ((str_ptr < str_end  &&  *str_ptr == '-')
             ||  (flags & fMandatorySign))
    {
        S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    }

    const char* num_start = str_ptr;
    bool have_dot = false;
    bool allow_commas = (flags & fAllowCommas) != 0;
    bool allow_dot = (flags & fDS_ProhibitFractions) == 0;
    Uint4 digs_pre_dot = 0, digs_post_dot = 0;

    for (; str_ptr < str_end; ++str_ptr) {
        char c = *str_ptr;
        if (isdigit(c)) {
            if (have_dot)
                ++digs_post_dot;
            else
                ++digs_pre_dot;
        }
        else if (c == '.'  &&  allow_dot) {
            if (have_dot  ||  str_ptr == num_start)
                break;
            if (*(str_ptr - 1) == ',') {
                --str_ptr;
                break;
            }
            have_dot = true;
        }
        else if (c == ','  &&  allow_commas) {
            if (have_dot  ||  str_ptr == num_start)
                break;
            if (*(str_ptr - 1) == ',') {
                --str_ptr;
                break;
            }
        }
        else
            break;
    }
    if (have_dot  &&  digs_post_dot == 0)
        --str_ptr;
    else if (str_ptr > num_start  &&  *(str_ptr - 1) == ',')
        --str_ptr;

    const char* num_end = str_ptr;
    if (num_start == num_end) {
        S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    }
    if (str_ptr < str_end  &&  *str_ptr == ' '
        &&  !(flags & fDS_ProhibitSpaceBeforeSuffix))
    {
        ++str_ptr;
    }
    char suff_c = 0;
    if (str_ptr < str_end)
        suff_c = (char)toupper(*str_ptr);

    static const char s_Suffixes[] = {'K', 'M', 'G', 'T', 'P', 'E'};
    static const char* const s_BinCoefs[] = {"1024", "1048576", "1073741824",
                                             "1099511627776",
                                             "1125899906842624",
                                             "1152921504606846976"};
    static const Uint4 s_NumSuffixes = (Uint4)(sizeof(s_Suffixes) / sizeof(s_Suffixes[0]));

    bool binary_suff = (flags & fDS_ForceBinary) != 0;
    Uint4 suff_idx = 0;
    for (; suff_idx < s_NumSuffixes; ++suff_idx) {
        if (suff_c == s_Suffixes[suff_idx])
            break;
    }
    if (suff_idx < s_NumSuffixes) {
        ++str_ptr;
        if (str_ptr + 1 < str_end  &&  toupper(*str_ptr) == 'I'
            &&  toupper(*(str_ptr + 1)) == 'B')
        {
            str_ptr += 2;
            binary_suff = true;
        }
        else if (str_ptr < str_end  &&  toupper(*str_ptr) == 'B')
            ++str_ptr;
    }
    else if (suff_c == 'B') {
        ++str_ptr;
    }
    else if (*(str_ptr - 1) == ' ')
        --str_ptr;

    if (flags & fAllowTrailingSymbols) {
        bool allow_all = (flags & fAllowTrailingSymbols) != fAllowTrailingSpaces;
        for (; str_ptr < str_end; ++str_ptr) {
            char c = *str_ptr;
            if (isspace(c))
                continue;
            if (!allow_all)
                break;
        }
    }
    if (str_ptr != str_end) {
        S2N_CONVERT_ERROR(Uint8, kEmptyStr, EINVAL, str_ptr - str.data());
    }

    Uint4 orig_digs = digs_pre_dot + digs_post_dot;
    //AutoArray<Uint1> orig_num(orig_digs);
    vector<Uint1> orig_num(orig_digs);
    str_ptr = num_start;
    for (Uint4 i = 0; str_ptr < num_end; ++str_ptr) {
        if (*str_ptr == ','  ||  *str_ptr == '.')
            continue;
        orig_num[i++] = Uint1(*str_ptr - '0');
    }

    Uint1* num_to_conv = orig_num.data();
    //Uint1* num_to_conv = orig_num.get();
    Uint4 digs_to_conv = digs_pre_dot;
    //AutoArray<Uint1> mul_num;
    vector<Uint1> mul_num;
    if (binary_suff  &&  suff_idx < s_NumSuffixes) {
        const char* coef = s_BinCoefs[suff_idx];
        Uint4 coef_size = Uint4(strlen(coef));
        mul_num.resize(orig_digs + coef_size);
        memset(mul_num.data(), 0, mul_num.size());
        //mul_num = new Uint1[orig_digs + coef_size];
        //memset(mul_num.get(), 0, orig_digs + coef_size);
        for (Uint4 coef_i = 0; coef_i < coef_size; ++coef_i) {
            Uint1 coef_d = Uint1(coef[coef_i] - '0');
            Uint1 carry = 0;
            Uint4 res_idx = orig_digs + coef_i;
            for (int orig_i = orig_digs - 1; orig_i >= 0; --orig_i, --res_idx) {
                Uint1 orig_d = orig_num[orig_i];
                Uint1 res_d = Uint1(coef_d * orig_d + carry + mul_num[res_idx]);
                carry = 0;
                while (res_d >= 10) {
                    res_d = (Uint1)(res_d - 10); // res_d -= 10;
                    ++carry;
                }
                mul_num[res_idx] = res_d;
            }
            _ASSERT(carry <= 9);
            for (; carry != 0; --res_idx) {
                Uint1 res_d = Uint1(mul_num[res_idx] + carry);
                carry = 0;
                while (res_d >= 10) {
                    res_d = (Uint1)(res_d - 10); // res_d -= 10;
                    ++carry;
                }
                mul_num[res_idx] = res_d;
            }
        }
        digs_to_conv = orig_digs + coef_size - digs_post_dot;
        //num_to_conv = mul_num.get();
        num_to_conv = mul_num.data();
        while (digs_to_conv > 1  &&  *num_to_conv == 0) {
            --digs_to_conv;
            ++num_to_conv;
        }
    }
    else if (suff_idx < s_NumSuffixes) {
        Uint4 coef_size = (suff_idx + 1) * 3;
        if (coef_size <= digs_post_dot) {
            digs_to_conv += coef_size;
            digs_post_dot -= coef_size;
        }
        else {
            digs_to_conv += digs_post_dot;
            coef_size -= digs_post_dot;
            digs_post_dot = 0;
            mul_num.resize(digs_to_conv + coef_size);
            memmove(mul_num.data(), num_to_conv, digs_to_conv);
            memset(mul_num.data() + digs_to_conv, 0, coef_size);
            num_to_conv = mul_num.data();
            //mul_num = new Uint1[digs_to_conv + coef_size];
            //memmove(mul_num.get(), num_to_conv, digs_to_conv);
            //memset(mul_num.get() + digs_to_conv, 0, coef_size);
            //num_to_conv = mul_num.get();
            digs_to_conv += coef_size;
        }
    }

    const Uint8 limdiv = kMax_UI8/10;
    const int   limoff = int(kMax_UI8 % 10);
    Uint8 n = 0;
    for (Uint4 i = 0; i < digs_to_conv; ++i) {
        Uint1 d = num_to_conv[i];
        if (n >= limdiv  &&  (n > limdiv  ||  d > limoff)) {
            S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, i);
        }
        n *= 10;
        n += d;
    }
    if (digs_post_dot != 0  &&  num_to_conv[digs_to_conv] >= 5) {
        if (n == kMax_UI8) {
            S2N_CONVERT_ERROR(Uint8, kEmptyStr, ERANGE, digs_to_conv);
        }
        ++n;
    }
    return n;
}

/// @internal
static void s_SignedToString(string&                 out_str,
                             unsigned long           value,
                             long                    svalue,
                             NStr::TNumToStringFlags flags,
                             int                     base)
{
    const SIZE_TYPE kBufSize = CHAR_BIT * sizeof(value);
    char  buffer[kBufSize];
    char* pos = buffer + kBufSize;
    
    if ( base == 10 ) {
        if ( svalue < 0 ) {
            value = static_cast<unsigned long>(-svalue);
        }
        
        if ( (flags & NStr::fWithCommas) ) {
            int cnt = -1;
            do {
                if (++cnt == 3) {
                    *--pos = ',';
                    cnt = 0;
                }
                unsigned long a = '0'+value;
                value /= 10;
                *--pos = char(a - value*10);
            } while ( value );
        }
        else {
            do {
                unsigned long a = '0'+value;
                value /= 10;
                *--pos = char(a - value*10);
            } while ( value );
        }

        if (svalue < 0)
            *--pos = '-';
        else if (flags & NStr::fWithSign)
            *--pos = '+';
    }
    else if ( base == 16 ) {
        do {
            *--pos = kDigit[value % 16];
            value /= 16;
        } while ( value );
    }
    else {
        do {
            *--pos = kDigit[value % base];
            value /= base;
        } while ( value );
    }

    out_str.assign(pos, buffer + kBufSize - pos);
}


void IntToString(string& out_str, int svalue,
                       TNumToStringFlags flags, int base)
{
    if ( base < 2  ||  base > 36 ) {
        errno = EINVAL;
        return;
    }
    unsigned int value = static_cast<unsigned int>(svalue);
    
    if ( base == 10  &&  svalue < 0 ) {
        value = static_cast<unsigned int>(-svalue);
    }
    s_SignedToString(out_str, value, svalue, flags, base);
    errno = 0;
}

void LongToString(string& out_str, long svalue,
                       TNumToStringFlags flags, int base)
{
    if ( base < 2  ||  base > 36 ) {
        errno = EINVAL;
        return;
    }
    unsigned long value = static_cast<unsigned long>(svalue);
    
    if ( base == 10  &&  svalue < 0 ) {
        value = static_cast<unsigned long>(-svalue);
    }
    s_SignedToString(out_str, value, svalue, flags, base);
    errno = 0;
}

void ULongToString(string&          out_str,
                        unsigned long     value,
                        TNumToStringFlags flags,
                        int               base)
{
    if ( base < 2  ||  base > 36 ) {
        errno = EINVAL;
        return;
    }
    const SIZE_TYPE kBufSize = CHAR_BIT * sizeof(value);
    char  buffer[kBufSize];
    char* pos = buffer + kBufSize;

    if ( base == 10 ) {
        if ( (flags & fWithCommas) ) {
            int cnt = -1;
            do {
                if (++cnt == 3) {
                    *--pos = ',';
                    cnt = 0;
                }
                unsigned long a = '0'+value;
                value /= 10;
                *--pos = char(a - value*10);
            } while ( value );
        }
        else {
            do {
                unsigned long a = '0'+value;
                value /= 10;
                *--pos = char(a - value*10);
            } while ( value );
        }

        if ( (flags & fWithSign) ) {
            *--pos = '+';
        }
    }
    else if ( base == 16 ) {
        do {
            *--pos = kDigit[value % 16];
            value /= 16;
        } while ( value );
    }
    else {
        do {
            *--pos = kDigit[value % base];
            value /= base;
        } while ( value );
    }

    out_str.assign(pos, buffer + kBufSize - pos);
    errno = 0;
}

// On some platforms division of Int8 is very slow,
// so will try to optimize it working with chunks.
// Works only for radix base == 10.

#define PRINT_INT8_CHUNK 1000000000
#define PRINT_INT8_CHUNK_SIZE 9

/// @internal
static char* s_PrintUint8(char*                   pos,
                          Uint8                   value,
                          NStr::TNumToStringFlags flags,
                          int                     base)
{
    if ( base == 10 ) {
        if ( (flags & NStr::fWithCommas) ) {
            int cnt = -1;
#ifdef PRINT_INT8_CHUNK
            // while n doesn't fit in Uint4 process the number
            // by 9-digit chunks within 32-bit Uint4
            while ( value & ~Uint8(Uint4(~0)) ) {
                Uint4 chunk = Uint4(value);
                value /= PRINT_INT8_CHUNK;
                chunk -= PRINT_INT8_CHUNK*Uint4(value);
                char* end = pos - PRINT_INT8_CHUNK_SIZE - 2; // 9-digit chunk should have 2 commas
                do {
                    if (++cnt == 3) {
                        *--pos = ',';
                        cnt = 0;
                    }
                    Uint4 a = '0'+chunk;
                    chunk /= 10;
                    *--pos = char(a-10*chunk);
                } while ( pos != end );
            }
            // process all remaining digits in 32-bit number
            Uint4 chunk = Uint4(value);
            do {
                if (++cnt == 3) {
                    *--pos = ',';
                    cnt = 0;
                }
                Uint4 a = '0'+chunk;
                chunk /= 10;
                *--pos = char(a-10*chunk);
            } while ( chunk );
#else
            do {
                if (++cnt == 3) {
                    *--pos = ',';
                    cnt = 0;
                }
                Uint8 a = '0'+value;
                value /= 10;
                *--pos = char(a - 10*value);
            } while ( value );
#endif
        }
        else {
#ifdef PRINT_INT8_CHUNK
            // while n doesn't fit in Uint4 process the number
            // by 9-digit chunks within 32-bit Uint4
            while ( value & ~Uint8(Uint4(~0)) ) {
                Uint4 chunk = Uint4(value);
                value /= PRINT_INT8_CHUNK;
                chunk -= PRINT_INT8_CHUNK*Uint4(value);
                char* end = pos - PRINT_INT8_CHUNK_SIZE;
                do {
                    Uint4 a = '0'+chunk;
                    chunk /= 10;
                    *--pos = char(a-10*chunk);
                } while ( pos != end );
            }
            // process all remaining digits in 32-bit number
            Uint4 chunk = Uint4(value);
            do {
                Uint4 a = '0'+chunk;
                chunk /= 10;
                *--pos = char(a-10*chunk);
            } while ( chunk );
#else
            do {
                Uint8 a = '0'+value;
                value /= 10;
                *--pos = char(a-10*value);
            } while ( value );
#endif
        }
    }
    else if ( base == 16 ) {
        do {
            *--pos = kDigit[value % 16];
            value /= 16;
        } while ( value );
    }
    else {
        do {
            *--pos = kDigit[value % base];
            value /= base;
        } while ( value );
    }
    return pos;
}


void Int8ToString(string& out_str, Int8 svalue,
                        TNumToStringFlags flags, int base)
{
    if ( base < 2  ||  base > 36 ) {
        //CNcbiError::SetErrno(errno = EINVAL);
        errno = EINVAL;
        return;
    }
    Uint8 value;
    if (base == 10) {
        value = static_cast<Uint8>(svalue<0?-svalue:svalue);
    } else {
        value = static_cast<Uint8>(svalue);
    }
    const SIZE_TYPE kBufSize = CHAR_BIT * sizeof(value);
    char  buffer[kBufSize];

    char* pos = s_PrintUint8(buffer + kBufSize, value, flags, base);

    if (base == 10) {
        if (svalue < 0)
            *--pos = '-';
        else if (flags & fWithSign)
            *--pos = '+';
    }
    out_str.assign(pos, buffer + kBufSize - pos);
    errno = 0;
}

void UInt8ToString(string& out_str, Uint8 value,
                         TNumToStringFlags flags, int base)
{
    if ( base < 2  ||  base > 36 ) {
        //CNcbiError::SetErrno(errno = EINVAL);
        errno = EINVAL;
        return;
    }
    const SIZE_TYPE kBufSize = CHAR_BIT  * sizeof(value);
    char  buffer[kBufSize];

    char* pos = s_PrintUint8(buffer + kBufSize, value, flags, base);

    if ( (base == 10)  &&  (flags & fWithSign) ) {
        *--pos = '+';
    }
    out_str.assign(pos, buffer + kBufSize - pos);
    errno = 0;
}

#ifdef HAVE_INTTYPES_H
#  define NCBI_CONST_INT8(v)     INT64_C(v)
#  define NCBI_CONST_UINT8(v)   UINT64_C(v)
#  define NCBI_INT8_FORMAT_SPEC  PRId64
#  define NCBI_UINT8_FORMAT_SPEC PRIu64
#elif (SIZEOF_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##L
#  define NCBI_CONST_UINT8(v)  v##UL
#  define NCBI_INT8_FORMAT_SPEC   "ld"
#  define NCBI_UINT8_FORMAT_SPEC  "lu"
#elif (SIZEOF_LONG_LONG == 8)
#  define NCBI_CONST_INT8(v)   v##LL
#  define NCBI_CONST_UINT8(v)  v##ULL
#  if defined(__MINGW32__)  ||  defined(__MINGW64__)
#    define NCBI_INT8_FORMAT_SPEC   "I64d"
#    define NCBI_UINT8_FORMAT_SPEC  "I64u"
#  else
#    define NCBI_INT8_FORMAT_SPEC   "lld"
#    define NCBI_UINT8_FORMAT_SPEC  "llu"
#  endif
#elif defined(NCBI_USE_INT64)
#  define NCBI_CONST_INT8(v)   v##i64
#  define NCBI_CONST_UINT8(v)  v##ui64
#  define NCBI_INT8_FORMAT_SPEC   "I64d"
#  define NCBI_UINT8_FORMAT_SPEC  "I64u"
#else
#  define NCBI_CONST_INT8(v)   v
#  define NCBI_CONST_UINT8(v)  v
#  define NCBI_INT8_FORMAT_SPEC   "d"
#  define NCBI_UINT8_FORMAT_SPEC  "u"
#endif
#if (SIZEOF_LONG_DOUBLE > SIZEOF_DOUBLE)
#  define NCBI_CONST_LONGDOUBLE(v)   v##L
#else
#  define NCBI_CONST_LONGDOUBLE(v)   v
#endif

void UInt8ToString_DataSize(string& out_str,
                                  Uint8 value,
                                  TNumToStringFlags flags /* = 0 */,
                                  unsigned int max_digits /* = 3 */)
{
    TNumToStringFlags allowed_flags = fWithSign +
                                      fWithCommas + 
                                      fDS_Binary + 
                                      fDS_NoDecimalPoint +
                                      fDS_PutSpaceBeforeSuffix +
                                      fDS_ShortSuffix +
                                      fDS_PutBSuffixToo;

    if ((flags & allowed_flags) != flags) {
        //NCBI_THROW2(CStringException, eConvert, "Wrong set of flags", 0);
        errno = EINVAL;
        cerr << "Wrong set of flags" << endl;
        return;
    }

    if (max_digits < 3)
        max_digits = 3;

    static const char s_Suffixes[] = {'K', 'M', 'G', 'T', 'P', 'E'};
    static const Uint4 s_NumSuffixes = Uint4(sizeof(s_Suffixes) / sizeof(s_Suffixes[0]));

    static const SIZE_TYPE kBufSize = 50;
    char  buffer[kBufSize];
    char* num_start;
    char* dot_ptr;
    char* num_end;
    Uint4 digs_pre_dot, suff_idx;

    if (!(flags &fDS_Binary)) {
        static const Uint8 s_Coefs[] = {1000, 1000000, 1000000000,
                                        NCBI_CONST_UINT8(1000000000000),
                                        NCBI_CONST_UINT8(1000000000000000),
                                        NCBI_CONST_UINT8(1000000000000000000)};
        suff_idx = 0;
        for (; suff_idx < s_NumSuffixes; ++suff_idx) {
            if (value < s_Coefs[suff_idx])
                break;
        }
        num_start = s_PrintUint8(buffer + kBufSize, value, 0, 10);
        num_start[-1] = '0';
        dot_ptr = buffer + kBufSize - 3 * suff_idx;
        digs_pre_dot = Uint4(dot_ptr - num_start);
        if (!(flags & fDS_NoDecimalPoint)) {
            num_end = min(buffer + kBufSize, dot_ptr + (max_digits - digs_pre_dot));
        }
        else {
            while (suff_idx > 0  &&  max_digits - digs_pre_dot >= 3) {
                --suff_idx;
                digs_pre_dot += 3;
                dot_ptr += 3;
            }
            num_end = dot_ptr;
        }
        char* round_dig = num_end - 1;
        if (num_end < buffer + kBufSize  &&  *num_end >= '5')
            ++(*round_dig);
        while (*round_dig == '0' + 10) {
            *round_dig = '0';
            --round_dig;
            ++(*round_dig);
        }
        if (round_dig < num_start) {
            _ASSERT(num_start - round_dig == 1);
            num_start = round_dig;
            ++digs_pre_dot;
            if (!(flags & fDS_NoDecimalPoint)) {
                if (digs_pre_dot > 3) {
                    ++suff_idx;
                    digs_pre_dot -= 3;
                    dot_ptr -= 3;
                }
                --num_end;
            }
            else {
                if (digs_pre_dot > max_digits) {
                    ++suff_idx;
                    digs_pre_dot -= 3;
                    dot_ptr -= 3;
                    num_end = dot_ptr;
                }
            }
        }
    }
    else {
        static const Uint8 s_Coefs[] = {1, 1024, 1048576, 1073741824,
                                        NCBI_CONST_UINT8(1099511627776),
                                        NCBI_CONST_UINT8(1125899906842624),
                                        NCBI_CONST_UINT8(1152921504606846976)};

        suff_idx = 1;
        for (; suff_idx < s_NumSuffixes; ++suff_idx) {
            if (value < s_Coefs[suff_idx])
                break;
        }
        bool can_try_another = true;
try_another_suffix:
        Uint8 mul_coef = s_Coefs[suff_idx - 1];
        Uint8 whole_num = value / mul_coef;
        if (max_digits == 3  &&  whole_num >= 1000) {
            ++suff_idx;
            goto try_another_suffix;
        }
        num_start = s_PrintUint8(buffer + kBufSize, whole_num, 0, 10);
        num_start[-1] = '0';
        digs_pre_dot = Uint4(buffer + kBufSize - num_start);
        if (max_digits - digs_pre_dot >= 3  &&  (flags & fDS_NoDecimalPoint)
            &&  suff_idx != 1  &&  can_try_another)
        {
            Uint4 new_suff = suff_idx - 1;
try_even_more_suffix:
            Uint8 new_num = value / s_Coefs[new_suff - 1];
            char* new_start = s_PrintUint8(buffer + kBufSize / 2, new_num, 0, 10);
            Uint4 new_digs = Uint4(buffer + kBufSize / 2 - new_start);
            if (new_digs <= max_digits) {
                if (max_digits - digs_pre_dot >= 3  &&  new_suff != 1) {
                    --new_suff;
                    goto try_even_more_suffix;
                }
                suff_idx = new_suff;
                can_try_another = false;
                goto try_another_suffix;
            }
            if (new_suff != suff_idx - 1) {
                suff_idx = new_suff + 1;
                can_try_another = false;
                goto try_another_suffix;
            }
        }
        memcpy(buffer, num_start - 1, digs_pre_dot + 1);
        num_start = buffer + 1;
        dot_ptr = num_start + digs_pre_dot;
        Uint4 cnt_more_digs = 1;
        if (!(flags & fDS_NoDecimalPoint))
            cnt_more_digs += min(max_digits - digs_pre_dot, 3 * (suff_idx - 1));
        num_end = dot_ptr;
        Uint8 left_val = value - whole_num * mul_coef;
        do {
            left_val *= 10;
            Uint1 d = Uint1(left_val / mul_coef);
            *num_end = char(d + '0');
            ++num_end;
            left_val -= d * mul_coef;
            --cnt_more_digs;
        }
        while (cnt_more_digs != 0);
        --num_end;

        char* round_dig = num_end - 1;
        if (*num_end >= '5')
            ++(*round_dig);
        while (*round_dig == '0' + 10) {
            *round_dig = '0';
            --round_dig;
            ++(*round_dig);
        }
        if (round_dig < num_start) {
            _ASSERT(round_dig == buffer);
            num_start = round_dig;
            ++digs_pre_dot;
            if (digs_pre_dot > max_digits) {
                ++suff_idx;
                goto try_another_suffix;
            }
            if (num_end != dot_ptr)
                --num_end;
        }
        if (!(flags & fDS_NoDecimalPoint)  &&  digs_pre_dot == 4
            &&  num_start[0] == '1'  &&  num_start[1] == '0'
            &&  num_start[2] == '2'  &&  num_start[3] == '4')
        {
            ++suff_idx;
            goto try_another_suffix;
        }

        --suff_idx;
    }

    out_str.erase();
    if (flags & fWithSign)
        out_str.append(1, '+');
    if (!(flags & fWithCommas)  ||  digs_pre_dot <= 3) {
        out_str.append(num_start, digs_pre_dot);
    }
    else {
        Uint4 digs_first = digs_pre_dot % 3;
        out_str.append(num_start, digs_first);
        char* left_ptr = num_start + digs_first;
        Uint4 digs_left = digs_pre_dot - digs_first;
        while (digs_left != 0) {
            out_str.append(1, ',');
            out_str.append(left_ptr, 3);
            left_ptr += 3;
            digs_left -= 3;
        }
    }
    if (num_end != dot_ptr) {
        out_str.append(1, '.');
        out_str.append(dot_ptr, num_end - dot_ptr);
    }

    if (suff_idx == 0) {
        if (flags & fDS_PutBSuffixToo) {
            if (flags & fDS_PutSpaceBeforeSuffix)
                out_str.append(1, ' ');
            out_str.append(1, 'B');
        }
    }
    else {
        --suff_idx;
        if (flags & fDS_PutSpaceBeforeSuffix)
            out_str.append(1, ' ');
        out_str.append(1, s_Suffixes[suff_idx]);
        if (!(flags & fDS_ShortSuffix)) {
            if (flags & fDS_Binary)
                out_str.append(1, 'i');
            out_str.append(1, 'B');
        }
    }
    errno = 0;
}

// A maximal double precision used in the double to string conversion
#if defined(NCBI_OS_MSWIN)
    const int kMaxDoublePrecision = 200;
#else
    const int kMaxDoublePrecision = 308;
#endif
// A maximal size of a double value in a string form.
// Exponent size + sign + dot + ending '\0' + max.precision
const int kMaxDoubleStringSize = 308 + 3 + kMaxDoublePrecision;


void DoubleToString(string& out_str, double value,
                          int precision, TNumToStringFlags flags)
{
    char buffer[kMaxDoubleStringSize];
    if (precision >= 0 ||
        ((flags & fDoublePosix) && (!isfinite(value) || value == 0.))) {
        SIZE_TYPE n = DoubleToString(value, precision, buffer,
                                     kMaxDoubleStringSize, flags);
        buffer[n] = '\0';
    } else {
        const char* format;
        switch (flags & fDoubleGeneral) {
            case fDoubleFixed:
                format = "%f";
                break;
            case fDoubleScientific:
                format = "%e";
                break;
            case fDoubleGeneral: // default
            default: 
                format = "%g";
                break;
        }
        ::sprintf(buffer, format, value);
        if (flags & fDoublePosix) {
            struct lconv* conv = localeconv();
            if ('.' != *(conv->decimal_point)) {
                char* pos = strchr(buffer, *(conv->decimal_point));
                if (pos) {
                    *pos = '.';
                }
            }
        }
    }
    out_str = buffer;
    errno = 0;
}


SIZE_TYPE DoubleToString(double value, unsigned int precision,
                               char* buf, SIZE_TYPE buf_size,
                               TNumToStringFlags flags)
{
    char buffer[kMaxDoubleStringSize];
    int n = 0;
    if ((flags & fDoublePosix) && (!isfinite(value) || value == 0.)) {
        if (value == 0.) {
            double zero = 0.;
            if (memcmp(&value, &zero, sizeof(double)) == 0) {
                strcpy(buffer, "0");
                n = 2;
            } else {
                strcpy(buffer, "-0");
                n = 3;
            }
        } else if (std::isnan(value)) {
            strcpy(buffer, "NaN");
            n = 4;
        } else if (value > 0.) {
            strcpy(buffer, "INF");
            n = 4;
        } else {
            strcpy(buffer, "-INF");
            n = 5;
        }
    } else {
        if (precision > (unsigned int)kMaxDoublePrecision) {
            precision = (unsigned int)kMaxDoublePrecision;
        }
        const char* format;
        switch (flags & fDoubleGeneral) {
            case fDoubleScientific:
                format = "%.*e";
                break;
            case fDoubleGeneral:
                format = "%.*g";
                break;
            case fDoubleFixed: // default
            default:
                format = "%.*f";
                break;
        }
        n = ::sprintf(buffer, format, (int)precision, value);
        if (n < 0) {
            n = 0;
        }
        if (flags & fDoublePosix) {
            struct lconv* conv = localeconv();
            if ('.' != *(conv->decimal_point)) {
                char* pos = strchr(buffer, *(conv->decimal_point));
                if (pos) {
                    *pos = '.';
                }
            }
        }
    }
    SIZE_TYPE n_copy = min((SIZE_TYPE) n, buf_size);
    memcpy(buf, buffer, n_copy);
    errno = 0;
    return n_copy;
}


char* s_ncbi_append_int2str(char* buffer, unsigned int value, size_t digits, bool zeros)
{
    char* buffer_start = buffer;
    char* buffer_end = (buffer += digits-1);
    if (zeros) {
        do {
            *buffer-- = (char)(48 + (value % 10));
            value /= 10;
        } while (--digits);
    } else {
        do {
            *buffer-- = (char)(48 + (value % 10));
        } while (value /= 10);

        if (++buffer != buffer_start) {
            memmove(buffer_start, buffer, buffer_end-buffer+1);
            buffer_end -= buffer - buffer_start;
        }
    }
    return ++buffer_end;
}


#define __NLG NCBI_CONST_LONGDOUBLE

SIZE_TYPE DoubleToString_Ecvt(double val, unsigned int precision,
                                    char* buffer, SIZE_TYPE bufsize,
                                    int* dec, int* sign)
{
    //errno = 0;
    *dec = *sign = 0;
    if (precision==0) {
        return 0;
    }
    if (precision > DBL_DIG) {
        precision = DBL_DIG;
    }
    if (val == 0.) {
        double zero = 0.;
        if (memcmp(&val, &zero, sizeof(double)) == 0) {
            *buffer='0';
            return 1;
        }
        *buffer='-';
        *(++buffer)='0';
        *sign = -1;
        return 2;
    }
    *sign = val < 0. ? -1 : 1;
    if (*sign < 0) {
        val = -val;
    }
    bool high_precision = precision > 9;

// calculate exponent
    unsigned int exp=0;
    bool exp_positive = val >= 1.;
    unsigned int first, second=0;
    long double mult = __NLG(1.);
    long double value = val;

    if (exp_positive) {
        while (value>=__NLG(1.e256))
            {value/=__NLG(1.e256); exp+=256;}
        if (value >= __NLG(1.e16)) {
            if      (value>=__NLG(1.e240)) {value*=__NLG(1.e-240); exp+=240;}
            else if (value>=__NLG(1.e224)) {value*=__NLG(1.e-224); exp+=224;}
            else if (value>=__NLG(1.e208)) {value*=__NLG(1.e-208); exp+=208;}
            else if (value>=__NLG(1.e192)) {value*=__NLG(1.e-192); exp+=192;}
            else if (value>=__NLG(1.e176)) {value*=__NLG(1.e-176); exp+=176;}
            else if (value>=__NLG(1.e160)) {value*=__NLG(1.e-160); exp+=160;}
            else if (value>=__NLG(1.e144)) {value*=__NLG(1.e-144); exp+=144;}
            else if (value>=__NLG(1.e128)) {value*=__NLG(1.e-128); exp+=128;}
            else if (value>=__NLG(1.e112)) {value*=__NLG(1.e-112); exp+=112;}
            else if (value>=__NLG(1.e96))  {value*=__NLG(1.e-96);  exp+=96;}
            else if (value>=__NLG(1.e80))  {value*=__NLG(1.e-80);  exp+=80;}
            else if (value>=__NLG(1.e64))  {value*=__NLG(1.e-64);  exp+=64;}
            else if (value>=__NLG(1.e48))  {value*=__NLG(1.e-48);  exp+=48;}
            else if (value>=__NLG(1.e32))  {value*=__NLG(1.e-32);  exp+=32;}
            else if (value>=__NLG(1.e16))  {value*=__NLG(1.e-16);  exp+=16;}
        }
        if      (value<   __NLG(1.)) {mult=__NLG(1.e+9); exp-= 1;}
        else if (value<  __NLG(10.)) {mult=__NLG(1.e+8);         }
        else if (value< __NLG(1.e2)) {mult=__NLG(1.e+7); exp+= 1;}
        else if (value< __NLG(1.e3)) {mult=__NLG(1.e+6); exp+= 2;}
        else if (value< __NLG(1.e4)) {mult=__NLG(1.e+5); exp+= 3;}
        else if (value< __NLG(1.e5)) {mult=__NLG(1.e+4); exp+= 4;}
        else if (value< __NLG(1.e6)) {mult=__NLG(1.e+3); exp+= 5;}
        else if (value< __NLG(1.e7)) {mult=__NLG(1.e+2); exp+= 6;}
        else if (value< __NLG(1.e8)) {mult=  __NLG(10.); exp+= 7;}
        else if (value< __NLG(1.e9)) {mult=   __NLG(1.); exp+= 8;}
        else if (value<__NLG(1.e10)) {mult=  __NLG(0.1); exp+= 9;}
        else if (value<__NLG(1.e11)) {mult=__NLG(1.e-2); exp+=10;}
        else if (value<__NLG(1.e12)) {mult=__NLG(1.e-3); exp+=11;}
        else if (value<__NLG(1.e13)) {mult=__NLG(1.e-4); exp+=12;}
        else if (value<__NLG(1.e14)) {mult=__NLG(1.e-5); exp+=13;}
        else if (value<__NLG(1.e15)) {mult=__NLG(1.e-6); exp+=14;}
        else if (value<__NLG(1.e16)) {mult=__NLG(1.e-7); exp+=15;}
        else                         {mult=__NLG(1.e-8); exp+=16;}
    } else {
        while (value<=__NLG(1.e-256))
            {value*=__NLG(1.e256); exp+=256;}
        if (value <= __NLG(1.e-16)) {
            if      (value<=__NLG(1.e-240)) {value*=__NLG(1.e240); exp+=240;}
            else if (value<=__NLG(1.e-224)) {value*=__NLG(1.e224); exp+=224;}
            else if (value<=__NLG(1.e-208)) {value*=__NLG(1.e208); exp+=208;}
            else if (value<=__NLG(1.e-192)) {value*=__NLG(1.e192); exp+=192;}
            else if (value<=__NLG(1.e-176)) {value*=__NLG(1.e176); exp+=176;}
            else if (value<=__NLG(1.e-160)) {value*=__NLG(1.e160); exp+=160;}
            else if (value<=__NLG(1.e-144)) {value*=__NLG(1.e144); exp+=144;}
            else if (value<=__NLG(1.e-128)) {value*=__NLG(1.e128); exp+=128;}
            else if (value<=__NLG(1.e-112)) {value*=__NLG(1.e112); exp+=112;}
            else if (value<=__NLG(1.e-96))  {value*=__NLG(1.e96);  exp+=96;}
            else if (value<=__NLG(1.e-80))  {value*=__NLG(1.e80);  exp+=80;}
            else if (value<=__NLG(1.e-64))  {value*=__NLG(1.e64);  exp+=64;}
            else if (value<=__NLG(1.e-48))  {value*=__NLG(1.e48);  exp+=48;}
            else if (value<=__NLG(1.e-32))  {value*=__NLG(1.e32);  exp+=32;}
            else if (value<=__NLG(1.e-16))  {value*=__NLG(1.e16);  exp+=16;}
        }
        if      (value<__NLG(1.e-15)) {mult=__NLG(1.e24); exp+=16;}
        else if (value<__NLG(1.e-14)) {mult=__NLG(1.e23); exp+=15;}
        else if (value<__NLG(1.e-13)) {mult=__NLG(1.e22); exp+=14;}
        else if (value<__NLG(1.e-12)) {mult=__NLG(1.e21); exp+=13;}
        else if (value<__NLG(1.e-11)) {mult=__NLG(1.e20); exp+=12;}
        else if (value<__NLG(1.e-10)) {mult=__NLG(1.e19); exp+=11;}
        else if (value<__NLG(1.e-9))  {mult=__NLG(1.e18); exp+=10;}
        else if (value<__NLG(1.e-8))  {mult=__NLG(1.e17); exp+=9;}
        else if (value<__NLG(1.e-7))  {mult=__NLG(1.e16); exp+=8;}
        else if (value<__NLG(1.e-6))  {mult=__NLG(1.e15); exp+=7;}
        else if (value<__NLG(1.e-5))  {mult=__NLG(1.e14); exp+=6;}
        else if (value<__NLG(1.e-4))  {mult=__NLG(1.e13); exp+=5;}
        else if (value<__NLG(1.e-3))  {mult=__NLG(1.e12); exp+=4;}
        else if (value<__NLG(1.e-2))  {mult=__NLG(1.e11); exp+=3;}
        else if (value<__NLG(1.e-1))  {mult=__NLG(1.e10); exp+=2;}
        else if (value<__NLG(1.))     {mult=__NLG(1.e9);  exp+=1;}
        else                          {mult=__NLG(1.e8);         }
    }

// get all digits
    long double t1 = value * mult;
    if (t1 >= __NLG(1.e9)) {
        first = 999999999;
    } else if (t1 < __NLG(1.e8)) {
        first = 100000000;
        t1 = first;
    } else {
        first = (unsigned int)t1;
    }
    if (high_precision) {
        long double t2 = (t1-first) * __NLG(1.e8);
        if (t2 >= __NLG(1.e8)) {
            second = 99999999;
        } else {
            second = (unsigned int)t2;
        }
    }

// convert them into string
    bool use_ext_buffer = bufsize > 20;
    char tmp[32];
    char *digits = use_ext_buffer ? buffer : tmp;
    char *digits_end = s_ncbi_append_int2str(digits,first,9,false);
    if (high_precision) {
        digits_end = s_ncbi_append_int2str(digits_end,second,8,true);
    }
    size_t digits_len = digits_end - digits;
    size_t digits_got = digits_len;
    size_t digits_expected = high_precision ? 17 : 9;

// get significant digits according to requested precision
    size_t pos = precision;
    if (digits_len > precision) {
        digits_len = precision;

        // this is questionable, but in fact,
        // improves the result (on average)
#if 1
        if (high_precision) {
            if (digits[pos] == '4') {
                size_t pt = pos-1;
                while (pt != 0 && digits[--pt] == '9')
                    ;
                if (pt != 0 && (pos-pt) > precision/2)
                    digits[pos]='5';
            } else if (digits[pos] == '5') {
                size_t pt = pos;
                while (pt != 0 && digits[--pt] == '0')
                    ;
                if (pt != 0 && (pos-pt) > precision/2)
                    digits[pos]='4';
            }
        }
#endif

        if (digits[pos] >= '5') {
            do {
                if (digits[--pos] < '9') {
                    ++digits[pos++];
                    break;
                }
                digits[pos]='0';
            } while (pos > 0);
            if (pos == 0) {
                if (digits_expected <= digits_got) {
                    if (exp_positive) {
                       ++exp; 
                    } else {
// exp cannot be 0, by design
                        exp_positive = --exp == 0;
                    }
                }
                *digits = '1';
                digits_len = 1;
            }
        }
    }

// truncate trailing zeros
    for (pos = digits_len; pos-- > 0 && digits[pos] == '0';)
        --digits_len;

    *dec = (int)exp;
    if (!exp_positive) {
        *dec = -*dec;
    }
    if (!use_ext_buffer) {
        if (digits_len <= bufsize) {
            strncpy(buffer,digits,digits_len);
        } else {
            //NCBI_THROW2(CStringException, eConvert,
            //            "Destination buffer too small", 0);
            cerr << "Destination buffer too small" << endl;
            errno = ERANGE;
            return 0;
        }
    }
    return digits_len;
}
#undef __NLG


SIZE_TYPE DoubleToStringPosix(double val, unsigned int precision,
                                    char* buffer, SIZE_TYPE bufsize)
{
    if (bufsize < precision+8) {
        //NCBI_THROW2(CStringException, eConvert,
        //            "Destination buffer too small", 0);
        errno = ERANGE;
        cerr << "Destination buffer too small" << endl;
        return 0;
    }
    int dec=0, sign=0;
    char digits[32];
    size_t digits_len = DoubleToString_Ecvt(
        val, precision, digits, sizeof(digits), &dec, &sign);
    if (digits_len == 0) {
        errno = 0;
        return 0;
    }
    if (val == 0.) {
        strncpy(buffer,digits, digits_len);
        return digits_len;
    }
    if (digits_len == 1 && dec == 0 && sign >=0) {
        *buffer = digits[0];
        errno = 0;
        return 1;
    }
    bool exp_positive = dec >= 0;
    unsigned int exp= (unsigned int)(exp_positive ? dec : (-dec));

    // assemble the result
    char *buffer_pos = buffer;
//    char *buffer_end = buffer + bufsize;
    char *digits_pos = digits;

    if (sign < 0) {
        *buffer_pos++ = '-';
    }
    // The 'e' format is used when the exponent of the value is less than -4
    //  or greater than or equal to the precision argument
    if ((exp_positive && exp >= precision) || (!exp_positive && exp > 4)) {
        *buffer_pos++ = *digits_pos++;
        --digits_len;
        if (digits_len != 0) {
            *buffer_pos++ = '.';
            strncpy(buffer_pos,digits_pos,digits_len);
            buffer_pos += digits_len;
        }
        *buffer_pos++ = 'e';
        *buffer_pos++ = exp_positive ? '+' : '-';

//#if defined(NCBI_OS_MSWIN)
#if NCBI_COMPILER_MSVC && _MSC_VER < 1900
        bool need_zeros = true;
        size_t need_digits = 3;
#else
        bool need_zeros = exp < 10 ? true : false;
        size_t need_digits = exp < 100 ? 2 : 3;
#endif
        // assuming exp < 1000
        buffer_pos = s_ncbi_append_int2str(buffer_pos, exp, need_digits,need_zeros);
    } else if (exp_positive) {
        *buffer_pos++ = *digits_pos++;
        --digits_len;
        if (digits_len > exp) {
            strncpy(buffer_pos,digits_pos,exp);
            buffer_pos += exp;
            *buffer_pos++ = '.';
            strncpy(buffer_pos,digits_pos+exp,digits_len-exp);
            buffer_pos += digits_len-exp;
        } else {
            strncpy(buffer_pos,digits_pos,digits_len);
            buffer_pos += digits_len;
            exp -= (unsigned int)digits_len;
            while (exp--) {
                *buffer_pos++ = '0';
            }
        }
    } else {
        *buffer_pos++ = '0';
        *buffer_pos++ = '.';
        for (--exp; exp--;) {
            *buffer_pos++ = '0';
        }
        strncpy(buffer_pos,digits_pos, digits_len);
        buffer_pos += digits_len;
    }
    errno = 0;
    return buffer_pos - buffer;
}

static const char* s_kTrueString  = "true";
static const char* s_kFalseString = "false";
static const char* s_kTString     = "t";
static const char* s_kFString     = "f";
static const char* s_kYesString   = "yes";
static const char* s_kNoString    = "no";
static const char* s_kYString     = "y";
static const char* s_kNString     = "n";


const string BoolToString(bool value)
{
    return value ? s_kTrueString : s_kFalseString;
}


bool StringToBool(const CTempString str)
{
    if ( AStrEquiv(str, s_kTrueString,  PNocase())  ||
         AStrEquiv(str, s_kTString,     PNocase())  ||
         AStrEquiv(str, s_kYesString,   PNocase())  ||
         AStrEquiv(str, s_kYString,     PNocase()) ) {
        errno = 0;
        return true;
    }
    if ( AStrEquiv(str, s_kFalseString, PNocase())  ||
         AStrEquiv(str, s_kFString,     PNocase())  ||
         AStrEquiv(str, s_kNoString,    PNocase())  ||
         AStrEquiv(str, s_kNString,     PNocase()) ) {
        errno = 0;
        return false;
    }
    HBN_ERR("String '%s' cannot be converted to bool", str.data());
}

END_HBNSTR_SCOPE