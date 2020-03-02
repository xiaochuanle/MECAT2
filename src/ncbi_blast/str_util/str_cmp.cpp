#include "str_cmp.hpp"

#include "../../corelib/hbn_aux.h"

BEGIN_HBNSTR_SCOPE

int strcmp(const char* s1, const size_t n1, const char* s2, const size_t n2, ECase use_case)
{
    if (!n1) return n2 ? -1 : 0;
    if (!n2) return 1;

    size_t n = hbn_min(n1, n2);
    if (use_case == ECase::eCase) {
        int res = memcmp(s1, s2, n);
        if (res) return res;
        return (n1 == n2) ? 0 : (n1 > n2 ? 1 : -1);
    } else {
        const char* p1 = s1;
        const char* p2 = s2;
        while (n) {
            int r = (*p1 == *p2 || tolower(*p1) == tolower(*p2));
            if (!r) break;
            ++p1;
            ++p2;
            --n;
        }
        if (!n) return (n1 == n2) ? 0 : (n1 > n2 ? 1 : -1);
        if (*p1 == *p2) return 0;
        return tolower(*p1) - tolower(*p2);
    }
}

bool IsBlank(const CTempString str, SIZE_TYPE pos)
{
    SIZE_TYPE len = str.length();
    for (SIZE_TYPE idx = pos; idx < len; ++idx) {
        if (!isspace((unsigned char) str[idx])) {
            return false;
        }
    }
    return true;
}

bool IsLower(const CTempString str)
{
    SIZE_TYPE len = str.length();
    for (SIZE_TYPE i = 0; i < len; ++i) {
        if (isalpha((unsigned char)str[i])  &&  !islower((unsigned char)str[i])) {
            return false;
        }
    }
    return true;
}

bool IsUpper(const CTempString str)
{
    SIZE_TYPE len = str.length();
    for (SIZE_TYPE i = 0; i < len; ++i) {
        if (isalpha((unsigned char)str[i])  &&  !isupper((unsigned char)str[i])) {
            return false;
        }
    }
    return true;
}

string& ToLower(string& str)
{
    NON_CONST_ITERATE (string, it, str) {
        *it = (char)tolower((unsigned char)(*it));
    }
    return str;
}

char* ToLower(char* str)
{
    char* s;
    for (s = str;  *str;  str++) {
        *str = (char)tolower((unsigned char)(*str));
    }
    return s;
}

string& ToUpper(string& str)
{
    NON_CONST_ITERATE (string, it, str) {
        *it = (char)toupper((unsigned char)(*it));
    }
    return str;
}

char* ToUpper(char* str)
{
    char* s;
    for (s = str;  *str;  str++) {
        *str = (char)toupper((unsigned char)(*str));
    }
    return s;
}

int CompareCase(const CTempStringEx s1, const CTempStringEx s2)
{
    SIZE_TYPE n1 = s1.length();
    SIZE_TYPE n2 = s2.length();
    if ( !n1 ) {
        return n2 ? -1 : 0;
    }
    if ( !n2 ) {
        return 1;
    }
    if (int res = memcmp(s1.data(), s2.data(), min(n1, n2))) {
        return res;
    }
    return (n1 == n2) ? 0 : (n1 > n2 ? 1 : -1);
}

int CompareNocase(const CTempStringEx s1, const CTempStringEx s2)
{
    SIZE_TYPE n1 = s1.length();
    SIZE_TYPE n2 = s2.length();

    if ( !n1 ) {
        return n2 ? -1 : 0;
    }
    if ( !n2 ) {
        return 1;
    }
    SIZE_TYPE n = min(n1, n2);
    const char* p1 = s1.data();
    const char* p2 = s2.data();

    while (n  &&  (*p1 == *p2  ||  
                   tolower((unsigned char)(*p1)) == tolower((unsigned char)(*p2))) ) {
        p1++;  p2++;  n--;
    }
    if ( !n ) {
        return (n1 == n2) ? 0 : (n1 > n2 ? 1 : -1);
    }
    if (*p1 == *p2) {
        return 0;
    }
    return tolower((unsigned char)(*p1)) - tolower((unsigned char)(*p2));
}


END_HBNSTR_SCOPE