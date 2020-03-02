#include "hbn_aux.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>
#include <ctype.h>
#include <locale.h>

#include "ksort.h"

KSORT_INIT_GENERIC(i32)
KSORT_INIT_GENERIC(idx)
KSORT_INIT_GENERIC(u32)
KSORT_INIT_GENERIC(u64)
KSORT_INIT_GENERIC(uidx)
KSORT_INIT_GENERIC(size_t)
#define int_pair_lt(a, b) (((a).first < (b).first) || ((a).first == (b).first && (a).second < (b).second))
KSORT_INIT(int_pair, IntPair, int_pair_lt)

static const char* hbn_log_level_name[] = {
    [LOG_EMERG]     = "EMERG",
    [LOG_ALERT]     = "ALERT",
    [LOG_CRIT]      = "CRIT",
    [LOG_ERR]       = "ERROR",
    [LOG_WARNING]   = "WARNING",
    [LOG_NOTICE]    = "NOTICE",
    [LOG_INFO]      = "INFO",
    [LOG_DEBUG]     = "DEBUG",
};

void
what_is_time_now(char now[])
{
    time_t ltime;
    time(&ltime);
    ctime_r(&ltime, now);
    size_t n = strlen(now);
    now[n-1] = '\0';
}

static void 
_hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, va_list args)
{
    char buf[HBN_MAX_PATH_LEN + 1];
    char time_str[256];
    what_is_time_now(time_str);
    vsnprintf(buf, sizeof(buf), fmt, args);
    fprintf(stderr, "[%s %s:%s:%d] <%s> %s\n", time_str, HBN_LOG_ARGS_GENERIC, hbn_log_level_name[level], buf);
}

void
hbn_dump_message(HBN_LOG_PARAMS_GENERIC, const int level, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    _hbn_dump_message(HBN_LOG_ARGS_GENERIC, level, fmt, args);
    va_end(args);
}

////////////////////////////////////////////////////////////////////////////////////////////

int safe_gzread(HBN_LOG_PARAMS_GENERIC, gzFile stream, void* buf, unsigned int len)
{
    int ret = gzread(stream, buf, len);
    if (ret < 0) {
        int errnum = 0;
        const char* msg = gzerror(stream, &errnum);
        const char* why = (Z_ERRNO == errnum) ? strerror(errno) : msg;
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

int err_gzread(gzFile stream, void* buf, unsigned int len)
{
    return safe_gzread(HBN_LOG_ARGS_DEFAULT, stream, buf, len);
}

gzFile safe_gzopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    gzFile stream;
    if (strcmp(path, "-") == 0) {
        stream = gzdopen(fileno(strstr(mode, "r") ? stdin : stdout), mode);
        /* according to zlib.h, this is the only reason gzdopen can fail */
        if (!stream) HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': out of memory", path, mode);
        return stream;
    }
    if ((stream = gzopen(path, mode)) == 0) {
        const char* why = errno ? strerror(errno) : "out of memory";
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, why);
    }
    return stream;
}

int safe_gzclose(HBN_LOG_PARAMS_GENERIC, gzFile stream)
{
    int ret = gzclose(stream);
    if (Z_OK != ret) {
        const char* why = (Z_ERRNO == ret) ? strerror(errno) : zError(ret);
        HBN_ERR_GENERIC("%s", why);
    }
    return ret;
}

FILE* safe_fopen(HBN_LOG_PARAMS_GENERIC, const char* path, const char* mode)
{
    FILE* stream = 0;
    if (strcmp(path, "-") == 0) return strstr(mode, "r") ? stdin : stdout;

    if ((stream = fopen(path, mode)) == 0) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to open file '%s' with mode '%s': %s", path, mode, y);
    }
    return stream;
}

size_t safe_fwrite(HBN_LOG_PARAMS_GENERIC, const void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fwrite(buf, size, nmemb, stream);
    if (ret != nmemb) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

size_t safe_fread(HBN_LOG_PARAMS_GENERIC, void* buf, size_t size, size_t nmemb, FILE* stream)
{
    size_t ret = fread(buf, size, nmemb, stream);
    if (ret != nmemb) {
        const char* y = ferror(stream) ? strerror(errno) : "Unexpected end of file";
        HBN_ERR_GENERIC("%s", y);
    }
    return ret;
}

int safe_fclose(HBN_LOG_PARAMS_GENERIC, FILE* stream)
{
    int ret = fclose(stream);
    if (ret != 0) {
        HBN_ERR_GENERIC("%s", strerror(errno));
    }
    return ret;
}

void hbn_exception(const char* expr, HBN_LOG_PARAMS_GENERIC, const char* fmt, ...)
{
    fprintf(stderr, "Assertion Failed At '%s:%s:%d'\n", file, func, line);
    fprintf(stderr, "\tExpression: '%s'\n", expr);
    if (!fmt) return;
    fprintf(stderr, "Context Information:\n");
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    abort();
}

off_t hbn_get_file_size(HBN_LOG_PARAMS_GENERIC, const char* path)
{
    struct stat sbuf;
    if (stat(path, &sbuf) == -1) {
        const char* y = strerror(errno);
        HBN_ERR_GENERIC("fail to stat file '%s': %s", path, y);
    }
    return sbuf.st_size;
}

double hbn_time_diff(const struct timeval* begin, const struct timeval* end)
{
    double d = end->tv_sec - begin->tv_sec;
    d += 1.0 * (end->tv_usec - begin->tv_usec) / 1e6;
    return d;
}

u8 nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

u8 nst_nt16_table[256] = {
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, //15
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, // 31
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 15, 16, 16, // 47
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, // 63
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 0, 10, 1, 11, 16, 16, 2, 12, 16, 16, 7, 16, 6, 14, 16,
	16, 16, 4, 9, 3, 16, 13, 8, 16, 5, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
	16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
};

void i64_to_string_with_comma_r(i64 n, char str_n[])
{
    static int comma = '\0';
    char rebuf[64];
    char* p = &rebuf[ sizeof(rebuf)-1 ];
    int i = 0;
    if (comma == '\0') {
        struct lconv* lcp = localeconv();
        if (lcp != NULL) {
            if (lcp->thousands_sep != NULL && *lcp->thousands_sep != '\0') {
                comma = *lcp->thousands_sep;
            } else {
                comma = ',';
            }
        }
    }
    int n_is_nonnegative = 1;
    if (n < 0) {
        n_is_nonnegative = 0;
        n *= -1;
    }
    hbn_assert(n >= 0);
    *p = '\0';
    do {
        if (i && i % 3 == 0) *--p = comma;
        *--p = '0' + n % 10;
        n /= 10;
        ++i;
    } while (n);
    if (!n_is_nonnegative) *--p = '-';

    char* q = str_n;
    while (*p != '\0') {
        *q = *p;
        ++q;
        ++p;
    }
    *q = '\0';
}

char* i64_to_string_with_comma(i64 n)
{
    static char n_str[64];
    i64_to_string_with_comma_r(n, n_str);
    return n_str;
}

void u64_to_fixed_width_string_r(u64 n, char n_str[], const int width)
{
    int n_digit = 0;
    size_t m = n;

    do {
        ++n_digit;
        n /= 10;
    } while (n);
    hbn_assert(n_digit <= width);

    char* p = n_str;
    int i = n_digit;
    while (i < width) {
        *p = '0';
        ++p;
        ++i;
    }
    sprintf(p, "%zu", m);
}

char*
u64_to_fixed_width_string(u64 n, const int width)
{
    static char n_str[64];
    u64_to_fixed_width_string_r(n, n_str, width);
    return n_str;
}

void
print_digit_with_comma(FILE* out, size_t num)
{
    char buf[64];
    i64_to_string_with_comma_r(num, buf);
    fprintf(out, "%s", buf);
}

void
print_fixed_width_string(FILE* out, const char* s, const int width)
{
    int n = strlen(s);
    fprintf(out, "%s", s);
    for (int i = n; i < width; ++i) fprintf(out, " ");
}

int
string_is_valid_number(const char * arg)
{
    char state = 1;
    for (char c = *arg; c; c = *(++arg)){
        switch (c) {
            case 'e':
                if (state != 2 && state != 14) return 0;
                state = 9;
                break;
            case '.':
                if (state >> 2) return 0;
                state = 12;
                break;
            case '+':
            case '-':
                if (!(state & 1)) return 0;
                state--;
                break;
            default:
                if (c >= '0' && c <= '9'){
                    if (state >> 2 == 1) return 0;
                    state = "\2\4\012\016"[state >> 2];
                } else {
                    return 0;
                }
                break;
        }
    }
    return state & 2;
}

#ifdef WIN32
#include "windows.h"
#else
#include "unistd.h"
#endif

int hbn_get_cpu_count()
{
#if WIN32
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    int allCPUNum_    = sysInfo.dwNumberOfProcessors;
    return allCPUNum_;
#else //linux
    //获取当前系统的所有CPU核数，包含禁用的
    int allCPUNum_    = sysconf(_SC_NPROCESSORS_CONF);
    //获取当前系统的可用CPU核数
    int enableCPUNum_ = sysconf(_SC_NPROCESSORS_ONLN);
    return enableCPUNum_;
#endif
}