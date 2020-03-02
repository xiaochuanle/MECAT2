#ifndef __HBN_PACKAGE_VERSION_H
#define __HBN_PACKAGE_VERSION_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define HBN_PACKAGE                         1
#define HBN_PACKAGE_NAME                    "hs-blastn"
#define HBN_PACKAGE_VERSION_MAJOR           2
#define HBN_PACKAGE_VERSION_MINOR           0
#define HBN_PACKAGE_VERSION_PATCH           0
#define HBN_PACKAGE_CONFIG                  ""

#define HBN_PACKAGE_VERSION_STRINGIFY(x)    #x
#define HBN_PACKAGE_VERSION_COMPOSE_STR(a, b, c)    \
    HBN_PACKAGE_VERSION_STRINGIFY(a) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(b) "."            \
    HBN_PACKAGE_VERSION_STRINGIFY(c)

#define HBN_PACKAGE_VERSION             \
    HBN_PACKAGE_VERSION_COMPOSE_STR     \
    (                                   \
        HBN_PACKAGE_VERSION_MAJOR,      \
        HBN_PACKAGE_VERSION_MINOR,      \
        HBN_PACKAGE_VERSION_PATCH       \
    )

typedef struct {
    const char* app_name;
    const char* app_desc;
} HbnAppInfo;

void
hbn_build_info(char build_info []);

void
hbn_dump_app_info(FILE* stream, const HbnAppInfo* app);

#ifdef __cplusplus
}
#endif

#endif // __HBN_PACKAGE_VERSION_H