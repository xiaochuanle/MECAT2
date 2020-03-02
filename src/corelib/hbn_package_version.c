#include "hbn_package_version.h"

void hbn_build_info(char build_info [])
{
    sprintf(build_info, "%s %s", __DATE__, __TIME__);
}

static void
hbn_dump_package_info(FILE* stream, const char* prolog)
{
    if (prolog) fprintf(stream, "%s", prolog);
    char build_info[512];
    hbn_build_info(build_info);
    fprintf(stream, "Package: %s %s, build %s",
        HBN_PACKAGE_NAME,
        HBN_PACKAGE_VERSION,
        build_info);
}

void
hbn_dump_app_info(FILE* stream, const HbnAppInfo* app)
{
    fprintf(stream, "%s %s v%s", app->app_desc, app->app_name, HBN_PACKAGE_VERSION);
    fprintf(stream, " [");
    hbn_dump_package_info(stream, NULL);
    fprintf(stream, "]\n");
}