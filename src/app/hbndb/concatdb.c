#include "../../corelib/seqdb.h"
#include "../../corelib/fasta.h"
#include "../../corelib/hbn_package_version.h"

#include <unistd.h>
#include <sys/stat.h>

static const HbnAppInfo sAppInfo = {
    "necat2concatdb",
    "Application to concatenate subjects in database into one sequence"
};

static void
print_usage(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s database\n", pn);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION:\n");
    hbn_dump_app_info(out, &sAppInfo);
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    FILE* out = stdout;
    fprintf(out, ">max_cns_reference\n");
    HbnFastaReader* reader = HbnFastaReaderNew(argv[1]);
    while (!HbnLineReaderAtEof(reader->line_reader)) {
        if (!HbnFastaReaderReadOneSeq(reader)) continue;
        hbn_fwrite(ks_s(reader->sequence), 1, ks_size(reader->sequence), out);
    }
    fprintf(out, "\n");
    HbnFastaReaderFree(reader);
}