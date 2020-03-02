#include "../../corelib/seqdb.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/seqdb_summary.h"
#include "../../corelib/db_format.h"

#include <sys/stat.h>

static const HbnAppInfo sAppInfo = {
    "necat2viewdb",
    "Application to view HBN database summary information"
};

typedef enum {
    eEmpty,
    eFasta,
    eFastq,
    eSeqDB,
    eContext,
    eViewDBTypeError,
} EViewDBType;

EViewDBType
guess_view_db_type(int argc, char* argv[], 
    const char** seqdb_dir,
    const char** seqdb_title,
    const char** fasta,
    const char** fastq,
    const char** context,
    int* min_seq_size)
{
    if (argc < 2) return eViewDBTypeError;

    if (argc >= 3) {
        if (argc > 4) {
            HBN_LOG("too many arguments");
            return eViewDBTypeError;
        }
        if (!string_is_valid_number(argv[2])) {
            *seqdb_dir = argv[1];
            *seqdb_title = argv[2];
            if (argc == 4) {
                if (!string_is_valid_number(argv[3])) {
                    HBN_LOG("'%s' seems not a plausible number", argv[3]);
                    return eViewDBTypeError;
                }
                *min_seq_size = atoi(argv[3]);
            }
            return eSeqDB;
        }

        if (!string_is_valid_number(argv[2])) {
            HBN_LOG("'%s' seems not a plausible number", argv[3]);
            return eViewDBTypeError;            
        }
        *min_seq_size = atoi(argv[2]);
    }

    EDbFormat fmt = hbn_guess_db_format(argv[1]);
    if (fmt == eDbFormatEmptyFile) return eEmpty;
    if (fmt == eDbFormatFasta) {
        *fasta = argv[1];
        return eFasta;
    }
    if (fmt == eDbFormatFastq) {
        *fastq = argv[1];
        return eFastq;
    }
    return eContext;
}

static void
print_usage(const char* prog)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");

    fprintf(out, "\n");
    fprintf(out, "For HBNDB:\n");
    fprintf(out, "%s hbndb_dir hbndb_title [min_seq_len]\n", prog);

    fprintf(out, "\n");
    fprintf(out, "For FASTA:\n");
    fprintf(out, "%s fasta [min_seq_len]\n", prog);

    fprintf(out, "\n");
    fprintf(out, "For FASTQ:\n");
    fprintf(out, "%s fastq [min_seq_len]\n", prog);

    fprintf(out, "\n");
    fprintf(out, "For context_file:\n");
    fprintf(out, "%s context_file [min_seq_len]\n", prog);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION:\n");
    hbn_dump_app_info(out, &sAppInfo);
}

int main(int argc, char* argv[])
{
    CSeqDBSummary* summary = CSeqDBSummaryBuild(argc, argv);
    if (!summary) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    CSeqDBSummaryView(summary);
    free(summary);

    return EXIT_SUCCESS;
}