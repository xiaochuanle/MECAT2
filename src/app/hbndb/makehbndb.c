#include "../../corelib/build_db.h"
#include "../../corelib/cmd_arg.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/hbn_package_version.h"

#include <sys/stat.h>

static const HbnAppInfo sAppInfo = {
    "necat2mkdb",
    "Application to create HBN database"
};

typedef struct {
    int min_seq_size;
    int rename_seq_name;
    int max_file_seqs;
    size_t max_file_res;
} HbnDBOptions;

static const HbnDBOptions 
sDefaultHbnDBOptions = {
    0,
    0,
    I32_MAX,
    4000000000
};

static void 
print_usage(const char* prog)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s [options] input hbndb_dir hbndb_title\n", prog);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION:\n");
    hbn_dump_app_info(out, &sAppInfo);

    fprintf(out, "\n");
    fprintf(out, "OPTIONAL ARGUMENTS:\n");

    fprintf(out, "  -min_size <Integer, >=0>\n");
    fprintf(out, "    Skip sequences shorter than this value\n");
    fprintf(out, "    Default = '%d'\n", sDefaultHbnDBOptions.min_seq_size);

    fprintf(out, "  -rename_seq\n");
    fprintf(out, "    Give each sequence a new header, start from 0\n");
    fprintf(out, "    Default: disabled\n");

    fprintf(out, "  -max_file_res <Integer, >=1>\n");
    fprintf(out, "    Maximum residues in one HBN database file\n");
    fprintf(out, "    Default = '%zu'\n", sDefaultHbnDBOptions.max_file_res);

    fprintf(out, "  -max_file_seqs <Integer, >=1s>\n");
    fprintf(out, "    Maximum number of sequences in one HBN database file\n");
    fprintf(out, "    Default = '%d'\n", sDefaultHbnDBOptions.max_file_seqs);
}

ECmdArgParseStatus
parse_cmd_args(int argc,
    char* argv[],
    HbnDBOptions* opts,
    const char** input,
    const char** seqdb_dir,
    const char** seqdb_title)
{
    *opts = sDefaultHbnDBOptions;
    *input = NULL;
    *seqdb_dir = NULL;
    *seqdb_title = NULL;

    int i = 1;
    int64_t ix;
    ECmdArgParseStatus status = eCmdArgParseSuccess;

    while (i < argc) {
        if (argv[i][0] != '-') break;
        if (strlen(argv[i]) == 1) break;

        if (strcmp(argv[i], "-min_size") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->min_seq_size = ix;
            i += 2;
            continue;
        }      

        if (strcmp(argv[i], "-rename_seq") == 0) {
            opts->rename_seq_name = 1;
            i += 1;
            continue;
        }   

        if (strcmp(argv[i], "-max_file_res") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->max_file_res = ix;
            i += 2;
            continue;
        }   

        if (strcmp(argv[i], "-max_file_seqs") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->max_file_seqs = ix;
            i += 2;
            continue;
        }   

        fprintf(stderr, "[%s] ERROR: unrecognised option: %s\n", __func__, argv[i]);
        return eCmdArgParseError;
    }
    if (status == eCmdArgParseExitNormally) return status;

    if (i + 3 != argc) {
        fprintf(stderr, "[%s] ERROR: sequence data, hbndb directory and hbndb title must be provided\n", __func__);
        return eCmdArgParseError;
    }
    *input = argv[i];
    *seqdb_dir = argv[i + 1];
    *seqdb_title = argv[i + 2];
    
    return eCmdArgParseSuccess;
}

#define SEQDB_BUILD_ARGS \
    input, \
    seqdb_dir, \
    seqdb_title, \
    opts.min_seq_size, \
    opts.max_file_seqs, \
    opts.max_file_res, \
    opts.rename_seq_name

int main(int argc, char* argv[])
{
    const char* input = NULL;
    const char* seqdb_dir = NULL;
    const char* seqdb_title = NULL;
    HbnDBOptions opts;
    ECmdArgParseStatus status = parse_cmd_args(argc, argv, &opts, &input, &seqdb_dir, &seqdb_title);
    if (status == eCmdArgParseExitNormally) return EXIT_SUCCESS;
    if (status == eCmdArgParseError) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (seqdb_dir) {
        if (access(seqdb_dir, F_OK) != 0) {
            if (mkdir(seqdb_dir , S_IRWXU) != 0) {
                HBN_ERR("failed to create foldef %s", seqdb_dir);
            }
        }
    }

    if (hbndb_is_built(SEQDB_BUILD_ARGS)) {
        HBN_LOG("HBN database is already built. Exit normally.");
        return EXIT_SUCCESS;
    }
    build_db(SEQDB_BUILD_ARGS);
    hbndb_make_built(SEQDB_BUILD_ARGS);

    return EXIT_SUCCESS;
}