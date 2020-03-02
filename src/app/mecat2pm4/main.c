#include "../../corelib/partition_aux.h"
#include "../../corelib/m4_record.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/cmd_arg.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/stat.h>

static const HbnAppInfo sAppInfo = {
    "mecat2pm4",
    "Application to partition m4 records into small parts"
};

typedef struct {
    int part_size;
    int part_count;
    int num_threads;
} PartCnsHitOptions;

static const PartCnsHitOptions init_pcan_opts = {
    .part_size = 100000,
    .part_count = 100,
    .num_threads = 1,
};

static void 
print_usage(const char* prog)
{
    FILE* out = stderr;
    fprintf(out, "USAGE:\n");
    fprintf(out, "%s [options] seqdb_dir pm4_dir m4_path\n", prog);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION:\n");
    hbn_dump_app_info(out, &sAppInfo);

    fprintf(out, "\n");
    fprintf(out, "OPTIONAL ARGUMENTS:\n");

    fprintf(out, "  -p <Integer, >0>\n");
    fprintf(out, "    Number of reads in each part\n");
    fprintf(out, "    Default = '%d'\n", init_pcan_opts.part_size);

    fprintf(out, "  -k <Integer, >0>\n");
    fprintf(out, "    Number of parts processed at each round\n");
    fprintf(out, "    (If <0, then it will be set to system limit value)\n");
    fprintf(out, "    Default = '%d'\n", init_pcan_opts.part_count);

    fprintf(out, "  -t <Integer, >0>\n");
    fprintf(out, "    Number of CPU threads\n");
    fprintf(out, "    Default = '%d'\n", init_pcan_opts.num_threads);
}

ECmdArgParseStatus
parse_cmd_args(int argc,
    char* argv[],
    PartCnsHitOptions* opts,
    const char** seqdb_dir,
    const char** pm4_dir,
    const char** m4_path)
{
    *opts = init_pcan_opts;
    *seqdb_dir = NULL;
    *pm4_dir = NULL;
    *m4_path = NULL;

    int i = 1;
    int64_t ix;

    while (i < argc) {
        if (argv[i][0] != '-') break;

        if (strcmp(argv[i], "-p") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->part_size = ix;
            i += 2;
            continue;
        }    

        if (strcmp(argv[i], "-k") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->part_count = ix;
            i += 2;
            continue;
        }   

        if (strcmp(argv[i], "-t") == 0) {
            if (validate_cmd_arg_cnt(__func__, argv[i], argc, i, 1) == eCmdArgParseError) return eCmdArgParseError;
            if (parse_int_arg(__func__, argv[i], argv[i+1], &ix) == eCmdArgParseError) return eCmdArgParseError;
            opts->num_threads = ix;
            i += 2;
            continue;
        }   

        fprintf(stderr, "[%s] ERROR: unrecognised option: %s\n", __func__, argv[i]);
        return eCmdArgParseError;
    }

    if (i + 3 != argc) {
        fprintf(stderr, "[%s] ERROR: seqdb dir, pcan dir, and consensus hits must be provided\n", __func__);
        return eCmdArgParseError;
    }

    *seqdb_dir = argv[i];
    *pm4_dir = argv[i+1];
    *m4_path = argv[i+2];

    hbn_assert(opts->part_size > 0);
    hbn_assert(opts->part_count > 0);
    hbn_assert(opts->num_threads > 0);
    
    const int kMaxThreads = 8;
    if (opts->num_threads > kMaxThreads) {
        HBN_LOG("You specified %d CPU threads (-num_threads), we will use %d only.", opts->num_threads, kMaxThreads);
        opts->num_threads = kMaxThreads;
    }
    
    return eCmdArgParseSuccess;
}

int m4_qid(void* r)
{
    return ((M4Record*)(r))->qid;
}

int m4_sid(void* r)
{
    return ((M4Record*)(r))->sid;
}

void
change_m4_roles(void* src, void* dst)
{
    M4Record* min = (M4Record*)(src);
    M4Record* mout = (M4Record*)(dst);
    
    mout->qid = min->sid;
    mout->qdir = min->sdir;
    mout->qoff = min->soff;
    mout->qend = min->send;
    mout->qsize = min->ssize;
    mout->sid = min->qid;
    mout->sdir = min->qdir;
    mout->soff = min->qoff;
    mout->send = min->qend;
    mout->ssize = min->qsize;
    mout->ident_perc = min->ident_perc;
    mout->score = min->score;
}

void 
normalise_m4_sdir(void* r)
{
    M4Record* m = (M4Record*)(r);
    if (m->sdir == REV) {
        m->sdir = FWD;
        m->qdir = 1 - m->qdir;
    }
}

void
sort_m4_by_sid(size_t n, void* a)
{
    ks_introsort_m4_sid_lt(n, (M4Record*)(a));
}

int
fix_part_counts(int num_files) {
    const int kMaxAllowedNumFiles = sysconf(_SC_OPEN_MAX) - 8;
	if (num_files < 0) {
		num_files = kMaxAllowedNumFiles;
        HBN_LOG("Set value of argument '-k' to its maximal allowed value %d", num_files);
	} else {
        if (num_files > kMaxAllowedNumFiles) {
            HBN_LOG("The value (%d) of argument '-k' exceeds the system limits, we reset it to %d", 
                num_files, kMaxAllowedNumFiles);
            num_files = kMaxAllowedNumFiles;
        }
    }
	return num_files;
}

static void
make_pm4_wrk_dir(const char* path)
{
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

static BOOL
pm4_is_done(const char* pm4_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/pcan.done", pm4_dir);
    return access(path, F_OK) == 0;
}

static void
pm4_make_done(const char* pm4_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/pcan.done", pm4_dir);
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);    
}

int main(int argc, char* argv[])
{
    const char* seqdb_dir = NULL;
    const char* pm4_dir = NULL;
    const char* m4_path = NULL;
    PartCnsHitOptions opts;
    if (parse_cmd_args(argc, argv, &opts, &seqdb_dir, &pm4_dir, &m4_path) == eCmdArgParseError) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (pm4_is_done(pm4_dir)) return 0;
    opts.part_count = fix_part_counts(opts.part_count);
    make_pm4_wrk_dir(pm4_dir);

    int num_reads = seqdb_load_num_reads(seqdb_dir, INIT_QUERY_DB_TITLE);
    int num_batches = (num_reads + opts.part_size - 1) / opts.part_size;
    dump_partition_count(pm4_dir, NULL, num_batches);
    part_record_main(pm4_dir,
        m4_path,
        num_batches,
        opts.part_size,
        opts.num_threads,
        opts.part_count,
        sizeof(M4Record),
        m4_qid,
        m4_sid,
        change_m4_roles,
        normalise_m4_sdir,
        sort_m4_by_sid);
    
    pm4_make_done(pm4_dir);
    return 0;
}