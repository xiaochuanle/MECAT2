#include "../../corelib/cmd_arg.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/gapped_candidate.h"
#include "../../corelib/partition_aux.h"
#include "../../corelib/seqdb.h"

#include <errno.h>
#include <sys/stat.h>

static const HbnAppInfo sAppInfo = {
    "mecat2pcan",
    "Application to partition consensus hits into small parts"
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
    fprintf(out, "%s [options] seqdb_dir pcan_dir cns_hits\n", prog);

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
    const char** pcan_dir,
    const char** cns_hits)
{
    *opts = init_pcan_opts;
    *seqdb_dir = NULL;
    *pcan_dir = NULL;
    *cns_hits = NULL;

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
    *pcan_dir = argv[i+1];
    *cns_hits = argv[i+2];

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

int cns_hit_qid(void* r)
{
    HbnConsensusInitHit* hit = (HbnConsensusInitHit*)(r);
    return hit->qid;
}

int cns_hit_sid(void* r)
{
    HbnConsensusInitHit* hit = (HbnConsensusInitHit*)(r);
    return hit->sid;
}

void cns_hit_change_roles(void* src, void* dst)
{
    HbnConsensusInitHit* s = (HbnConsensusInitHit*)(src);
    HbnConsensusInitHit* d = (HbnConsensusInitHit*)(dst);

    d->qid = s->sid;
    d->qoff = s->soff;
    d->sid = s->qid;
    d->soff = s->qoff;
    d->score = s->score;
    d->strand = s->strand;
}

void cns_hit_sort_by_sid(size_t n, void* a)
{
    HbnConsensusInitHit* hit_array = (HbnConsensusInitHit*)(a);
    ks_introsort_cns_hit_sid_lt(n, hit_array);
}

void cns_hit_normalise_sdir(void* r)
{
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
make_pcan_wrk_dir(const char* path)
{
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }
}

static BOOL
pcan_is_done(const char* pcan_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/pcan.done", pcan_dir);
    return access(path, F_OK) == 0;
}

static void
pcan_make_done(const char* pcan_dir)
{
    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/pcan.done", pcan_dir);
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);    
}

int main(int argc, char* argv[])
{
    const char* seqdb_dir = NULL;
    const char* pcan_dir = NULL;
    const char* cns_hits = NULL;
    PartCnsHitOptions opts;
    if (parse_cmd_args(argc, argv, &opts, &seqdb_dir, &pcan_dir, &cns_hits) == eCmdArgParseError) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }
    if (pcan_is_done(pcan_dir)) return 0;
    opts.part_count = fix_part_counts(opts.part_count);
    make_pcan_wrk_dir(pcan_dir);

    const int num_reads = seqdb_load_num_reads(seqdb_dir, INIT_QUERY_DB_TITLE);
    const int num_parts = (num_reads + opts.part_size - 1) / opts.part_size;
    dump_partition_count(pcan_dir, NULL, num_parts);
    part_record_main(pcan_dir,
        cns_hits,
        num_parts,
        opts.part_size,
        opts.num_threads,
        opts.part_count,
        sizeof(HbnConsensusInitHit),
        cns_hit_qid,
        cns_hit_sid,
        cns_hit_change_roles,
        cns_hit_normalise_sdir,
        cns_hit_sort_by_sid);

    pcan_make_done(pcan_dir);
    return EXIT_SUCCESS;
}