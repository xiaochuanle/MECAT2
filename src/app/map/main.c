#include "hbn_job_control.h"
#include "cmdline_args.h"
#include "hbn_build_seqdb.h"
#include "hbn_task_struct.h"
#include "hbn_find_subseq_hit.h"
#include "hbn_align_one_volume.h"
#include "mecat_results.h"
#include "../../corelib/hbn_package_version.h"

static void
merge_qi_vs_sj_results(const HbnProgramOptions* opts,
    const int qi,
    const int sj,
    FILE* out)
{
    char path[HBN_MAX_PATH_LEN];
    make_qi_vs_sj_results_path(opts->db_dir, kBackupAlignResultsDir, qi, sj, path);
    hbn_dfopen(in, path, "r");
    const int kBufLen = 2048;
    char buffer[kBufLen];
    while (1) {
        int r = fread(buffer, 1, kBufLen, in);
        fwrite(buffer, 1, r, out);
        if (r < kBufLen) break;
    }
    hbn_fclose(in);
}

static void
merge_all_vs_sj_results(const HbnProgramOptions* opts,
    const int qid_start,
    const int num_query_vols,
    const int sj,
    FILE* out)
{
    for (int qi = qid_start; qi < num_query_vols; qi += opts->num_nodes) {
        merge_qi_vs_sj_results(opts, qi, sj, out);
    }
}

int main(int argc, char* argv[])
{
    HbnProgramOptions* opts = (HbnProgramOptions*)calloc(1, sizeof(HbnProgramOptions));
    ParseHbnProgramCmdLineArguments(argc, argv, opts);
    hbn_build_seqdb(opts, INIT_QUERY_DB_TITLE, INIT_SUBJECT_DB_TITLE);

    char path[HBN_MAX_PATH_LEN];
    sprintf(path, "%s/%s", opts->db_dir, kBackupAlignResultsDir);
    if ((access(path, F_OK) != 0)
        &&
        (mkdir(path, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", path, strerror(errno));
    }

    hbn_task_struct* task_struct = hbn_task_struct_new(opts);
    if (opts->outfmt == eSAM) {
        print_sam_prolog(task_struct->out, 
            kSamVersion, 
            HBN_PACKAGE_VERSION, 
            argc, 
            argv);
    }
    const int num_query_vols = seqdb_load_num_volumes(opts->db_dir, task_struct->query_db_title);
    const int num_subject_vols = seqdb_load_num_volumes(opts->db_dir, task_struct->subject_db_title);
    const int query_vol_stride = opts->num_nodes;
    const int subject_vol_stride = 1;
    char job_name[256];

    for (int svid = 0; svid < num_subject_vols; svid += subject_vol_stride) {
        HBN_LOG("Searching against S%s", u64_to_fixed_width_string(svid, HBN_DIGIT_WIDTH));
        int qvid = (task_struct->query_and_subject_are_the_same ? svid : 0) + opts->node_id;
        if (all_vs_sj_is_mapped(opts->db_dir, 
                kBackupAlignResultsDir,
                qvid,
                num_query_vols,
                svid,
                opts->num_nodes)) {
            merge_all_vs_sj_results(task_struct->opts, qvid, num_query_vols, svid, task_struct->out);
            continue;
        }
        
        hbn_task_struct_build_subject_vol_context(task_struct, svid);
        for (; qvid < num_query_vols; qvid += query_vol_stride) {
            if (qi_vs_sj_is_mapped(opts->db_dir, kBackupAlignResultsDir, qvid, svid)) {
                merge_qi_vs_sj_results(task_struct->opts, qvid, svid, task_struct->out);
                continue;
            }
            char qibuf[64], sjbuf[64];
            u64_to_fixed_width_string_r(qvid, qibuf, HBN_DIGIT_WIDTH);
            u64_to_fixed_width_string_r(svid, sjbuf, HBN_DIGIT_WIDTH);
            sprintf(job_name, "Q%s_vs_S%s", qibuf, sjbuf);
            hbn_timing_begin(job_name);
            hbn_task_struct_build_query_vol_context(task_struct, qvid);
            hbn_align_one_volume(task_struct);
            qi_vs_sj_make_mapped(opts->db_dir, kBackupAlignResultsDir, qvid, svid);
            hbn_timing_end(job_name);
        }
    }

    task_struct = hbn_task_struct_free(task_struct);
    opts = HbnProgramOptionsFree(opts);
    return 0;
}