#include "cmdline_args.h"
#include "cns_one_part.h"
#include "../../corelib/partition_aux.h"

#include <errno.h>
#include <sys/stat.h>

static BOOL
partition_is_corrected(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".corrected");
    return access(path, F_OK) == 0;
}

static void
partition_make_corrected(const char* can_dir, const int pid)
{
    char path[HBN_MAX_PATH_LEN];
    make_partition_name(can_dir, DEFAULT_PART_PREFIX, pid, path);
    strcat(path, ".corrected");
    hbn_dfopen(out, path, "w");
    hbn_fclose(out);
}

int main(int argc, char* argv[])
{
    HbnProgramOptions opts;
    ParseHbnProgramCmdLineArguments(argc, argv, &opts);
    hbn_task_struct* ht_struct = hbn_task_struct_new(&opts);
    int num_parts = load_partition_count(opts.can_dir, NULL);
    char job_name[256];
    char pid_str[64];
    for (int i = opts.node_id; i < num_parts; i += opts.num_nodes) {
        if (partition_is_corrected(opts.can_dir, i)) continue;
        u64_to_fixed_width_string_r(i, pid_str, HBN_DIGIT_WIDTH);
        sprintf(job_name, "correcting part %s", pid_str);
        hbn_timing_begin(job_name);
        cns_one_part(ht_struct, i);
        hbn_timing_end(job_name);
        partition_make_corrected(opts.can_dir, i);
    }
    hbn_task_struct_free(ht_struct);
    return 0;
}