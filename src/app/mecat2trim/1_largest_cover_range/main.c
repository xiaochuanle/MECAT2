#include "largest_cover_range.h"
#include "cmdline_args.h"
#include "../../../corelib/seqdb.h"
#include "../../../corelib/partition_aux.h"
#include "../../../ncbi_blast/c_ncbi_blast_aux.h"

int main(int argc, char* argv[])
{
    HbnProgramOptions opts;
    ParseHbnProgramCmdLineArguments(argc, argv, &opts);
    const int num_seqs = seqdb_load_num_reads(opts.db_dir, INIT_QUERY_DB_TITLE);
    const int num_parts = load_partition_count(opts.pm4_dir, NULL);
    ClippedRange* cr_array = (ClippedRange*)calloc(num_seqs, sizeof(ClippedRange));
    CSeqInfo* seqinfo_array = load_seq_infos(opts.db_dir, INIT_QUERY_DB_TITLE, 0, num_seqs);
    for (int i = 0; i < num_seqs; ++i) cr_array[i].size = seqinfo_array[i].seq_size;
    sfree(seqinfo_array);
    char job_name[256];
    char part_buf[64];
    for (int i = 0; i < num_parts; ++i) {
        u64_to_fixed_width_string_r(i, part_buf, HBN_DIGIT_WIDTH);
        sprintf(job_name, "Find largest cover ranges for part %s", part_buf);
        hbn_timing_begin(job_name);
        get_largest_cover_range_for_one_partition(opts.pm4_dir,
            i,
            opts.perc_identity,
            opts.ovlp_cov_res,
            opts.min_cov,
            cr_array,
            opts.num_threads);
        hbn_timing_end(job_name);
    }

    hbn_dfopen(out, opts.output, "w");
    for (int i = 0; i < num_seqs; ++i) {
        ClippedRange* cr = cr_array + i;
        if (cr->right - cr->left < opts.min_size) {
            cr->left = 0;
            cr->right = 0;
            cr->size = 0;
        }
        fprintf(out, "%d\t%d\t%d\t%d\n", i, cr->left, cr->right, cr->size);
    }
    hbn_fclose(out);
    free(cr_array);

    return 0;
}