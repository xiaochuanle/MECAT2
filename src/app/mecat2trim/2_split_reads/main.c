#include "split_reads.h"
#include "cmdline_args.h"
#include "../../../corelib/seqdb.h"
#include "../../../corelib/partition_aux.h"
#include "../../../corelib/line_reader.h"
#include "../../../ncbi_blast/c_ncbi_blast_aux.h"

static ClippedRange*
load_clipped_ranges(const char* lcr_path, const int num_seqs)
{
    ClippedRange* cr_array = (ClippedRange*)calloc(num_seqs, sizeof(ClippedRange));
    HbnLineReader* reader = HbnLineReaderNew(lcr_path);
    kstring_t* line = &reader->line;
    for (int i = 0; i < num_seqs; ++i) {
        hbn_assert(!HbnLineReaderAtEof(reader));
        HbnLineReaderReadOneLine(reader);
        kputc('\0', line);
        int id;
        HBN_SCANF(sscanf, ks_s(*line), 4, "%d%d%d%d", &id, &cr_array[i].left, &cr_array[i].right, &cr_array[i].size);
        hbn_assert(id == i);
    }
    HbnLineReaderFree(reader);
    return cr_array;
}

int main(int argc, char* argv[])
{
    HbnProgramOptions opts;
    ParseHbnProgramCmdLineArguments(argc, argv, &opts);
    const int num_seqs = seqdb_load_num_reads(opts.db_dir, INIT_QUERY_DB_TITLE);
    const int num_parts = load_partition_count(opts.pm4_dir, NULL);
    ClippedRange* cr_array = load_clipped_ranges(opts.lcr_path, num_seqs);
    ClippedRange* sr_array = (ClippedRange*)calloc(num_seqs, sizeof(ClippedRange));
    CSeqInfo* seqinfo_array = load_seq_infos(opts.db_dir, INIT_QUERY_DB_TITLE, 0, num_seqs);
    for (int i = 0; i < num_seqs; ++i) sr_array[i].size = seqinfo_array[i].seq_size;
    sfree(seqinfo_array);
    char job_name[256];
    char part_buf[64];
    for (int i = 0; i < num_parts; ++i) {
        u64_to_fixed_width_string_r(i, part_buf, HBN_DIGIT_WIDTH);
        sprintf(job_name, "Find clear ranges for part %s", part_buf);
        hbn_timing_begin(job_name);
        split_reads_for_one_partition(opts.pm4_dir, i, opts.min_size, cr_array, sr_array, opts.num_threads);
        hbn_timing_end(job_name);
    }

    hbn_dfopen(out, opts.output, "w");
    for (int i = 0; i < num_seqs; ++i) {
        ClippedRange* sr = sr_array + i;
        if (sr->right - sr->left < opts.min_size) {
            sr->left = 0;
            sr->right = 0;
            sr->size = 0;
        }
        fprintf(out, "%d\t%d\t%d\t%d\n", i, sr->left, sr->right, sr->size);
    }
    hbn_fclose(out);
    free(cr_array);
    free(sr_array);

    return 0;
}