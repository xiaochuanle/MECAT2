#include "../../corelib/seqdb.h"
#include "../../corelib/build_db.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/ksort.h"

#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>

#define seqinfo_size_gt(lhs, rhs) ((lhs).seq_size > (rhs).seq_size)
KSORT_INIT(seqinfo_size_gt, CSeqInfo, seqinfo_size_gt);

static void
print_usage(const char* prog)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s genomesize coverage fasta/fastq/filelist\n", prog);
}

static void
make_random_seqdb_dir_name(char path[])
{
    sprintf(path, "ext_seqs_db_");
    char buf[11];
    int i = 0;
    srand((long)time(0));
    for (; i < 5; ++i) {
        int n = rand() % 26;
        char c = 'a' + n;
        buf[i] = c;
    }
    for (; i < 10; ++i) {
        int n = rand() % 10;
        char c = '0' + n;
        buf[i] = c;
    }
    buf[i] = '\0';
    strcat(path, buf);
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }
    const size_t kGenomeSize = datasize_to_u64(argv[1]);
    const int kCov = string_to_int(argv[2]);
    const char* input = argv[3];

    char seqdb_path[HBN_MAX_PATH_LEN];
    make_random_seqdb_dir_name(seqdb_path);
    if (mkdir(seqdb_path, S_IRWXU) != 0) HBN_ERR("Failed to create directory %s: %s", seqdb_path, strerror(errno));
    
    const size_t kVolSize = datasize_to_u64("4g");
    build_db(input,
        seqdb_path,
        INIT_QUERY_DB_TITLE,
        0,
        INT32_MAX,
        kVolSize,
        FALSE);

    const int num_seqs = seqdb_load_num_reads(seqdb_path, INIT_QUERY_DB_TITLE);
    CSeqInfo* seqinfo_array = load_seq_infos(seqdb_path, INIT_QUERY_DB_TITLE, 0, num_seqs);
    ks_introsort_seqinfo_size_gt(num_seqs, seqinfo_array);
    const size_t kTargetRes = kGenomeSize * kCov;
    size_t curr = 0;
    int min_size = 0;
    for (int i = 0; i < num_seqs; ++i) {
        curr += seqinfo_array[i].seq_size;
        if (curr >= kTargetRes) {
            min_size = seqinfo_array[i].seq_size;
            break;
        }
    }
    free(seqinfo_array);

    const int num_vols = seqdb_load_num_volumes(seqdb_path, INIT_QUERY_DB_TITLE);
    kv_dinit(vec_u8, seq);
    int ext_seqs = 0;
    size_t ext_res = 0;
    for (int vid = 0; vid < num_vols; ++vid) {
        CSeqDB* vol = seqdb_load(seqdb_path, INIT_QUERY_DB_TITLE, vid);
        for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
            if (vol->seq_info_list[i].seq_size < min_size) continue;
            ++ext_seqs;
            ext_res += vol->seq_info_list[i].seq_size;
            seqdb_extract_sequence(vol, i, FWD, &seq);
            for (size_t p = 0; p < kv_size(seq); ++p) {
                u8 c = kv_A(seq, p);
                kv_A(seq, p) = DECODE_RESIDUE(c);
            }
            fprintf(stdout, ">%s\n", seqdb_seq_name(vol, i));
            hbn_fwrite(kv_data(seq), 1, kv_size(seq), stdout);
            fprintf(stdout, "\n");
        }
        CSeqDBFree(vol);
    }
    kv_destroy(seq);

    char cmd[HBN_MAX_PATH_LEN];
    sprintf(cmd, "rm -r %s", seqdb_path);
    hbn_system(cmd);

    const int kWidth = 20;
    print_fixed_width_string(stderr, "genomeSize:", kWidth);
    u64_to_string_datasize(kGenomeSize, cmd);
    fprintf(stderr, "%s\n", cmd);

    print_fixed_width_string(stderr, "coverage:", kWidth);
    u64_to_string_comma(kCov, cmd);
    fprintf(stderr, "%sx\n", cmd);

    print_fixed_width_string(stderr, "minimum length:", kWidth);
    u64_to_string_comma(min_size, cmd);
    fprintf(stderr, "%s\n", cmd);

    print_fixed_width_string(stderr, "# seqs:", kWidth);
    u64_to_string_comma(ext_seqs, cmd);
    fprintf(stderr, "%s\n", cmd);

    print_fixed_width_string(stderr, "# res:", kWidth);
    u64_to_string_datasize(ext_res, cmd);
    fprintf(stderr, "%s\n", cmd);

    return 0;
}