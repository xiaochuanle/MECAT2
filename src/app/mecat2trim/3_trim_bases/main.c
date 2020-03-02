#include "../common/range_list.h"
#include "../../../corelib/fasta.h"
#include "../../../corelib/seqdb.h"

#include <assert.h>
#include <stdio.h>

static void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s db_dir sr_path use_new_numeric_header output\n", prog);
}

static ClippedRange*
load_clipped_ranges(const char* sr_path, const int num_seqs)
{
    ClippedRange* cr_array = (ClippedRange*)calloc(num_seqs, sizeof(ClippedRange));
    HbnLineReader* reader = HbnLineReaderNew(sr_path);
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

static void
dump_one_volume(CSeqDB* vol, ClippedRange* sr_range, FILE* out, int* id)
{
    kv_dinit(vec_u8, seq);
    for (int i = 0; i < vol->dbinfo.num_seqs; ++i) {
        int gid = i + vol->dbinfo.seq_start_id;
        if (sr_range[gid].size == 0) continue;
        hbn_assert(sr_range[gid].size == vol->seq_info_list[i].seq_size);
        seqdb_extract_subsequence(vol, i, sr_range[i].left, sr_range[i].right, FWD, &seq);
        for (size_t p = 0; p < kv_size(seq); ++p) {
            u8 c = kv_A(seq, p);
            kv_A(seq, p) = DECODE_RESIDUE(c);
        }
        if (id == NULL) {
            fprintf(out, ">%s [From:To:OrgSeqSize:TrimSeqSize] = [%d:%d:%d:%d]\n",
                seqdb_seq_name(vol, i),
                sr_range[i].left,
                sr_range[i].right,
                sr_range[i].size,
                sr_range[i].right - sr_range[i].left);
        } else {
            fprintf(out, ">%d %s [From:To:OrgSeqSize:TrimSeqSize] = [%d:%d:%d:%d]\n",
                *id,
                seqdb_seq_name(vol, i),
                sr_range[i].left,
                sr_range[i].right,
                sr_range[i].size,
                sr_range[i].right - sr_range[i].left);
            ++(*id);
        }
        hbn_fwrite(kv_data(seq), 1, kv_size(seq), out);
        fprintf(out, "\n");
    }
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		print_usage(argv[0]);
		return 1;
	}
	
    const char* db_dir = argv[1];
    const char* sr_path = argv[2];
    const int use_new_numeric_header = atoi(argv[3]);
    const char* output = argv[4];
	
    const int num_seqs = seqdb_load_num_reads(db_dir, INIT_QUERY_DB_TITLE);
    ClippedRange* sr_array = load_clipped_ranges(sr_path, num_seqs);
	const int num_vols = seqdb_load_num_volumes(db_dir, INIT_QUERY_DB_TITLE);
    hbn_dfopen(out, output, "w");
    int id_ = 1;
    int* id = (use_new_numeric_header) ? (&id_) : NULL;
    for (int i = 0; i < num_vols; ++i) {
        CSeqDB* vol = seqdb_load(db_dir, INIT_QUERY_DB_TITLE, i);
        dump_one_volume(vol, sr_array, out, id);
        CSeqDBFree(vol);
    }
    hbn_fclose(out);
    free(sr_array);
}
