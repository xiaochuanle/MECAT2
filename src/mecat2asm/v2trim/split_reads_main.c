#include "pm4_aux.h"
#include "split_reads_aux.h"
#include "range_list.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s m4 packed_reads_dir clear_range min_read_size output num_threads\n", prog);
}

static void
load_clear_ranges(const char* path, const int num_reads, vec_ClippedRange* clear_ranges)
{
	char line[2048];
	ClippedRange clr;
	int id = 0;
	DFOPEN(in, path, "r");
	for (int i = 0; i < num_reads; ++i) {
		ontcns_getline(line, 2048, in);
		SAFE_SCANF(sscanf, line, 4, "%d%d%d%d", &id, &clr.left, &clr.right, &clr.size);
		kv_push(ClippedRange, *clear_ranges, clr);
	}
	FCLOSE(in);
}

int main(int argc, char* argv[])
{
	if (argc != 7) {
		print_usage(argv[0]);
		return 1;
	}

	const char* m4_path = argv[1];
	const char* reads_dir = argv[2];
	const char* clear_range_path = argv[3];
	const int min_size = atoi(argv[4]);
	const char* output = argv[5];
	const int num_threads = atoi(argv[6]);

	int num_reads = load_num_reads(reads_dir);
	new_kvec(vec_ClippedRange, clipped_ranges);
	load_clear_ranges(clear_range_path, num_reads, &clipped_ranges);
	new_kvec(vec_ClippedRange, split_ranges);
	kv_resize(ClippedRange, split_ranges, num_reads);
	for (int i = 0; i < num_reads; ++i) {
		kv_A(split_ranges, i).left = 0;
		kv_A(split_ranges, i).right = 0;
		kv_A(split_ranges, i).size = 0;
	}
	
	int num_partitions = load_num_partitions(m4_path);
	for (int i = 0; i < num_partitions; ++i) {
		DFOPEN(log_out, "split_reads_log.txt", "a+");
		fprintf(log_out, "processing partition %d\n", i);
		FCLOSE(log_out);
		split_reads_for_one_partition(m4_path, i, min_size, kv_data(clipped_ranges), kv_data(split_ranges), num_threads);
	}
	
	DFOPEN(out, output, "w");
	for (int i = 0; i < num_reads; ++i) {
		ClippedRange range = kv_A(split_ranges, i);
		if (range.right - range.left < min_size) {
			range.left = 0;
			range.right = 0;
			range.size = 0;
		}
		fprintf(out, "%d\t%d\t%d\t%d\n", i, range.left, range.right, range.size);
	}
	free_kvec(clipped_ranges);
	free_kvec(split_ranges);
	FCLOSE(out);
}
