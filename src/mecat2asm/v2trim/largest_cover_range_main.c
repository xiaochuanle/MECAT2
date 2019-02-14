#include "pm4_aux.h"
#include "largest_cover_range.h"
#include "range_list.h"
#include "../common/ontcns_aux.h"

#include <stdio.h>

void
print_usage(const char* prog)
{
	FILE* out = stdout;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s m4 packed_reads_dir error_cutoff min_ovlp_size min_cov min_read_size output num_threads\n", prog);
}

int main(int argc, char* argv[])
{
	if (argc != 9) {
		print_usage(argv[0]);
		return 1;
	}

	const char* m4_path = argv[1];
	const char* reads_dir = argv[2];
	const double min_ident_perc = 100.0 - 100.0 * atof(argv[3]);
	const int min_ovlp_size = atoi(argv[4]);
	const int min_cov = atoi(argv[5]);
	const int min_size = atoi(argv[6]);
	const char* output = argv[7];
	const int num_threads = atoi(argv[8]);

	int num_reads = load_num_reads(reads_dir);
	new_kvec(vec_ClippedRange, clipped_ranges);
	kv_resize(ClippedRange, clipped_ranges, num_reads);
	for (int i = 0; i < num_reads; ++i) {
		kv_A(clipped_ranges, i).left = 0;
		kv_A(clipped_ranges, i).right = 0;
		kv_A(clipped_ranges, i).size = 0;
	}
	
	int num_partitions = load_num_partitions(m4_path);
printf("number of partitions: %d\n", num_partitions);
	for (int i = 0; i < num_partitions; ++i) {
		DFOPEN(log_out, "lcr_log.txt", "a+");
		fprintf(log_out, "processing partition %d\n", i);
		FCLOSE(log_out);
		get_largest_cover_range_for_one_partition(m4_path, 
				i, 
				min_ident_perc,
				min_ovlp_size,
				min_cov,
				kv_data(clipped_ranges),
				num_threads);
	}
	
	DFOPEN(out, output, "w");
	for (int i = 0; i < num_reads; ++i) {
		ClippedRange range = kv_A(clipped_ranges, i);
		if (range.right - range.left < min_size) {
			range.left = 0;
			range.right = 0;
			range.size = 0;
		}
		fprintf(out, "%d\t%d\t%d\t%d\n", i, range.left, range.right, range.size);
	}
	free_kvec(clipped_ranges);
	FCLOSE(out);
}
