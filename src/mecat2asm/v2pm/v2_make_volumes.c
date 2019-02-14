#include <stdio.h>

#include "../klib/kstring.h"

#include <assert.h>

void
print_usage(const char* prog)
{
	FILE* out = stderr;
	fprintf(out, "USAGE:\n");
	fprintf(out, "%s wrk_dir fasta_input\n", prog);
}

static FILE* open_file(const char* path, const char* mode)
{
	FILE* file = fopen(path, mode);
	if (!file) {
		fprintf(stderr, "Error: failed to open file '%s' with mode '%s'\n", path, mode);
		exit(1);
	}
	return file;
}

static void output_kstring(FILE* out, kstring_t* s)
{
	size_t i;
	for (i = 0; i < kstr_size(*s); ++i) {
		char c = kstr_A(*s, i);
		fprintf(out, "%c", c);
	}
	fprintf(out, "\n");
}

int main(int argc, char* argv[])
{
	if (argc != 3) {
		print_usage(argv[0]);
		return 1;
	}

	const char* wrk_dir = argv[1];
	const char* input = argv[2];
	int min_id = 1, max_id = 0;
	int num_reads = 0;
	int s = 0;
	int volume = 1;
	new_kstring(hdr);
	new_kstring(seq);
	char path[2048];
	const int vs = 2000000000;

	FILE* in = open_file(input, "r");
	sprintf(path, "%s/%06d.fasta", wrk_dir, volume);
	FILE* out = open_file(path, "w");
	sprintf(path, "%s/ovlprep", wrk_dir);
	FILE* vi_out = open_file(path, "w");
	while (1) {
		kstr_clear(hdr);
		int r = kgetline(&hdr, fgets, in);
		if (r == EOF) break;
		kstr_clear(seq);
		r = kgetline(&seq, fgets, in);
		assert(r != EOF);
		int ss = (int)kstr_size(seq);
		if (s + ss > vs) {
			fprintf(vi_out, "-allreads -allbases -b %d -e %d\n", min_id, max_id);
			min_id = max_id + 1;
			fclose(out);
			++volume;
			sprintf(path, "%s/%06d.fasta", wrk_dir, volume);
			out = open_file(path, "w");
			s = 0;
		}
		output_kstring(out, &hdr);
		output_kstring(out, &seq);
		++num_reads;
		++max_id;
		s += ss;
	}

	if (s > 0) {
		fprintf(vi_out, "-allreads -allbases -b %d -e %d\n", min_id, max_id);
	}

	fclose(out);
	fclose(vi_out);
	fclose(in);
	free_kstring(hdr);
	free_kstring(seq);

	sprintf(path, "%s/num_reads.txt", wrk_dir);
	out = open_file(path, "w");
	fprintf(out, "%d\n", num_reads);
	fclose(out);
	sprintf(path, "%s/num_volumes.txt", wrk_dir);
	out = open_file(path, "w");
	fprintf(out, "%d\n", volume);
	fclose(out);
	sprintf(path, "%s/reads_info.txt", wrk_dir);
	out = open_file(path, "w");
	fprintf(out, "%d\t%d\n", volume, num_reads);
	fclose(out);
}
