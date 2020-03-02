#include "hbn_traceback_aux.h"

double
calc_ident_perc(const char* query_mapped_string, 
				const char* target_mapped_string,
			    const int align_size,
                int* dist,
				int* score)
{
	if (align_size == 0) return 0.0;
	
	int n = 0;
	for (int i = 0; i < align_size; ++i) {
		if (query_mapped_string[i] == target_mapped_string[i]) ++n;
	}
    if (dist) *dist = align_size - n;
	if (score) *score = align_size * MATCH_REWARD + n * MISMATCH_PENALTY;
	return 100.0 * n / align_size;
}

void
validate_aligned_string(const char* source_file,
						const char* source_func,
						const int source_line,
						int qid,
						const u8* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const u8* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	//return;
	//fwrite(query_mapped_string, 1, align_size, stderr);
	//fprintf(stderr, "\n");
	//fwrite(target_mapped_string, 1, align_size, stderr);
	//fprintf(stderr, "\n");

	//fprintf(stderr, "validating [%s, %s, %d]\n", source_file, source_func, source_line);

	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
            int c = right_extend ? query[x] : query[-x];
			hbn_assert(c >= 0 && c < 4);
			const char qc1 = DECODE_RESIDUE(c);
            if (qc != qc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu\n",
					  source_file,
					  source_func,
					  source_line,
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
                abort();
            }		  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
            int c = right_extend ? target[y] : target[-y];
			const char tc1 = DECODE_RESIDUE(c);
            if (tc != tc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d\n",
						  source_func,
						  source_func,
						  source_line,
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);
                abort();
            }
			++y;
		}
	}

	//fprintf(stderr, "validating [%s, %s, %d] done.\n", source_file, source_func, source_line);
}