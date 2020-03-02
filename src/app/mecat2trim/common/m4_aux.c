#include "m4_aux.h"

M4Record*
load_partition_m4(const char* pm4_dir, const int pid, size_t* m4_count, vec_int* idx_range)
{
	char path[HBN_MAX_PATH_LEN];
	make_partition_name(pm4_dir, DEFAULT_PART_PREFIX, pid, path);
	M4Record* m4_array = load_part_records(path, sizeof(M4Record), m4_count);
	size_t nm4 = *m4_count;
	ks_introsort_m4_sid_lt(nm4, m4_array);

	int i = 0;
	kv_clear(*idx_range);
	while (i < nm4) {
		int j = i + 1;
		while (j < nm4 && m4_array[j].sid == m4_array[i].sid) ++j;
		kv_push(int, *idx_range, i);
		i = j;
	}
	kv_push(int, *idx_range, nm4);
	return m4_array;
}

M4Record*
get_next_range(M4Record* m4v, 
    pthread_mutex_t* range_get_lock, 
    int* idx_range, 
    int* next_range_id, 
    int nrange, 
    int* nm4)
{
	M4Record* m4s = NULL;
	int i = 0;
	pthread_mutex_lock(range_get_lock);
	i = (*next_range_id);
	++(*next_range_id);
	pthread_mutex_unlock(range_get_lock);
	
	if (i < nrange) {
		m4s = m4v + idx_range[i];
		*nm4 = idx_range[i+1] - idx_range[i];
	}
	return m4s;
}