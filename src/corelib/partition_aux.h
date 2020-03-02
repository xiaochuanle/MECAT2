#ifndef __PARTITION_AUX_H
#define __PARTITION_AUX_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_PART_PREFIX "p"

void
make_partition_name(const char* data_dir, const char* prefix, const int pid, char path[]);

void
dump_partition_count(const char* data_dir, const char* prefix, const int np);

int
load_partition_count(const char* data_dir, const char* prefix);

typedef int (*qid_extract_func)(void* r);
typedef int (*sid_extract_func)(void* r);
typedef void (*change_record_roles_func)(void* src, void* dst);
typedef void (*normolize_sdir_func)(void* r);
typedef void (*record_sort_func)(size_t n, void* a);

void* load_part_records(const char* path, const size_t record_size, size_t* n_record);

void
part_record_main(const char* part_wrk_dir,
    const char* record_path,
    const int num_batches,
    const int batch_size,
    const int num_threads,
    const int num_dumpped_files,
    const size_t record_size,
    qid_extract_func get_qid,
    sid_extract_func get_sid,
    change_record_roles_func change_roles,
    normolize_sdir_func normalise_sdir,
    record_sort_func sort_records);

#ifdef __cplusplus
}
#endif

#endif // __PARTITION_AUX_H