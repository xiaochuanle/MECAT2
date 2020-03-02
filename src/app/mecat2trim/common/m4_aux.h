#ifndef __M4_AUX_H
#define __M4_AUX_H

#include "../../../corelib/m4_record.h"
#include "../../../corelib/partition_aux.h"

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

M4Record*
load_partition_m4(const char* pm4_dir, const int pid, size_t* m4_count, vec_int* idx_range);

M4Record*
get_next_range(M4Record* m4v, 
    pthread_mutex_t* range_get_lock, 
    int* idx_range, 
    int* next_range_id, 
    int nrange, 
    int* nm4);

#ifdef __cplusplus
}
#endif

#endif // __M4_AUX_H