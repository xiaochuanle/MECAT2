#ifndef __CNS_ONE_PART_H
#define __CNS_ONE_PART_H

#include "hbn_task_struct.h"

#ifdef __cplusplus
extern "C" {
#endif

void
cns_one_part(hbn_task_struct* ht_struct, const int pid);

#ifdef __cplusplus
}
#endif

#endif // __CNS_ONE_PART_H