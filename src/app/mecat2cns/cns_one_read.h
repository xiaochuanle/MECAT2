#ifndef __CNS_ONE_READ_H
#define __CNS_ONE_READ_H

#include "../../corelib/hbn_aux.h"
#include "cns_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
consensus_one_read(CnsThreadData* data, const int raw_read_id);

#ifdef __cplusplus
}
#endif

#endif // __CNS_ONE_READ_H