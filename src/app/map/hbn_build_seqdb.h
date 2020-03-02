#ifndef __HBN_BUILD_SEQDB_H
#define __HBN_BUILD_SEQDB_H

#include "hbn_options.h"

#ifdef __cplusplus
extern "C" {
#endif

void 
hbn_build_seqdb(const HbnProgramOptions* opts, 
    const char* query_db_title, 
    const char* subject_db_title);

#ifdef __cplusplus
}
#endif

#endif // __HBN_BUILD_SEQDB_H