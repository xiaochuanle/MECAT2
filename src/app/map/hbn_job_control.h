#ifndef __HBN_JOB_CONTROL_H
#define __HBN_JOB_CONTROL_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char* kBackupAlignResultsDir;

const char*
make_qi_vs_sj_results_path(const char* wrk_dir, const char* stage, const int qi, const int sj, char path[]);

FILE*
open_qi_vs_sj_results_file(const char* wrk_dir, const char* stage, const int qi, const int sj);

BOOL 
qi_vs_sj_is_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj);

void
qi_vs_sj_make_mapped(const char* wrk_dir, const char* stage, const int qi, const int sj);

BOOL
all_vs_sj_is_mapped(const char* wrk_dir, 
    const char* stage, 
    const int qid_start,
    const int num_query_vols,
    const int sj,
    const int num_nodes);

#ifdef __cplusplus
}
#endif

#endif // __HBN_JOB_CONTROL_H