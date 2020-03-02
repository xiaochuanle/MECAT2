#ifndef __LCR_OPTIONS_H
#define __LCR_OPTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    const char* db_dir;
    const char* pm4_dir;
    const char* output;
    double perc_identity;
    int ovlp_cov_res;
    int min_cov;
    int min_size;
    int num_threads;
} HbnProgramOptions;

#ifdef __cplusplus
}
#endif

#endif // __LCR_OPTIONS_H