#ifndef __DIFF_GAPALIGN_H
#define __DIFF_GAPALIGN_H

#include "../corelib/hbn_aux.h"
#include "ksw2_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int segment_size;
    int row_size;
    int column_size;
    int segment_align_size;
    int d_path_size;
    int aln_path_size;
} DiffAlignParams;

void
DiffAlignParamsInit(DiffAlignParams* params);

typedef struct {
    int aln_str_size;
    int dist;
    int aln_q_s;
    int aln_q_e;
    int aln_t_s;
    int aln_t_e;
    char* q_aln_str;
    char* t_aln_str;
} DiffAlignment;

DiffAlignment*
DiffAlignmentNew(DiffAlignParams* params);

DiffAlignment*
DiffAlignmentFree(DiffAlignment* align);

void
DiffAlignmentClear(DiffAlignment* align);

typedef struct {
    int pre_k, x1, y1, x2, y2;
} DPathData;

typedef struct {
    int d, k, pre_k, x1, y1, x2, y2;
} DPathData2;

typedef struct {
    int x, y;
} PathPoint;

typedef struct {
    int qid, qoff, qend, qsize;
    int sid, soff, send, ssize;
    int score;
    int dist;
    double ident_perc;
    kstring_t qabuf;
    kstring_t sabuf;
    char* qas;
    char* qae;
    char* sas;
    char* sae;
    DiffAlignParams params;
    int* dynq;
    int* dynt;
    DiffAlignment* align;
    DPathData2* d_path;
    PathPoint*  aln_path;
    Ksw2Data* ksw;
} DiffGapAlignData;

#define DiffGapAlignDataDump(output_func, out, data) \
    output_func(out, "[%d, %d, %d] x [%d, %d, %d], %g\n", \
        (data)->qoff, (data)->qend, (data)->qsize, \
        (data)->soff, (data)->send, (data)->ssize, \
        (data)->ident_perc)

DiffGapAlignData*
DiffGapAlignDataNew();

DiffGapAlignData*
DiffGapAlignDataFree(DiffGapAlignData* data);

void
DiffGapAlignDataInit(DiffGapAlignData* data,
    int qoff,
    int qsize,
    int soff,
    int ssize);

int
diff_align(DiffGapAlignData* data,
    const u8* query,
    int qoff,
    int qsize,
    const u8* target,
    int toff,
    int tsize,
    const int min_align_size,
    const double min_ident_perc,
    const BOOL process_over_hang);

#ifdef __cplusplus
}
#endif

#endif // __DIFF_GAPALIGN_H