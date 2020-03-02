#ifndef __SEQDB_SUMMARY_H
#define __SEQDB_SUMMARY_H

#include "seqdb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DbLenStatN 11

typedef struct {
    size_t  ambig_subseq_count;
    size_t  ambig_res_count;
    size_t  seq_count;
    size_t  res_count;
    size_t  max_len;
    size_t  min_len;
    size_t  avg_len;
    size_t  med_len;
    IntPair len_stats[DbLenStatN];
} CSeqDBSummary;

CSeqDBSummary*
CSeqDBSummaryFree(CSeqDBSummary* summary);

void
CSeqDBSummaryView(const CSeqDBSummary* summary);

CSeqDBSummary*
CSeqDBSummaryBuild(int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

#endif // __SEQDB_SUMMARY_H