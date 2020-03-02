#ifndef __SEQ_TAG_H
#define __SEQ_TAG_H

#include "name2id_map.h"
#include "gapped_candidate.h"
#include "string2hsp.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SEQ_TAG_PERFECT         0
#define SEQ_TAG_UNPERFECT       1
#define SEQ_TAG_REPEAT          2
#define SEQ_TAG_CHIMERIC        3
#define SEQ_TAG_CIRCULAR        4
#define SEQ_TAG_TANDEM_REPEAT   5
#define SEQ_TAG_RCDUP           6
#define SEQ_TAG_MAX             7

typedef int seq_tag_t;

typedef struct {
    int qid;
    int qsize;
    int sid;
    double ident_perc;
    int tag;
    size_t qoff, qend;
    size_t soff, send;
} SeqTag;

const char* GetSeqTagName(int tag);
int SeqNameToTag(const char* name);

seq_tag_t TagSeqFromHSPArray(HbnHSP* hsp_array, int hsp_count, double min_ident_perc, SeqTag* seq_tag);

BOOL query_tag_is_valid(const SeqTag* tag_array, const HbnConsensusInitHit* hit);

BOOL examine_ovlp_quality(const SeqTag* tag_array,
    int qoff,
    int qend,
    int soff,
    int send,
    double min_ovlp_frac,
    const int read_id,
    const int read_size,
    const int subject_id,
    const int subject_size);

#ifdef __cplusplus
}
#endif

#endif // __SEQ_TAG_H