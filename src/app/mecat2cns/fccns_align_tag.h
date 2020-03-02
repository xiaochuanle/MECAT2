#ifndef __FCCNS_ALIGN_TAG_H
#define __FCCNS_ALIGN_TAG_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DEFAULT_CNS_WEIGHT  1.0

typedef u16 align_tag_delta_t;
#define ALIGN_TAG_MAX_DELTA U16_MAX

typedef struct {
    double weight;
    int t_pos;
    int p_t_pos;
    align_tag_delta_t delta;
    align_tag_delta_t p_delta;
    char q_base;
    char p_q_base;
} AlignTag;

#define align_tag_plink_eq(a, b) \
    ((a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base == (b).p_q_base)

typedef kvec_t(AlignTag) vec_align_tag;

void ks_introsort_align_tag_lt(size_t n, AlignTag* a);

void
make_align_tags_from_ovlp(const char* qaln,
    const char* taln,
    const size_t aln_size,
    const int qoff,
    const int qend,
    const int toff,
    const int tend,
    const double weight,
    vec_align_tag* align_tag_list);

#ifdef __cplusplus
}
#endif

#endif // __FCCNS_ALIGN_TAG_H