#ifndef __FCCNS_AUX_H
#define __FCCNS_AUX_H

#include "../../corelib/small_object_alloc.h"

#include "fccns_align_tag.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double weight;
    int p_t_pos;
    align_tag_delta_t p_delta;
    char p_q_base;
    int link_count;
} LinkInfo;

typedef struct {
    int n_link;
    int coverage;
    LinkInfo* plinks;
    int best_p_t_pos;
    align_tag_delta_t best_p_delta;
    u8 best_p_q_base;
    double score;
} BaseLinks;

#define init_base_links(blnk) ( \
    (blnk).n_link = 0, \
    (blnk).coverage = 0, \
    (blnk).best_p_t_pos = -1, \
    (blnk).best_p_delta = ALIGN_TAG_MAX_DELTA, \
    (blnk).best_p_q_base = '.', \
    (blnk).score = .0 \
)

typedef struct {
    BaseLinks links[5];
} DeltaCovInfo;

typedef struct {
    int n_delta;
    DeltaCovInfo* delta;
} BackboneItem;

typedef kvec_t(BackboneItem) vec_backbone_item;

#define clear_backbone_item(item) ((item).n_delta = 0, (item).delta = NULL)

void
build_backbone(AlignTag* tag_array,
    const int tag_count,
    const int template_size,
    SmallObjectAlloc* dci_alloc,
    SmallObjectAlloc* li_alloc,
    vec_backbone_item* item_list,
    vec_int* cov_list);

void
consensus_backbone_segment(BackboneItem* backbone,
        int from,
        int to,
        int* coverage,
        kstring_t* cns_seq,
        int* cns_from,
        int* cns_to);

#ifdef __cplusplus
}
#endif

#endif // __FCCNS_AUX_H