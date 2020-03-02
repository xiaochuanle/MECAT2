#ifndef __FCCNS_H
#define __FCCNS_H

#include "fccns_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    vec_align_tag tag_list;
    vec_backbone_item item_list;
    vec_int cov_list;
    SmallObjectAlloc* li_alloc;
    SmallObjectAlloc* dci_alloc;
    int template_size;
} FCCnsData;

FCCnsData* FCCnsDataNew();

FCCnsData* FCCnsDataFree(FCCnsData* data);

void FCCnsDataClear(FCCnsData* data);

#ifdef __cplusplus
}
#endif

#endif // __FCCNS_H