#include "fccns.h"

FCCnsData*
FCCnsDataNew()
{
    FCCnsData* data = (FCCnsData*)calloc(1, sizeof(FCCnsData));
    kv_init(data->tag_list);
    kv_init(data->item_list);
    kv_init(data->cov_list);
    data->li_alloc = SmallObjectAllocNew(sizeof(LinkInfo));
    data->dci_alloc = SmallObjectAllocNew(sizeof(DeltaCovInfo));
    return data;
}

FCCnsData*
FCCnsDataFree(FCCnsData* data)
{
    kv_destroy(data->tag_list);
    kv_destroy(data->item_list);
    kv_destroy(data->cov_list);
    data->li_alloc = SmallObjectAllocFree(data->li_alloc);
    data->dci_alloc = SmallObjectAllocFree(data->dci_alloc);
    free(data);
    return NULL;
}

void
FCCnsDataClear(FCCnsData* data)
{
    kv_clear(data->tag_list);
    kv_clear(data->item_list);
    kv_clear(data->cov_list);
    SmallObjectAllocClear(data->li_alloc);
    SmallObjectAllocClear(data->dci_alloc);
}