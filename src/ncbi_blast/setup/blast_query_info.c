#include "blast_query_info.h"

BlastQueryInfo*
BlastQueryInfoNew(int max_num_contexts)
{
    BlastQueryInfo* query_info = (BlastQueryInfo*)calloc(1, sizeof(BlastQueryInfo));
    if (max_num_contexts > 0) {
        query_info->contexts = (BlastContextInfo*)calloc(max_num_contexts, sizeof(BlastContextInfo));
    }
    return query_info;
}

BlastQueryInfo*
BlastQueryInfoFree(BlastQueryInfo* query_info)
{
    if (query_info->contexts) free(query_info->contexts);
    free(query_info);
    return NULL;
}