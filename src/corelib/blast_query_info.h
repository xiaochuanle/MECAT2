#ifndef BLAST_QUERY_INFO_H
#define BLAST_QUERY_INFO_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/** The context related information */
typedef struct BlastContextInfo {
    Int4 query_offset;      /**< Offset of this query, strand or frame in the
                               concatenated super-query. */
    Int4 query_length;      /**< Length of this query, strand or frame */
    Int8 eff_searchsp;      /**< Effective search space for this context. */
    Int4 length_adjustment; /**< Length adjustment for boundary conditions */
    Int4 query_index;       /**< Index of query (same for all frames) */
    Int1 frame;             /**< Frame number (-1, -2, -3, 0, 1, 2, or 3) */
    Boolean is_valid;       /**< Determine if this context is valid or not.
                              This field should be set only by the setup code
                              and read by subsequent stages of the BLAST search
                              */
    Int4 segment_flags;     /**< Flags describing segments for paired reads */
} BlastContextInfo;

/** Forward declaration of SPHIQueryInfo */
//struct SPHIQueryInfo;

#ifndef SPHIQueryInfo
typedef struct SPHIQueryInfo {
} SPHIQueryInfo;
#endif

/** The query related information 
 */
typedef struct BlastQueryInfo {
    Int4 first_context;  /**< Index of the first element of the context array */
    Int4 last_context;   /**< Index of the last element of the context array */
    int num_queries;     /**< Number of query sequences */
    BlastContextInfo * contexts; /**< Information per context */
    Uint4 max_length;    /**< Length of the longest among the concatenated
                            queries */
    Uint4 min_length;    /**< Length of the shortest among the concatenated
                            queries */
    struct SPHIQueryInfo* pattern_info; /**< Counts of PHI BLAST pattern
                                      occurrences, used in PHI BLAST only. */
} BlastQueryInfo;

BlastQueryInfo*
BlastQueryInfoNew(int max_num_contexts);

BlastQueryInfo*
BlastQueryInfoFree(BlastQueryInfo* query_info);

#ifdef __cplusplus
}
#endif
#endif // BLAST_QUERY_INFO_H