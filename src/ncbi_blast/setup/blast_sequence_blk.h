#ifndef __BLAST_SEQUENCE_BLK_H
#define __BLAST_SEQUENCE_BLK_H

#include "../c_ncbi_blast_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/** A structure containing two integers, used e.g. for locations for the 
 * lookup table.
 */
typedef struct SSeqRange {
   Int4 left;  /**< left endpoint of range (zero based) */
   Int4 right;  /**< right endpoint of range (zero based) */
} SSeqRange;

/** Used to hold a set of positions, mostly used for filtering. 
 * oid holds the index of the query sequence.
*/
typedef struct BlastSeqLoc {
        struct BlastSeqLoc *next;  /**< next in linked list */
        SSeqRange *ssr;            /**< location data on the sequence. */
} BlastSeqLoc;

/** Structure for keeping the query masking information */
typedef struct BlastMaskLoc {
   /** Total size of the BlastSeqLoc array below. This is always the number 
     of queries times the number of contexts. Note that in the case of 
     translated query searches, these locations must be provided in protein 
     coordinates to BLAST_MainSetUp.
     @sa BLAST_GetNumberOfContexts 
     @sa BlastMaskLocDNAToProtein
    */
   Int4 total_size; 

   /** Array of masked locations. 
     Every query is allocated the number of contexts associated with the 
     program. In the case of nucleotide searches, the strand(s) to search 
     dictatate which elements of the array for a given query are filled. For 
     translated searches, this should also be the same (by design) but the 
     C toolkit API does NOT implement this, it rather fills all elements 
     for all queries with masked locations in protein coordinates (if any). 
     The C++ API does follow the convention which populates each element, only
     if so dictated by the strand(s) to search for each query.
     @sa BLAST_GetNumberOfContexts
    */
   BlastSeqLoc** seqloc_array; 
} BlastMaskLoc;

/** Define the possible subject masking types */
typedef enum ESubjectMaskingType {
    eNoSubjMasking,
    eSoftSubjMasking,
    eHardSubjMasking
} ESubjectMaskingType;

/** Structure to hold a sequence. */
typedef struct BLAST_SequenceBlk {
   Uint1* sequence; /**< Sequence used for search (could be translation). */
   Uint1* sequence_start; /**< Start of sequence, usually one byte before 
                               sequence as that byte is a NULL sentinel byte.*/
   Int4 length;         /**< Length of sequence. */
   Int2 frame; /**< Frame of the query, needed for translated searches */
   Int2 subject_strand; /**< Strand of the subject sequence for translated searches. 
                          Uses the same values as ENa_strand. */
   Int4 oid; /**< The ordinal id of the current sequence */
   Boolean sequence_allocated; /**< TRUE if memory has been allocated for 
                                  sequence */
   Boolean sequence_start_allocated; /**< TRUE if memory has been allocated 
                                        for sequence_start */
   Uint1* sequence_start_nomask; /**< Query sequence without masking. */
   Uint1* sequence_nomask; /**< Start of query sequence without masking. */
   Boolean nomask_allocated; /**< If false the two above are just pointers to
                                   sequence and sequence_start. */
   Uint1* oof_sequence; /**< Mixed-frame protein representation of a
                             nucleotide sequence for out-of-frame alignment */
   Boolean oof_sequence_allocated; /**< TRUE if memory has been allocated 
                                        for oof_sequence */
   Uint1* compressed_nuc_seq; /**< 4-to-1 compressed version of sequence */
   Uint1* compressed_nuc_seq_start; /**< start of compressed_nuc_seq */
   BlastMaskLoc* lcase_mask; /**< Locations to be masked from operations on 
                                this sequence: lookup table for query; 
                                scanning for subject. */
   Boolean lcase_mask_allocated; /**< TRUE if memory has been allocated for 
                                    lcase_mask */
   Int4 chunk;  /**< Used for indexing only: the chunk number within the 
                     subject sequence. */
   Uint1 *gen_code_string;  /**< for nucleotide subject sequences (tblast[nx]),
                              the genetic code used to create a translated
                              protein sequence (NULL if not applicable). This
                              field is NOT owned by this data structure, it's
                              owned by the genetic code singleton. 
                              @sa gencode_singleton.h
                              */
   /* BEGIN: Data members needed for masking subjects from a BLAST database */
   SSeqRange* seq_ranges;   /**< Ranges of the sequence to search */
   Uint4 num_seq_ranges;    /**< Number of elements in seq_ranges */
   Boolean seq_ranges_allocated;   /**< TRUE if memory has been allocated for
                                      seq_ranges */
   ESubjectMaskingType mask_type;          /**< type of subject masking */
   /* END: Data members needed for masking subjects from a BLAST database */

   Uint1 bases_offset; /* Bases offset in first byte for SRA seq */

} BLAST_SequenceBlk;

BLAST_SequenceBlk*
BLAST_SequenceBlkNew();

BLAST_SequenceBlk*
BLAST_SequenceBlkFree(BLAST_SequenceBlk* seq_blk);

#ifdef __cplusplus
}
#endif

#endif // __BLAST_SEQUENCE_BLK_H