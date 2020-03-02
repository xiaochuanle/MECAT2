#ifndef __C_NCBI_BLAST_AUX_H
#define __C_NCBI_BLAST_AUX_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NCBI_XBLAST_EXPORT
#define NCBI_TABLES_EXPORT

#define HBN_MASKS_THESE_CODES 1

/** Codons are always of length 3 */
#ifndef CODON_LENGTH
#define CODON_LENGTH 3
#endif

/** For translated gapped searches, this is the default value in
 * nucleotides of longest_intron (for ungapped translated searches,
 * the default value of longest_intron is zero, which causes a legacy
 * method of HSP linking that does not use longest_intron to be
 * invoked).
 *
 * The value 122 corresponds to 40 amino acids: 40 codons * 3
 * nucleotides per codon + up to 2 frame shifts.  40 amino acids is
 * the maximum gap size in the untranslated sequence, so
 * DEFAULT_LONGEST_INTRON makes these two gap sizes equal.
 */ 
#ifndef DEFAULT_LONGEST_INTRON
#define DEFAULT_LONGEST_INTRON 122
#endif

/** Compression ratio of nucleotide bases (4 bases in 1 byte) */
#ifndef COMPRESSION_RATIO
#define COMPRESSION_RATIO 4
#endif

/** Number of frames to which we translate in translating searches */
#ifndef NUM_FRAMES
#define NUM_FRAMES 6
#endif

/** Number of frames in a nucleotide sequence */
#ifndef NUM_STRANDS
#define NUM_STRANDS 2
#endif

/** Length of the genetic code string */
#ifndef GENCODE_STRLEN
#define GENCODE_STRLEN 64
#endif

#ifndef ABS
/** returns absolute value of a (|a|) */
#define ABS(a)	((a)>=0?(a):-(a))
#endif

#ifndef SIGN
/** return +1 for a > 0, -1 for a < 0 */
#define SIGN(a)	((a)>0?1:((a)<0?-1:0))
#endif

#ifndef NCBIMATH_LN2
/** natural log of 2. */
#define NCBIMATH_LN2      0.69314718055994530941723212145818
#endif

#ifndef DIM
/** dimension of an array. */
#define DIM(A) (sizeof(A)/sizeof((A)[0]))
#endif

#ifndef NULLB
/** terminating byte of a char* string. */
#define NULLB '\0'
#endif

/** Maximal unpacked subject sequence length for which full translation is 
 * performed up front. 
 */
#define MAX_FULL_TRANSLATION 2100

/** This sentry value is used as a 'fence' around the valid portions
 * of partially decoded sequences.  If an alignment finds this value
 * in a subject sequence, the fence_hit flag should be used to request
 * a refetch of the whole sequence, and the alignment restarted.
 * @note this value is repeated in seqdbgeneral.hpp
 */
#define FENCE_SENTRY 201

/** Safe free a pointer: belongs to a higher level header. */
#ifndef sfree
#define sfree(x) __sfree((void**)(void*)&(x))
#endif

/** Implemented in blast_util.c. @sa sfree */
NCBI_XBLAST_EXPORT
void __sfree(void** x);

/**
 * A macro expression that returns 1, 0, -1 if a is greater than,
 * equal to or less than b, respectively.  This macro evaluates its
 * arguments more than once.
 */
#ifndef BLAST_CMP
#define BLAST_CMP(a,b) ((a)>(b) ? 1 : ((a)<(b) ? -1 : 0))
#endif

/** TRUE if c is between a and b; f between d and e.  Determines if the
 * coordinates are already in an HSP that has been evaluated. 
*/
#define CONTAINED_IN_HSP(a,b,c,d,e,f) \
    (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)

#ifndef ASSERT
#define ASSERT hbn_assert
#endif

#ifndef MAX
#define MAX hbn_max
#define MIN hbn_min
#endif

#ifndef INT1_MAX
#define INT1_MAX I8_MAX
#define INT2_MIN I16_MIN
#define INT2_MAX I16_MAX
#define INT4_MIN I32_MIN
#define INT4_MAX I32_MAX
#endif

/****************************** Constants *********************************/

NCBI_XBLAST_EXPORT extern const int kDustLevel;  /**< Level parameter used by dust. */
NCBI_XBLAST_EXPORT extern const int kDustWindow; /**< Window parameter used by dust. */
NCBI_XBLAST_EXPORT extern const int kDustLinker; /**< Parameter used by dust to link together close low-complexity segments. */

NCBI_XBLAST_EXPORT extern const int kSegWindow;  /**< Window that SEG examines at once. */
NCBI_XBLAST_EXPORT extern const double kSegLocut;   /**< Locut parameter for SEG. */
NCBI_XBLAST_EXPORT extern const double kSegHicut;   /**< Hicut parameter for SEG. */

/** Enumeration for the stages in the BLAST search */
typedef enum EBlastStage {
    /** None specified */
    eNone               = 0x0,
    /** Preliminary stage */
    ePrelimSearch       = 0x1 << 0,
    /** Traceback stage */
    eTracebackSearch    = 0x1 << 1,
    /** Both preliminary and traceback stages */
    eBoth               = (ePrelimSearch | eTracebackSearch)
} EBlastStage;

/* NOTE: Please keep these comments in sync with argument descriptions in
 * CCompositionBasedStatsArgs::SetArgumentDescriptions()
 */

/** An collection of constants that specify all permissible
 * modes of composition adjustment */
typedef enum ECompoAdjustModes {
    /** Don't use composition based statistics */
    eNoCompositionBasedStats       = 0, 
    /** Composition-based statistics as in NAR 29:2994-3005, 2001 */
    eCompositionBasedStats         = 1, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, conditioned on sequence properties. Cannot be applied to PSSMs. */
    eCompositionMatrixAdjust       = 2, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, unconditionally. Cannot be applied to PSSMs. */
    eCompoForceFullMatrixAdjust    = 3,
    eNumCompoAdjustModes
} ECompoAdjustModes;

/** Copies memory using memcpy and malloc
 * @param orig memory to be copied [in]
 * @param size amount to be copied [in]
 * @return pointer to newly allocated memory. NULL if orig NULL, size is zero,
 *   or allocation fails.
 */
NCBI_XBLAST_EXPORT
void* BlastMemDup (const void *orig, size_t size);


/******************************************************************************/

/** A generic linked list node structure */
typedef struct ListNode {
	Uint1 choice;   /**< to pick a choice */
	void *ptr;              /**< attached data */
	struct ListNode *next;  /**< next in linked list */
} ListNode;

/** Create a new list node  
 * @param vnp Pointer to the start of the list, may be NULL [in]
 * @return newly allocated node 
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeNew (ListNode* vnp);

/** Add a node to the list.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @return New node
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeAdd (ListNode** head);

/** Add a node to the list with a given choice and data pointer.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @param choice Choice value for the new node. [in]
 * @param value Data pointer for the new node. [in]
 * @return New node
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeAddPointer (ListNode** head, Uint1 choice, void *value);

/** Free all list's nodes, does not attempt to free data. 
 * @param vnp objects to be freed [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeFree (ListNode* vnp);

/** Free nodes as well as data (vnp->ptr) assuming it is one contiguous chunk.
 * @param vnp objects to be freed [in] 
 * @return NULL
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeFreeData (ListNode* vnp);

/** Add a node to the list with a provided choice, and attached data 
 * pointing to a provided string.
 * @param head Pointer to the start of the list, if *head is NULL will
 *  be Pointer to new node. [in] [out]
 * @param choice sets the "choice" field in ListNode [in]
 * @param str char* buffer to be copied [in]
 * @return newly allocated node 
 */
NCBI_XBLAST_EXPORT
ListNode* ListNodeCopyStr (ListNode** head, Uint1 choice, const char* str);

typedef int TNCBIScore; /* historically signed char */
typedef struct SNCBIPackedScoreMatrix {
    const char*       symbols;  /**< order of residues */
    const TNCBIScore* scores;   /**< strlen(symbols) x strlen(symbols) */
    TNCBIScore        defscore; /**< score for unknown residues */
} SNCBIPackedScoreMatrix;

/** Returns a copy of the input string with all its characters turned to
 * uppercase. Useful for saving score matrix names. Caller is responsible for
 * deallocating return value.
 * @param string string to copy [in]
 * @return newly allocated string in upper case or NULL if string is NULL or
 * out of memory
 */
NCBI_XBLAST_EXPORT
char*
BLAST_StrToUpper(const char* string);

/// Default value for overhang
#define kBestHit_OverhangDflt 0.1
/// Minimum value for overhang
#define kBestHit_OverhangMin 0.0
/// Maximum value for overhang
#define kBestHit_OverhangMax 0.5

/// Default value for score_edge
#define kBestHit_ScoreEdgeDflt 0.1
/// Minimum value for score_edge
#define kBestHit_ScoreEdgeMin  0.0
/// Maximum value for score_edge
#define kBestHit_ScoreEdgeMax  0.5

/****************************************************************************/
/* Matrix utility functions */

/** Generic 2 dimensional matrix allocator.
 * Allocates a ncols by nrows matrix with cells of size data_type_sz. Must be
 * freed using x_DeallocateMatrix
 * @param   ncols number of columns in matrix [in]
 * @param   nrows number of rows in matrix [in]
 * @param   data_type_sz size of the data type (in bytes) to allocate for each
 *          element in the matrix [in]
 * @return pointer to allocated memory or NULL in case of failure
 */
NCBI_XBLAST_EXPORT 
void**
_PSIAllocateMatrix(unsigned int ncols, unsigned int nrows, 
                   unsigned int data_type_sz);

/** Generic 2 dimensional matrix deallocator.
 * Deallocates the memory allocated by x_AllocateMatrix
 * @param matrix matrix to deallocate   [in]
 * @param ncols number of columns in the matrix [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT 
void**
_PSIDeallocateMatrix(void** matrix, unsigned int ncols);

/** Copies src matrix into dest matrix, both of which must be int matrices with
 * dimensions ncols by nrows
 * @param dest Destination matrix           [out]
 * @param src Source matrix                 [in]
 * @param ncols Number of columns to copy   [in]
 * @param nrows Number of rows to copy      [in]
 */
NCBI_XBLAST_EXPORT 
void
_PSICopyMatrix_int(int** dest, int** src,
                   unsigned int ncols, unsigned int nrows);

/** Copies src matrix into dest matrix, both of which must be double matrices 
 * with dimensions ncols by nrows
 * @param dest Destination matrix           [out]
 * @param src Source matrix                 [in]
 * @param ncols Number of columns to copy   [in]
 * @param nrows Number of rows to copy      [in]
 */
NCBI_XBLAST_EXPORT 
void
_PSICopyMatrix_double(double** dest, double** src,
                      unsigned int ncols, unsigned int nrows);

/////////////

/** Stores the frequency ratios along with their bit scale factor */
typedef struct SFreqRatios {

    /** The actual frequency ratios */
    double**   data;

    /** Used to multiply the values in the above matrix to obtain scores in bit
     * units */
    int        bit_scale_factor;

} SFreqRatios;

//////////////

/** Information about a single pattern occurence in the query. */
typedef struct SPHIPatternInfo {
    Int4 offset;  /**< Starting offset of this pattern occurrence. */
    Int4 length;  /**< Length of this pattern occurrence. */
} SPHIPatternInfo;

/** In PHI BLAST, structure containing information about all pattern 
 * occurrences in query.
 */
typedef struct SPHIQueryInfo {
    Int4 num_patterns;  /**< Number of pattern occurrences in query. */
    SPHIPatternInfo *occurrences; /**< Array of pattern occurrence information
                                        structures. */
    Int4 allocated_size; /**< Allocated size of the occurrences array. */
    double probability; /**< Estimated probability of the pattern */
    char* pattern;   /**< Pattern used, saved here for formatting purposes. */
} SPHIQueryInfo;

#ifdef __cplusplus
}
#endif

#endif // __C_NCBI_BLAST_AUX_H