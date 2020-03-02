#include "blast_sequence_blk.h"

BLAST_SequenceBlk*
BLAST_SequenceBlkNew()
{
    BLAST_SequenceBlk* seq_blk = (BLAST_SequenceBlk*)calloc(1, sizeof(BLAST_SequenceBlk));
    return seq_blk;
}

BLAST_SequenceBlk*
BLAST_SequenceBlkFree(BLAST_SequenceBlk* seq_blk)
{
    if (seq_blk->sequence) free(seq_blk->sequence);
    free(seq_blk);
    return NULL;
}