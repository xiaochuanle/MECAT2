/*******************************************************************************************
 i
 *  Dextractor: pullls requested info out of .bax.h5 files produced by Pacbio
 *
 *
 *  Author:  Martin Pippel
 *  Date  :  Dec 12, 2013
 *
 *  Author:  Gene Myers
 *  Date:    Jan 8, 2014, redesign of the modes of operation and flags, and also the
 *               logic for extraction in writeBaxReads
 *
 ********************************************************************************************/

#ifndef _BAX_H5
#define _BAX_H5

#include <hdf5.h>
#include "DB.h"

typedef struct
  { int     fastq;         // if non-zero get fastq quality values (obsolete)
    int     quivqv;        // if non-zero get quiver file
    int     arrow;         // if non-zero get arrow data

    hsize_t numMV;
    char   *movieName;     // name of movie

    hsize_t numBP;         // sum of all raw read lengths
    char   *baseCall;      // 8 streams that may be extracted dependent on flag settings
    char   *delQV;
    char   *delTag;
    char   *insQV;
    char   *mergeQV;
    char   *subQV;
    char   *fastQV;
    uint16 *pulseW;

    hsize_t numZMW;        // number of wells/holes
    int    *readLen;       // length of each read in events
    char   *holeType;      // Hole type, only SEQUENCING holes are extracted
    float  *snrVec;        // 4 floats per ZMW
    int     chan[4];       // Channel order (typically T, G, A, C => 3, 2, 0, 1)

    hsize_t numHQR;        // number of regions
    int    *regions;       // region information (5 ints per entry)

    int     delLimit;     //  The Del QV associated with N's in the Del Tag

  } BaxData;

typedef struct
  { int data_off;   //  Offset of stream data into vectors
    int zmw_off;    //  Offset of well-indexed data

    int well;
    int fpulse;
    int lpulse;
    int qv;
  } SubRead;

void initBaxData(BaxData *b, int fastq, int quivqv, int arrow);
void freeBaxData(BaxData *b);

int      getBaxData(BaxData *b, char *fname);
void     printBaxError(int ecode);
SubRead *nextSubread(BaxData *b, int prime);

#endif // _BAX_H5
