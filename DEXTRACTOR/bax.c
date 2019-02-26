/*******************************************************************************************
 *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>

#include <hdf5.h>
#include "DB.h"
#include "bax.h"

// Exception codes

#define CANNOT_OPEN_BAX_FILE   1
#define BAX_BASECALL_ERR       2
#define BAX_DEL_ERR            3
#define BAX_TAG_ERR            4
#define BAX_INS_ERR            5
#define BAX_MRG_ERR            6
#define BAX_SUB_ERR            7
#define BAX_QV_ERR             8
#define BAX_PULSE_ERR          9
#define BAX_NR_EVENTS_ERR     10
#define BAX_REGION_ERR        11
#define BAX_HOLESTATUS_ERR    12
#define BAX_MOVIENAME_ERR     13
#define BAX_SNR_ERR           14
#define BAX_CHIPSET_ERR       15

//  Initialize *the* BaxData structure

void initBaxData(BaxData *b, int fastq, int quivqv, int arrow)
{ b->fastq     = fastq;
  b->quivqv    = quivqv;
  b->arrow     = arrow;
  b->movieName = NULL;
  b->baseCall  = NULL;
  b->delQV     = NULL;
  b->delTag    = NULL;
  b->insQV     = NULL;
  b->mergeQV   = NULL;
  b->subQV     = NULL;
  b->fastQV    = NULL;
  b->pulseW    = NULL;
  b->readLen   = NULL;
  b->holeType  = NULL;
  b->snrVec    = NULL;
  b->regions   = NULL;
  b->delLimit  = 0;
}

//  Check if memory needed is above highwater mark, and if so allocate

static void ensureMovie(BaxData *b, hsize_t len)
{ static hsize_t mmax = 0;

  b->numMV = len;
  if (mmax < len)
    { mmax = 1.2*len + 500;
      b->movieName = (char *) Realloc(b->movieName, mmax+1, "Allocating movie name");
    }
}

//  Check if memory needed is above highwater mark, and if so allocate

static void ensureBases(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numBP = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->baseCall = (char *) Realloc(b->baseCall, smax, "Allocating basecall vector");
      if (b->fastq)
        b->fastQV = (char *) Realloc(b->fastQV, smax, "Allocating fastq vector");
      if (b->quivqv)
        { b->delQV   = (char *) Realloc(b->delQV, 5ll*smax, "Allocating 5 QV vectors");
          b->delTag  = b->delQV   + smax;
          b->insQV   = b->delTag  + smax;
          b->mergeQV = b->insQV   + smax;
          b->subQV   = b->mergeQV + smax;
        }
      if (b->arrow)
        b->pulseW = (uint16 *) Realloc(b->pulseW, 2ll*smax, "Allocating arrow vector");
    }
}

static void ensureZMW(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numZMW = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->holeType = (char *) Realloc(b->holeType, smax, "Allocating hole vector");
      b->readLen  = (int *) Realloc(b->readLen , smax * sizeof(int), "Allocating event vector");
      if (b->arrow)
        b->snrVec = (float *) Realloc(b->snrVec,(4ll*smax+1)*sizeof(float),"Allocating snr vector");
    }
}

static void ensureHQR(BaxData *b, hsize_t len)
{ static hsize_t smax = 0;

  b->numHQR = len;
  if (smax < len)
    { smax = 1.2*len + 10000;
      b->regions = (int *) Realloc(b->regions,(5ll*smax+1)*sizeof(int),"Allocating region vector");
    }
}

// Fetch the relevant contents of the current bax.h5 file and return the H5 file id.

static char  DNA_2_NUMBER[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 0, 0, 1, 0, 0, 0, 2,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
  };

int getBaxData(BaxData *b, char *fname)
{ hid_t   field_space;
  hid_t   field_set;
  hsize_t field_len[2];
  hid_t   file_id;
  herr_t  stat;
  int     ecode;
  hid_t   type;
  hid_t   attr;
  char   *name;

  H5Eset_auto(H5E_DEFAULT,0,0); // silence hdf5 error stack

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0)
    return (CANNOT_OPEN_BAX_FILE);

#ifdef DEBUG
  printf("PROCESSING %s, file_id: %d\n", baxFileName, file_id);
#endif

#define GET_SIZE(path,error)									\
  { ecode = error;										\
    if ((field_set = H5Dopen2(file_id, path, H5P_DEFAULT)) < 0) goto exit0;			\
    if ((field_space = H5Dget_space(field_set)) < 0) goto exit1;				\
    H5Sget_simple_extent_dims(field_space, field_len, NULL);					\
  }

#define FETCH(field,type)									\
  { stat = H5Dread(field_set, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, b->field);			\
    H5Sclose(field_space);									\
    H5Dclose(field_set);									\
    if (stat < 0) goto exit0;									\
  }

#define CHECK_FETCH(path,error,field,type,cntr)							\
  { GET_SIZE(path,error)									\
    if (b->cntr != field_len[0]) goto exit2;							\
    FETCH(field,type)										\
  }

  ecode = BAX_MOVIENAME_ERR;
  if ((field_set = H5Gopen2(file_id,"/ScanData/RunInfo",H5P_DEFAULT)) < 0) goto exit0;
  if ((attr = H5Aopen(field_set,"MovieName",H5P_DEFAULT)) < 0) goto exit3;
  if ((field_space = H5Aget_space(attr)) < 0) goto exit4;
  if ((type = H5Aget_type(attr)) < 0) goto exit5;

  H5Aread(attr,type,&name);
  ensureMovie(b,strlen(name));
  strcpy(b->movieName,name);

  H5Tclose(type);
  H5Sclose(field_space);
  H5Aclose(attr);
  H5Gclose(field_set);

  GET_SIZE("/PulseData/BaseCalls/Basecall",BAX_BASECALL_ERR)
  ensureBases(b,field_len[0]);
  FETCH(baseCall,H5T_NATIVE_UCHAR)
  if (b->arrow)
    CHECK_FETCH("/PulseData/BaseCalls/WidthInFrames",BAX_PULSE_ERR,pulseW,H5T_NATIVE_USHORT,numBP)
  if (b->fastq)
    CHECK_FETCH("/PulseData/BaseCalls/QualityValue",BAX_QV_ERR,fastQV,H5T_NATIVE_UCHAR,numBP)
  if (b->quivqv)
    { CHECK_FETCH("/PulseData/BaseCalls/DeletionQV",    BAX_DEL_ERR,delQV,  H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/DeletionTag",   BAX_TAG_ERR,delTag, H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/InsertionQV",   BAX_INS_ERR,insQV,  H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/MergeQV",       BAX_MRG_ERR,mergeQV,H5T_NATIVE_UCHAR,numBP)
      CHECK_FETCH("/PulseData/BaseCalls/SubstitutionQV",BAX_SUB_ERR,subQV,  H5T_NATIVE_UCHAR,numBP)
    }

  GET_SIZE("/PulseData/BaseCalls/ZMW/HoleStatus",BAX_HOLESTATUS_ERR)
  ensureZMW(b,field_len[0]);
  FETCH(holeType,H5T_NATIVE_UCHAR)
  CHECK_FETCH("/PulseData/BaseCalls/ZMW/NumEvent",BAX_NR_EVENTS_ERR,readLen,H5T_NATIVE_INT,numZMW)
  if (b->arrow)
    { CHECK_FETCH("/PulseData/BaseCalls/ZMWMetrics/HQRegionSNR",BAX_SNR_ERR,snrVec,
                  H5T_NATIVE_FLOAT,numZMW)

      ecode = BAX_CHIPSET_ERR;
      if ((field_set = H5Gopen2(file_id,"/ScanData/DyeSet",H5P_DEFAULT)) < 0) goto exit0;
      if ((attr = H5Aopen(field_set,"BaseMap",H5P_DEFAULT)) < 0) goto exit3;
      if ((field_space = H5Aget_space(attr)) < 0) goto exit4;
      if ((type = H5Aget_type(attr)) < 0) goto exit5;

      H5Aread(attr,type,&name);

      { int i;

        for (i = 0; i < 4; i++)
          b->chan[i] = DNA_2_NUMBER[(int) name[i]];
      }

      H5Tclose(type);
      H5Sclose(field_space);
      H5Aclose(attr);
      H5Gclose(field_set);

    }

  GET_SIZE("/PulseData/Regions",BAX_REGION_ERR)
  ensureHQR(b,field_len[0]);
  FETCH(regions,H5T_NATIVE_INT)

  //  Find the Del QV associated with N's in the Del Tag

  if (b->quivqv)
    { hsize_t  i;
  
      for (i = 0; i < b->numBP; i++)
        if (b->delTag[i] == 'N')
          { b->delLimit = b->delQV[i];
            break;
          }
    }

  H5Fclose(file_id);
  return (0);

exit5:
  H5Sclose(field_space);
exit4:
  H5Aclose(attr);
exit3:
  H5Gclose(field_set);
  H5Fclose(file_id);
  return (ecode);

exit2:
  H5Sclose(field_space);
exit1:
  H5Dclose(field_set);
exit0:
  H5Fclose(file_id);
  return (ecode);
}

// Find the good read invervals of the baxfile b(FileID), output the reads of length >= minLen and
//   score >= minScore to output (for the fasta or fastq part) and qvquiv (if b->quivqv is set)

#define HOLE   0
#define TYPE   1
#define    ADAPTER_REGION 0
#define    INSERT_REGION  1
#define    HQV_REGION     2
#define START  2
#define FINISH 3
#define SCORE  4

#ifdef OBSOLETE

static void writeBaxReads(BaxData *b, int minLen, int minScore, FILE *output, FILE* qvquiv)
{ int   nreads, *rlen;
  int   roff, *hlen, *cur, h, w;
  int   tolower;
  char *header;

  char   *baseCall;
  char   *delQV;
  char   *delTag;
  char   *insQV;
  char   *mergeQV;
  char   *subQV;
  char   *fastQV;

  baseCall = b->baseCall;
  delQV    = b->delQV;
  delTag   = b->delTag;
  insQV    = b->insQV;
  mergeQV  = b->mergeQV;
  subQV    = b->subQV;
  fastQV   = b->fastQV;

#ifdef DEBUG
  printf("printSubreadFields\n");
#endif

  //  Find the HQV regions and output as reads according to the various output options

  tolower = isupper(b->baseCall[0]);
  if (b->fastq)
    header = fastq_header;
  else
    header = fasta_header;

  rlen    = b->readLen;
  roff    = 0;
  cur     = b->regions;
  nreads  = b->numZMW + cur[HOLE];
  hlen    = rlen - cur[HOLE];
  cur[5*b->numHQR] = nreads; 

  for (h = cur[HOLE], w = 0; h < nreads; h++, w++)
    { int *bot, *top, *hqv, *r;
      int hbeg, hend, qv;
      int ibeg, iend;

      if (hlen[h] >= minLen)
        { while (cur[HOLE] < h)
            cur += 5;
          bot = hqv = cur;
          while (cur[HOLE] <= h)
            { if (cur[TYPE] == HQV_REGION)
                hqv = cur;
              cur += 5;
            }
          top = cur-5;

          qv = hqv[SCORE];
          if (qv >= minScore)
            { hbeg = hqv[START];
              hend = hqv[FINISH];
              for (r = bot; r <= top; r += 5)
                { if (r[TYPE] != INSERT_REGION)
                    continue;

                  ibeg = r[START];
                  iend = r[FINISH];

                  if (ibeg < hbeg)
                    ibeg = hbeg;
                  if (iend > hend)
                    iend = hend;
                  if (iend - ibeg < minLen || b->holeType[w] > 0)
                    continue;

                  fprintf(output,header,b->movieName,h,ibeg,iend,qv);

                  ibeg += roff;
                  iend += roff;

                  if (tolower)
                    { int a;

                      for (a = ibeg; a < iend; a++)
                        baseCall[a] += LOWER_OFFSET;
                      if (b->quivqv)
                        for (a = ibeg; a < iend; a++)
                          delTag[a] += LOWER_OFFSET;
                    }

                  if (b->fastq)
                    { int a;

                      fprintf(output,"%.*s\n", iend-ibeg, baseCall + ibeg);
                      fprintf(output,"+\n");
                      for (a = ibeg; a < iend; a++)
                        fputc(fastQV[a]+PHRED_OFFSET,output);
                      fputc('\n',output);
                    }
                  else
                    { int a;

                      for (a = ibeg; a < iend; a += 80)
                        if (a+80 > iend)
                          fprintf(output,"%.*s\n", iend-a, baseCall + a);
                        else
                          fprintf(output,"%.80s\n", baseCall + a);
                    }

                  if (b->quivqv)
                    { int a, d;

                      fprintf(qvquiv,"@%s/%d/%d_%d RQ=0.%d\n",
                                     b->movieName,h,ibeg-roff,iend-roff,qv);

                      d = b->delLimit;
                      for (a = ibeg; a < iend; a++)
                        { if (delQV[a] == d)
                            delTag[a] = 'n';
                          delQV[a]   += PHRED_OFFSET;
                          insQV[a]   += PHRED_OFFSET;
                          mergeQV[a] += PHRED_OFFSET;
                          subQV[a]   += PHRED_OFFSET;
                        }

                      iend -= ibeg;
                      fprintf (qvquiv, "%.*s\n", iend, delQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, delTag + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, insQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, mergeQV + ibeg);
                      fprintf (qvquiv, "%.*s\n", iend, subQV + ibeg);
                    }
                }
            }
        }
      roff += hlen[h];
    }
}

#endif

SubRead *nextSubread(BaxData *b, int prime)
{ static SubRead sub;

  static int roff, *cur, h, w;
  static int nreads, *hlen;

  static int *r, *top;
  static int  hbeg, hend;

#ifdef DEBUG
  printf("printSubreadFields\n");
#endif

  //  Find the HQV regions and output as reads according to the various output options

  if (prime)
    { cur     = b->regions;
      h       = cur[HOLE]-1;
      w       = -1;

      nreads  = b->numZMW;
      hlen    = b->readLen;
      roff    = 0;

      cur[5*b->numHQR] = nreads + cur[HOLE]; 

      r   = cur;
      top = cur;
      return (NULL);
    }

  for (r += 5 ; r <= top; r += 5)
    { int ibeg, iend;

      if (r[TYPE] != INSERT_REGION)
        continue;

      ibeg = r[START];
      iend = r[FINISH];

      if (ibeg < hbeg)
        ibeg = hbeg;
      if (iend > hend)
        iend = hend;
      if (iend-ibeg <= 0 || b->holeType[w] > 0)
        continue;

      sub.fpulse = ibeg;
      sub.lpulse = iend;
      return (&sub);
    }

  if (w >= 0)
    roff += hlen[w];
  for (h++, w++; w < nreads; h++, w++)
    { int *bot, *hqv, qv;
      int ibeg, iend;

      while (cur[HOLE] < h)
        cur += 5;
      bot = hqv = cur;
      while (cur[HOLE] <= h)
        { if (cur[TYPE] == HQV_REGION)
            hqv = cur;
          cur += 5;
        }
      top = cur-5;

      qv = hqv[SCORE];
      if (qv > 0)
        { hbeg = hqv[START];
          hend = hqv[FINISH];
          for (r = bot; r <= top; r += 5)
            { if (r[TYPE] != INSERT_REGION)
                continue;

              ibeg = r[START];
              iend = r[FINISH];

              if (ibeg < hbeg)
                ibeg = hbeg;
              if (iend > hend)
                iend = hend;
	      if (iend-ibeg <= 0 || b->holeType[w] > 0)
                continue;

              sub.well     = h;
              sub.fpulse   = ibeg;
              sub.lpulse   = iend;
              sub.qv       = qv;
              sub.data_off = roff;
              sub.zmw_off  = w;
              return (&sub);
            }
        }
      roff += hlen[w];
    }

  return (NULL);
}

//  Print an error message

void printBaxError(int ecode)
{ fprintf(stderr,"   ");
  switch (ecode)
    { case CANNOT_OPEN_BAX_FILE:
        fprintf(stderr,"Cannot open bax file\n");
        break;
      case BAX_BASECALL_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/Basecall\n");
        break;
      case BAX_DEL_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionQV\n");
        break;
      case BAX_TAG_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/DeletionTag\n");
        break;
      case BAX_INS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/InsertionQV\n");
        break;
      case BAX_MRG_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/MergeQV\n");
        break;
      case BAX_SUB_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/SubstitutionQV\n");
        break;
      case BAX_QV_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/QualityValue\n");
        break;
      case BAX_PULSE_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/WidthInFrames\n");
        break;
      case BAX_NR_EVENTS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/ZMW/NumEvent\n");
        break;
      case BAX_REGION_ERR:
        fprintf(stderr,"Cannot parse /PulseData/Regions\n");
        break;
      case BAX_HOLESTATUS_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/ZMW/HoleStatus\n");
        break;
      case BAX_MOVIENAME_ERR:
        fprintf(stderr,"Cannot parse /ScanData/RunInfo\n");
        break;
      case BAX_SNR_ERR:
        fprintf(stderr,"Cannot parse /PulseData/BaseCalls/ZMWMetrics/HQRegionSNR\n");
        break;
      case BAX_CHIPSET_ERR:
        fprintf(stderr,"Cannot parse /ScanData/DyeSet\n");
        break;
      default: 
        fprintf(stderr,"Cannot parse bax file\n");
        break;
    }
  fflush(stderr);
}

//  Free *the* bax data structure

void freeBaxData(BaxData *b)
{ free(b->baseCall);
  free(b->delQV);
  free(b->fastQV);
  free(b->holeType);
  free(b->readLen);
  free(b->regions);
  free(b->pulseW);
  free(b->snrVec);
}
