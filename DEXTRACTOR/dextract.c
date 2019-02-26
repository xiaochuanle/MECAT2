/*******************************************************************************************
 *
 *  Dextract: pullls requested info out of subreads.[bs]am and .bax.h5 files produced by
 *               Pacbio sequencing instruments and software
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 9, 2016
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "DB.h"
#include "sam.h"
#include "bax.h"
#include "expr.h"

#define LOWER_OFFSET 32
#define PHRED_OFFSET 33

static char *Usage = "[-vfaq] [-o[<path>]] [-e<expr(ln>=500 && rq>=750)> <input:pacbio> ...";

  //  Write subreads s from bax data set b to non-NULL file types

static void writeSubread(BaxData *b, SubRead *s, FILE *fas, FILE *arr, FILE* qvs)
{ int ibeg, iend, roff, len;

  ibeg = s->fpulse;
  iend = s->lpulse;
  roff = s->data_off + ibeg;
  len  = iend - ibeg;

  if (arr != NULL)  //  .arrow
    { int     a;
      uint16 *pulse;
      float  *snr;

      pulse = b->pulseW + roff;
      snr   = b->snrVec + 4*s->zmw_off;

      fprintf(arr,">%s SN=%.2f",b->movieName,snr[b->chan[0]]);
      for (a = 1; a < 4; a++)
        fprintf(arr,",%.2f",snr[b->chan[a]]);
      fprintf(arr,"\n");

      for (a = 0; a < len; a++)
        { if (pulse[a] >= 4)
            fputc('4',arr);
          else
            fputc(pulse[a]+'0',arr);
          if (a % 80 == 79)
            fputc('\n',arr);
        }
      if (a % 80 != 80)
        fputc('\n',arr);
    }

  if (fas != NULL)   //   .fasta
    { int   a;
      char *baseCall;

      baseCall = b->baseCall + roff;

      fprintf(fas,">%s/%d/%d_%d RQ=0.%0d\n",b->movieName,s->well,ibeg,iend,s->qv);

      if (isupper(baseCall[ibeg]))
        for (a = 0; a < len; a++)
          baseCall[a] += LOWER_OFFSET;

      for (a = 0; a < len; a += 80)
        if (a+80 > len)
          fprintf(fas,"%.*s\n", len-a, baseCall + a);
        else
          fprintf(fas,"%.80s\n", baseCall + a);
    }

  if (qvs != NULL)    //   .quiva
    { int   a, d;
      char *delQV, *delTag, *insQV, *mergeQV, *subQV;

      delQV   = b->delQV + roff;
      delTag  = b->delTag + roff;
      insQV   = b->insQV + roff;
      mergeQV = b->mergeQV + roff;
      subQV   = b->subQV + roff;

      fprintf(qvs,"@%s/%d/%d_%d RQ=0.%0d\n",b->movieName,s->well,ibeg,iend,s->qv);

      if (isupper(delTag[ibeg]))
        for (a = 0; a < len; a++)
          delTag[a] += LOWER_OFFSET;
      d = b->delLimit;
      if (isupper(d))
        d += LOWER_OFFSET;

      for (a = 0; a < len; a++)
        { if (delQV[a] == d)
            delTag[a] = 'n';
          if (delQV[a] > 93)
            delQV[a] = 126;
          else
            delQV[a] += PHRED_OFFSET;
          if (insQV[a] > 93)
            insQV[a] = 126;
          else
            insQV[a] += PHRED_OFFSET;
          if (mergeQV[a] > 93)
            mergeQV[a] = 126;
          else
            mergeQV[a] += PHRED_OFFSET;
          if (subQV[a] > 93)
            subQV[a] = 126;
          else
            subQV[a] += PHRED_OFFSET;
        }

      fprintf (qvs, "%.*s\n", len, delQV);
      fprintf (qvs, "%.*s\n", len, delTag);
      fprintf (qvs, "%.*s\n", len, insQV);
      fprintf (qvs, "%.*s\n", len, mergeQV);
      fprintf (qvs, "%.*s\n", len, subQV);
    }
}

  //  Write subread data in samRecord rec to non-NULL file types

static void writeSamRecord(samRecord *rec, FILE *fas, FILE *arr, FILE* qvs)
{ int i;

  if (fas != NULL)
    { fprintf(fas,">%s/%d/%d_%d RQ=0.%d\n",rec->header,rec->well,rec->beg,
                                               rec->end,(int) (rec->qual*1000.));
      for (i = 0; i < rec->len; i += 80)
        if (i+80 <= rec->len)
          fprintf(fas,"%.80s\n",rec->seq+i);
        else
          fprintf(fas,"%.*s\n",rec->len-i,rec->seq+i);
    }

  if (arr != NULL)
    { fprintf(arr,">%s SN=%.2f,%.2f,%.2f,%.2f\n",rec->header,
                      rec->snr[0],rec->snr[1],rec->snr[2],rec->snr[3]);
      for (i = 0; i < rec->len; i += 80)
        if (i+80 <= rec->len)
          fprintf(arr,"%.80s\n",rec->arr+i);
        else
          fprintf(arr,"%.*s\n",rec->len-i,rec->arr+i);
    }

  if (qvs != NULL)
    { fprintf(qvs,">%s/%d/%d_%d RQ=0.%d\n",rec->header,rec->well,rec->beg,
                                               rec->end,(int) (rec->qual*1000.));
      fprintf(qvs,"%.*s\n",rec->len,rec->qv[0]);
      fprintf(qvs,"%.*s\n",rec->len,rec->qv[1]);
      fprintf(qvs,"%.*s\n",rec->len,rec->qv[2]);
      fprintf(qvs,"%.*s\n",rec->len,rec->qv[3]);
      fprintf(qvs,"%.*s\n",rec->len,rec->qv[4]);
    }
}

  //  Main

int main(int argc, char* argv[])
{ char *output;
  char *path, *core;
  FILE *fileFas;
  FILE *fileArr;
  FILE *fileQvs;

  int     ARROW;
  int     QUIVA;
  int     FASTA;
  int     VERBOSE;
  Filter *EXPR;

  //  Process command line arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("dextract")

    path   = NULL;
    core   = NULL;
    output = NULL;
    EXPR   = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vfaq")
            break;
          case 'o':
            output = argv[i]+2;
            break;
          case 'e':
            EXPR = parse_filter(argv[i]+2);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    ARROW   = flags['a'];
    QUIVA   = flags['q'];
    FASTA   = flags['f'];
    if ( ! (ARROW || FASTA || QUIVA))
      FASTA = 1;

    if (EXPR == NULL)
      EXPR = parse_filter("ln>=500 && rq>=750");

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: extract a .fasta file with Pacbio-style line headers.\n");
        fprintf(stderr,"      -a: extract a .arrow file with SNR encoded in line headers.\n");
        fprintf(stderr,"      -q: extract a .quiva file with Pacbio-style line headers.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -o: If absent, output files use root name of input .bax or .bam.\n");
        fprintf(stderr,"        : If no path given, output sent to standard output.\n");
        fprintf(stderr,"        : If path given, output files use path name as root name.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: subread selection expression.  Possible variables are:\n");
        fprintf(stderr,"           zm  - well number\n");
        fprintf(stderr,"           ln  - length of subread\n");
        fprintf(stderr,"           rq  - quality value of subread (normalized to [0,1000])\n");
        fprintf(stderr,"           bc1 - # of first barcode\n");
        fprintf(stderr,"           bc2 - # of second barcode\n");
        fprintf(stderr,"           bq  - quality of barcode detection (normalized to [0,100])\n");
        fprintf(stderr,"           np  - number of passes producing subread\n");
        fprintf(stderr,"           qs  - start pulse of subread\n");
        exit (1);
      }
  }

  //  If -o set then set up output file streams

  fileFas = NULL;
  fileArr = NULL;
  fileQvs = NULL;
  if (output != NULL)
    { if (*output != '\0')
        { path   = PathTo(output);
          output = Root(output,NULL);

          if (FASTA)
            { fileFas = Fopen(Catenate(path,"/",output,".fasta"), "w");
              if (fileFas == NULL)
                goto error;
            }
          if (ARROW)
            { fileArr = Fopen(Catenate(path,"/",output,".arrow"), "w");
              if (fileArr == NULL)
                goto error;
            }
          if (QUIVA)
            { fileQvs = Fopen(Catenate(path,"/",output,".quiva"), "w");
              if (fileQvs == NULL)
                goto error;
            }
          free(path);
        }
      else
        { if (ARROW + FASTA + QUIVA > 1)
            { fprintf(stderr,"%s: Cannot send more than one stream to standard output\n",
                             Prog_Name);
              exit (1);
            }
          if (FASTA)
            fileFas = stdout;
          if (ARROW)
            fileArr = stdout;
          if (QUIVA)
            fileQvs = stdout;
        }
    }
 
  //  Process each input file

  { int      i;
    BaxData  b, *bp = &b;
    samFile *in;

    initBaxData(bp,0,QUIVA,ARROW);

    for (i = 1; i < argc; i++)
      { FILE *file;
        int   status, intype;

        //  Determine file type

#define IS_BAX 0
#define IS_BAM 1
#define IS_SAM 2

        path  = PathTo(argv[i]);
        core  = Root(argv[i],".subreads.bam");
        if ((file = fopen(Catenate(path,"/",core,".subreads.bam"),"r")) == NULL)
          { core  = Root(argv[i],".subreads.sam");
            if ((file = fopen(Catenate(path,"/",core,".subreads.sam"),"r")) == NULL)
              { core  = Root(argv[i],".bax.h5");
                if ((file = fopen(Catenate(path,"/",core,".bax.h5"),"r")) == NULL)
                  { fprintf(stderr,"%s: Cannot find %s/%s with a Pacbio extension\n",
                                   Prog_Name,path,core);
                    goto error;
                  }
                intype = IS_BAX;
              }
            else
              intype = IS_SAM;
          }
        else
          intype = IS_BAM;
        fclose(file);

        //  If -o not set then setup output file streams for this input

        if (output == NULL)
          { if (FASTA)
              { fileFas = Fopen(Catenate(path,"/",core,".fasta"), "w");
                if (fileFas == NULL)
                  goto error;
              }
            if (ARROW)
              { fileArr = Fopen(Catenate(path,"/",core,".arrow"), "w");
                if (fileArr == NULL)
                  goto error;
              }
            if (QUIVA)
              { fileQvs = Fopen(Catenate(path,"/",core,".quiva"), "w");
                if (fileQvs == NULL)
                  goto error;
              }
          }

        //  Extract from a .bax.h5

        if (intype == IS_BAX)
          { SubRead *s;

            if (VERBOSE)
              { fprintf(stderr, "Fetching file : %s ...\n", core); fflush(stderr); }

            if ((status = getBaxData(bp,Catenate(path,"/",core,".bax.h5"))) != 0)
              { fprintf(stderr, "%s: ", Prog_Name);
                printBaxError(status);
                goto error;
              }

            if (VERBOSE)
              { fprintf(stderr, "Extracting subreads ...\n"); fflush(stderr); }

            nextSubread(bp,1);
            while (1)
              { s = nextSubread(bp,0);
                if (s == NULL)
                  break;

                if ( ! evaluate_bax_filter(EXPR,bp,s))
                  continue;

                writeSubread(&b,s,fileFas,fileArr,fileQvs);
              }
          }

        //  Extract from a .bam or .sam

        else
          { if (VERBOSE)
              { fprintf(stderr, "Processing file : %s ...\n", core); fflush(stderr); }

            if (intype == IS_BAM)
              { if ((in = sam_open(Catenate(path,"/",core,".subreads.bam"))) == NULL)
                  { fprintf(stderr, "%s: can't open %s as a Bam file\n", Prog_Name, argv[i]);
                    goto error;
                  }
              }
            else
              { if ((in = sam_open(Catenate(path,"/",core,".subreads.sam"))) == NULL)
                  { fprintf(stderr, "%s: can't open %s as a Sam file\n", Prog_Name, argv[i]);
                    goto error;
                  }
              }

            status = sam_header_process(in,0);
            if (status < 0)
              goto error;
            else if ((status & HASPW) == 0 && ARROW)
              { fprintf(stderr, "%s: %s does not have Arrow information\n", Prog_Name, argv[i]);
                goto error;
              }
            else if ((status & HASQV) == 0 && QUIVA)
              { fprintf(stderr, "%s: %s does not have Quiver information\n", Prog_Name, argv[i]);
                goto error;
              }
            else
              { samRecord *rec;
  
                while (1)
                  { rec = sam_record_extract(in, status);
                    if (rec == NULL)
                      goto error;
                    if (rec == SAM_EOF)
                      break;

                    if ( ! evaluate_bam_filter(EXPR,rec))
                      continue;

                    writeSamRecord(rec,fileFas,fileArr,fileQvs);
                  }
              }

            if (sam_close(in))
              { fprintf(stderr, "%s: Error closing file %s\n", Prog_Name, core);
                goto error;
              }
          }

        //  If -o not set, close outputs for input file and free name strings

        if (output == NULL)
          { if (FASTA)
              fclose(fileFas);
            if (ARROW)
              fclose(fileArr);
            if (QUIVA)
              fclose(fileQvs);
            fileFas = NULL;
            fileQvs = NULL;
            fileArr = NULL;
          }

        free(path);
        free(core);

        if (VERBOSE)
          { fprintf(stderr, "Done\n"); fflush(stdout); }
      }
  }

  //  If -o<name> then close named outputs

  if (output != NULL && *output != '\0')
    { if (fileFas != NULL)
        fclose(fileFas);
      if (fileArr != NULL)
        fclose(fileArr);
      if (fileQvs != NULL)
        fclose(fileQvs);
      free(output);
    }

  exit (0);

  //  An error occured, carefully undo any files in progress

error:
  if (output == NULL)
    { if (fileFas != NULL)
        { fclose(fileFas);
          unlink(Catenate(path,"/",core,".fasta"));
        }
      if (fileQvs != NULL)
        { fclose(fileQvs);
          unlink(Catenate(path,"/",core,".quiva"));
        }
      if (fileArr != NULL)
        { fclose(fileArr);
          unlink(Catenate(path,"/",core,".arrow"));
        }
    }
  else if (*output != '\0')
    { if (fileFas != NULL)
        { fclose(fileFas);
          unlink(Catenate("","",output,".fasta"));
        }
      if (fileQvs != NULL)
        { fclose(fileQvs);
          unlink(Catenate("","",output,".quiva"));
        }
      if (fileArr != NULL)
        { fclose(fileArr);
          unlink(Catenate("","",output,".arrow"));
        }
      free(output);
    }
  free(path);
  free(core);

  exit (1);
}
