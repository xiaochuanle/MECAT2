/*******************************************************************************************
 *
 *  Compresses a .fasta file into a 2-bit per base .dexta file
 *
 *  Author:  Gene Myers
 *  Date  :  January 12, 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vk] ( -i | <path:fasta> ... )";

#define MAX_BUFFER 100000

//  Compress read into 2-bits per base (from [0-3] per byte representation

int main(int argc, char *argv[])
{ int     VERBOSE;
  int     KEEP;
  int     PIPE;

  { int i, j, k;
    int flags[128];

    ARG_INIT("dexta")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vki") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];
    PIPE    = flags['i'];

    if ((PIPE && argc > 1) || (!PIPE && argc <= 1))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -i: source is on standard input.\n");
        fprintf(stderr,"      -k: do *not* remove the .fasta file on completion.\n");
        fprintf(stderr,"      -w: line width for sequence lines.\n");
        exit (1);
      }
    if (PIPE)
      { KEEP = 1;
        argc = 2;
      }
  }

  // For each fasta file do:

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Malloc(rmax+1,"Allocating read buffer");
    if (read == NULL)
      exit (1);

    for (i = 1; i < argc; i++)

      { char *pwd, *root;
        FILE *input, *output;
        int   eof;

        // Open fasta file

        if (PIPE)
          { input  = stdin;
            output = stdout;
            pwd    = NULL;
            root   = Strdup("Standard Input","Allocaing string");
          }
        else
          { pwd   = PathTo(argv[i]);
            root  = Root(argv[i],".fasta");
            input = Fopen(Catenate(pwd,"/",root,".fasta"),"r");
            if (input == NULL)
              exit (1);
            output = Fopen(Catenate(pwd,"/",root,".dexta"),"w");
            if (output == NULL)
              exit (1);
          }

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        // Read the first header and output the endian key and short name

        { char  *slash;
          uint16 half;
          int    x;

          eof = (fgets(read,MAX_BUFFER,input) == NULL);
          if (read[strlen(read)-1] != '\n')
            { fprintf(stderr,"Line 1: Fasta line is too long (> %d chars)\n",MAX_BUFFER-2);
              exit (1);
            }
          if (!eof && read[0] != '>')
            { fprintf(stderr,"Line 1: First header in fasta file is missing\n");
              exit (1);
            }

          slash = index(read,'/');
          if (slash == NULL)
            { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
              exit (1);
            }

          half = 0x55aa;
          fwrite(&half,sizeof(uint16),1,output);

          x = slash-read;
          fwrite(&x,sizeof(int),1,output);
          fwrite(read,1,slash-read,output);
        }

        //  For each read do

        { int  nline, rlen, lwell;

          nline = 1;
          rlen  = 0;
          lwell = 0;
          while (!eof)
            { int    well, beg, end, qv, x;
              char  *slash;
              uint8  byte;

              //  Next header is always at read+(rlen+1).  Interpret its fields

              slash = index(read+(rlen+1),'/');
              if (slash == NULL)
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
    		    exit (1);
                }
              x = sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv);
              if (x < 3)
                { fprintf(stderr,"%s: Header line incorrectly formatted ?\n",Prog_Name);
                  exit (1);
                }
              else if (x == 3)
                qv = 0;

              //  Read fasta sequence (@read) and stop at eof or after having read next header

              rlen = 0;
              while (1)
                { int x;

                  eof = (fgets(read+rlen,MAX_BUFFER,input) == NULL);
                  nline += 1;
                  x      = strlen(read+rlen)-1;
                  if (read[rlen+x] != '\n')
                    { fprintf(stderr,"Line %d: Fasta line is too long (> %d chars)\n",
                                     nline,MAX_BUFFER-2);
                      exit (1);
                    }
                  if (eof || read[rlen] == '>')
                    break;
                  rlen += x;
                  if (rlen + MAX_BUFFER > rmax)
                    { rmax = ((int) (1.2 * rlen)) + 1000 + MAX_BUFFER;
                      read = (char *) Realloc(read,rmax+1,"Reallocaing read buffer");
                      if (read == NULL)
                        exit (1);
                    }
                }
              read[rlen] = '\0';

              //  Compress the header fields and output (except for short name, only output once)

              while (well - lwell >= 255)
                { byte = 0xff;
                  fwrite(&byte,1,1,output);
                  lwell += 255;
                }
              byte = (uint8) (well-lwell);
              fwrite(&byte,1,1,output);
              lwell = well;

              fwrite(&beg,sizeof(int),1,output);
              fwrite(&end,sizeof(int),1,output);
              fwrite(&qv,sizeof(int),1,output);

              //  Compress read and output

              Number_Read(read);
              Compress_Read(rlen,read);
              fwrite(read,1,COMPRESSED_LEN(rlen),output);
            }
        }

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".fasta"));
        free(root);
        free(pwd);

        if (VERBOSE)
          { fprintf(stderr,"Done\n");
            fflush(stderr);
          }
      }

    free(read);
  }

  exit (0);
}
