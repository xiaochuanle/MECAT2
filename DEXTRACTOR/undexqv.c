/*******************************************************************************************
 *
 *  Uncompressor for .dexqv files
 *
 *  Author:  Gene Myers
 *  Date:    Jan 18, 2014
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "DB.h"

static char *Usage = "[-vkU] <path:dexqv> ...";

static void flip_short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[1];
  v[1] = x;
}

static void flip_long(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[3];
  v[3] = x;
  x    = v[1];
  v[1] = v[2];
  v[2] = x;
}

int main(int argc, char* argv[])
{ int VERBOSE;
  int KEEP;
  int UPPER;

  { int i, j, k;
    int flags[128];

    ARG_INIT("undexqv")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vkU") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];
    UPPER   = flags['U'];

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: do *not* remove the .dexqv file on completion.\n");
        fprintf(stderr,"      -U: use uppercase letters (default is lower case).\n");
        exit (1);
      }
  }

  //  For each .dexqv file to be decompressed

  { int   i;
    char *entry[5] = { NULL, NULL, NULL, NULL, NULL };
    int   emax     = -1;
    
    for (i = 1; i < argc; i++)
      { char     *pwd, *root;
        uint16    half;
        int       newv;
        FILE     *input, *output;
        QVcoding *coding;

        //   Open it and the appropriately named .quiva file

        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".dexqv");
        input = Fopen(Catenate(pwd,"/",root,".dexqv"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".quiva"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        // Read in compression scheme

        if (fread(&half,sizeof(uint16),1,input) != 1)
          SYSTEM_READ_ERROR
        if (half == 0x55aa || half == 0xaa55)
          newv = 1;
        else
          { newv = 0;
            rewind(input);
          }

        coding = Read_QVcoding(input);

        //  For each compressed entry do

        { int well;

          well = 0;
          while (1)
            { int    beg, end, qv, rlen;
              uint16 half;
              uint8  byte;
              int    e;

              //  Decode the compressed header and write it out

              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  if (fread(&byte,1,1,input) != 1)
                    SYSTEM_READ_ERROR
                }
              well += byte;

              if (newv)
                if (coding->flip)
                  { if (fread(&beg,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_long(&beg);
                    if (fread(&end,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_long(&end);
                    if (fread(&qv,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_long(&qv);
                  }
                else
                  { if (fread(&beg,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                    if (fread(&end,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                    if (fread(&qv,sizeof(int),1,input) != 1)
                      SYSTEM_READ_ERROR
                  }
              else
                if (coding->flip)
                  { if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_short(&half);
                    beg = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_short(&half);
                    end = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    flip_short(&half);
                    qv = half;
                  }
                else
                  { if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    beg = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    end = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    qv = half;
                  }

              fprintf(output,"%s/%d/%d_%d RQ=0.%d\n",coding->prefix,well,beg,end,qv);

              //  Decode the QV entry and write it out

              rlen = end-beg;
              if (rlen > emax)
                { emax = ((int) (1.2*rlen)) + 1000;
                  entry[0] = (char *) Realloc(entry[0],5*emax,"Reallocating QV entry buffer");
                  if (entry[0] == NULL)
                    exit (1);
                  for (e = 1; e < 5; e++)
                    entry[e] = entry[e-1] + emax;
                }

              Uncompress_Next_QVentry(input,entry,coding,rlen);

              if (UPPER)
                { char *deltag = entry[1];
                  int   j;

                  for (j = 0; j < rlen; j++)
                    deltag[j] -= 32;
                }

              for (e = 0; e < 5; e++)
                fprintf(output,"%.*s\n",rlen,entry[e]);
            }
	}

        //  Clean up for the next file

	Free_QVcoding(coding);

        fclose(input);
        fclose(output);

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".dexqv"));
        free(root);
        free(pwd);

        if (VERBOSE)
          { fprintf(stderr,"Done\n");
            fflush(stderr);
          }
      }
  }

  free(QVentry());

  exit (0);
}
