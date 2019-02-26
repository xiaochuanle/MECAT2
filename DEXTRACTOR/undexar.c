/*******************************************************************************************
 *
 *  Uncompresses a .dexar file (2-bit per pulse width compression) back to an .arrow file
 *
 *  Author:  Gene Myers
 *  Date  :  October 12, 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vk] [-w<int(80)>] ( -i | <path:dexar> ... )";

#define MAX_BUFFER 100000

//  Uncompress read from 2-bits per base into [0-3] per byte representation

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

static void flip_short(void *w)
{ uint8 *v = (uint8 *) w;
  uint8  x;

  x    = v[0];
  v[0] = v[1];
  v[1] = x;
}

int main(int argc, char *argv[])
{ int     VERBOSE;
  int     KEEP;
  int     WIDTH;
  int     PIPE;

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("undexar")

    WIDTH   = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vki")
            break;
          case 'w':
            ARG_NON_NEGATIVE(WIDTH,"Line width")
            break;
        }
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
        fprintf(stderr,"      -k: do *not* remove the .dexar file on completion.\n");
        fprintf(stderr,"      -w: line width for arrow lines.\n");
        exit (1);
      }
    if (PIPE)
      { KEEP = 1;
        argc = 2;
      }
  }

  // For each .dexta file do

  { char   *read;
    int     rmax;
    int     i;

    rmax  = MAX_BUFFER + 30000;
    read  = (char *) Malloc(rmax+1,"Allocating read buffer");
    for (i = 1; i < argc; i++)
      { char *pwd, *root;
        FILE *input, *output;

        // Open dexta file

        if (PIPE)
          { input  = stdin;
            output = stdout;
            pwd    = NULL;
            root   = Strdup("Standard Input","Allocaing string");
          }
        else
          { pwd   = PathTo(argv[i]);
            root  = Root(argv[i],".dexar");
            input = Fopen(Catenate(pwd,"/",root,".dexar"),"r");
            if (input == NULL)
              exit (1);
            output = Fopen(Catenate(pwd,"/",root,".arrow"),"w");
            if (output == NULL)
              exit (1);
          }

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        { char *name;
          int   well, flip;

          // Read endian key and short name common to all headers

          { uint16 half;

            if (fread(&half,sizeof(uint16),1,input) != 1)
              SYSTEM_READ_ERROR
            if (half == 0x55aa)
              flip = 0;
            else if (half == 0xaa55)
              flip = 1;
            else
              { fprintf(stderr,"%s: Not a .dexar file, endian key invalid\n",Prog_Name);
                exit (1);
              }

            if (fread(&well,sizeof(int),1,input) != 1)
              SYSTEM_READ_ERROR
            if (flip) flip_long(&well);
            name = (char *) Malloc(well+1,"Allocating header prefix");
            if (well > 0)
              { if (fread(name,well,1,input) != 1)
                  SYSTEM_READ_ERROR
              }
            name[well] = '\0';
          }

          // For each encoded entry do

          well = 0;
          while (1)
            { int    rlen, beg, end, x;
              float  snr[4];
              uint16 cnr[4];
              int    clen;
              uint8  byte;

              //  Read and decompress header and output

              if (fread(&byte,1,1,input) < 1) break;
              while (byte == 255)
                { well += 255;
                  if (fread(&byte,1,1,input) != 1)
                    SYSTEM_READ_ERROR
                }
              well += byte;

              if (flip)
                { if (fread(&beg,sizeof(int),1,input) != 1)
                    SYSTEM_READ_ERROR
                  flip_long(&beg);
                  if (fread(&end,sizeof(int),1,input) != 1)
                    SYSTEM_READ_ERROR
                  flip_long(&end);
                  if (fread(cnr,sizeof(uint16),4,input) != 4)
                    SYSTEM_READ_ERROR
                  for (x = 0; x < 4; x++)
                    flip_short(cnr+x);
                }
              else
                { if (fread(&beg,sizeof(int),1,input) != 1)
                    SYSTEM_READ_ERROR
                  if (fread(&end,sizeof(int),1,input) != 1)
                    SYSTEM_READ_ERROR
                  if (fread(cnr,sizeof(uint16),4,input) != 4)
                    SYSTEM_READ_ERROR
                }

              for (x = 0; x < 4; x++)
                snr[x] = cnr[x]/100.;

              fprintf(output,"%s/%d/%d_%d SN=%.2f,%.2f,%.2f,%.2f\n",name,well,beg,end,
                                                                    snr[0],snr[1],snr[2],snr[3]);

              //  Read compressed sequence (into buffer big enough for uncompressed sequence)
              //  Uncompress and output WIDTH symbols to a line

              rlen = end-beg;
              if (rlen > rmax)
                { rmax = ((int) (1.2 * rlen)) + 1000 + MAX_BUFFER;
                  read = (char *) Realloc(read,rmax+1,"Allocating read buffer");
                }
              clen = COMPRESSED_LEN(rlen);
              if (clen > 0)
                { if (fread(read,clen,1,input) != 1)
                    SYSTEM_READ_ERROR
                }
              Uncompress_Read(rlen,read);
              Letter_Arrow(read);

              { int j;

                for (j = 0; j < rlen; j += WIDTH)
                  if (j+WIDTH > rlen)
                    fprintf(output,"%.*s\n", rlen-j, read+j);
                  else
                    fprintf(output,"%.*s\n", WIDTH, read+j);
              }
            }

          free(name);
        }

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".dexar"));
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
