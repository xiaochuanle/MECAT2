/*******************************************************************************************
 *
 *  Uncompresses a .dexta file (2-bit per base compression) back to a .fasta file
 *
 *  Author:  Gene Myers
 *  Date  :  January 12, 2014
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>

#include "DB.h"

static char *Usage = "[-vkU] [-w<int(80)>] ( -i | <path:dexta> ... )";

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
  int     UPPER;
  int     WIDTH;
  int     PIPE;

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("undexta")

    WIDTH   = 80;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vkiU")
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
    UPPER   = flags['U'];
    PIPE    = flags['i'];

    if ((PIPE && argc > 1) || (!PIPE && argc <= 1))
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -i: source is on standard input.\n");
        fprintf(stderr,"      -k: do *not* remove the .dexta file on completion.\n");
        fprintf(stderr,"      -U: use uppercase letters (default is lower case).\n");
        fprintf(stderr,"      -w: line width for sequence lines.\n");
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
            root  = Root(argv[i],".dexta");
            input = Fopen(Catenate(pwd,"/",root,".dexta"),"r");
            if (input == NULL)
              exit (1);
            output = Fopen(Catenate(pwd,"/",root,".fasta"),"w");
            if (output == NULL)
              exit (1);
          }

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        { char *name;
          int   well, flip, newv;

          // Read endian key and short name common to all headers

          { uint16 half;

            if (fread(&half,sizeof(uint16),1,input) != 1)
              SYSTEM_READ_ERROR
            if (half == 0x33cc)
              { flip = 0;
                newv = 0;
              }
            else if (half == 0xcc33)
              { flip = 1;
                newv = 0;
              }
            else if (half == 0x55aa)
              { flip = 0;
                newv = 1;
              }
            else if (half == 0xaa55)
              { flip = 1;
                newv = 1;
              }
            else
              { fprintf(stderr,"%s: Not a .dexta file, endian key invalid\n",Prog_Name);
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
            { int    rlen, beg, end, qv;
              int    clen;
              uint8  byte;

              //  Read and decompress header and output

              if (fread(&byte,1,1,input) < 1)
                break;
              while (byte == 255)
                { well += 255;
                  if (fread(&byte,1,1,input) != 1)
                    SYSTEM_READ_ERROR
                }
              well += byte;

              if (newv)
                if (flip)
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
                if (flip)
                  { uint16 half;

                    if (fread(&half,sizeof(uint16),1,input) != 1)
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
                  { uint16 half;

                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    beg = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    end = half;
                    if (fread(&half,sizeof(uint16),1,input) != 1)
                      SYSTEM_READ_ERROR
                    qv = half;
                  }

              fprintf(output,"%s/%d/%d_%d RQ=0.%d\n",name,well,beg,end,qv);

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
              if (UPPER)
                Upper_Read(read);
              else
                Lower_Read(read);

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
          unlink(Catenate(pwd,"/",root,".dexta"));
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
