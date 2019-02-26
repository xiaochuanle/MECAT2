/*******************************************************************************************
 *
 *  Compressor for .quiv files, customized Huffman codes for each stream based on the
 *    histogram of values occuring in the given file.  The two low complexity streams
 *    (deletionQV and substitutionQV) use a Huffman coding of the run length of the prevelant
 *    character.
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
#include <stdint.h>

#include "DB.h"

static char *Usage = "[-vkl] <path:quiva> ...";

int main(int argc, char* argv[])
{ int        VERBOSE;
  int        KEEP;
  int        LOSSY;

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("dexqv")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        { ARG_FLAGS("vkl") }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    KEEP    = flags['k'];
    LOSSY   = flags['l'];

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -k: do *not* remove the .quiva file on completion.\n");
        fprintf(stderr,"      -l: use lossy compression (not recommended).\n");
        exit (1);
      }
  }

  // For each .quiva file to be compressed:

  { int       i;

    for (i = 1; i < argc; i++)
      { char     *pwd, *root;
        uint16    half;
        FILE     *input, *output;
        QVcoding *coding;
        
        pwd   = PathTo(argv[i]);
        root  = Root(argv[i],".quiva");
        input = Fopen(Catenate(pwd,"/",root,".quiva"),"r");
        if (input == NULL)
          exit (1);
        output = Fopen(Catenate(pwd,"/",root,".dexqv"),"w");
        if (output == NULL)
          exit (1);

        if (VERBOSE)
          { fprintf(stderr,"Processing '%s' ...\n",root);
            fflush(stderr);
          }

        //  Scan the file collecting statistics for Huffman schemes
        
        Set_QV_Line(0);
        QVcoding_Scan(input,INT32_MAX,NULL);

        //  Create and output the encoding schemes

        coding = Create_QVcoding(LOSSY);

        { char *slash, *read;   //  Get header line prefix from first line

          rewind (input);
          Read_Lines(input,1);
          read = QVentry();

          slash = index(read+1,'/');
          coding->prefix = (char *) malloc((slash-read)+1);
          if (coding->prefix == NULL)
             { fprintf(stderr,"%s: Out of memory (Allocating header prefix)\n",Prog_Name);
               exit (1);
             }
          *slash = '\0';
          strcpy(coding->prefix,read);
          *slash = '/';
        }

        half = 0x55aa;
        fwrite(&half,sizeof(uint16),1,output);

        Write_QVcoding(output,coding);

        //  For each entry do

        { int lwell;

          rewind (input);
          Set_QV_Line(0);

          lwell = 0;
          while (Read_Lines(input,1) > 0)
            { int    well, beg, end, qv;
              char  *slash;
              uint8  byte;

              //  Interpret the header, encode and write out the fields

              slash = index(QVentry(),'/');
              sscanf(slash+1,"%d/%d_%d RQ=0.%d\n",&well,&beg,&end,&qv);

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

              Compress_Next_QVentry(input,output,coding,LOSSY);
            }
        }

        //  Clean up for the next file

        Free_QVcoding(coding);

        fclose(input);
        fclose(output);

        if (!KEEP)
          unlink(Catenate(pwd,"/",root,".quiva"));
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
