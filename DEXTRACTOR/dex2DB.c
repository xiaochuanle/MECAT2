/*******************************************************************************************
 *
 *  Add PacBio .bax.h5 or .subreads.[bs]am files to a DB:
 *     Adds the given HDF5 or BAM files in the given order to <path>.db.  If the db does not exist
 *     then it is created.  All the HDF5 files added to a given data base must be Pacbio .bax.h5
 *     files, and all the BAM files added to a given data base must be Pacbio .subreads.bam or .sam
 *     files.  A file cannot be added twice and this is enforced.  Initially one must select
 *     to create either an Arrow (-A) or Quiver (-Q) database, and on subsequent adds the type
 *     flag if set must agree with the type of the database.  The command either builds
 *     or appends to the .<path>.idx, .<path>.bps, and <path>.(qvs|arw) files, where the index
 *     file (.idx) contains information about each read and the offsets to a 2-bit encoded
 *     sequence in the base-pair file (.bps) and also the arrow file (.arw) if it is an Arrow
 *     database.  If a Quiver database, then the .qvs file is the concatentation of the Huffman
 *     table and the compressed quality streams for each added file.  All the files are hidden
 *     by virtue of their names beginning with a '.'.  <path>.db is effectively a stub file
 *     with given name that contains an ASCII listing of the files added to the DB and possibly
 *     the block partitioning for the DB if DBsplit has been called upon it.
 *
 *  Author:  Gene Myers
 *  Date  :  October 2016
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>

#include "DB.h"
#include "sam.h"
#include "bax.h"
#include "expr.h"

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage[] =
         { "[-vlaq] [-e<expr(ln>=500 && rq>=750)>]",
           "  <path:string> ( -f<file> | <input:pacbio> ... )"
         };

typedef struct
  { int    argc;
    char **argv;
    FILE  *input;
    int    count;
    char  *name;
  } File_Iterator;

File_Iterator *init_file_iterator(int argc, char **argv, FILE *input, int first)
{ File_Iterator *it;

  it = Malloc(sizeof(File_Iterator),"Allocating file iterator");
  if (it == NULL)
    return (NULL);
  it->argc  = argc;
  it->argv  = argv;
  it->input = input;
  if (input == NULL)
    it->count = first;
  else
    { it->count = 1;
      rewind(input);
    }
  return (it);
}

int next_file(File_Iterator *it)
{ static char nbuffer[MAX_NAME+8];

  if (it->input == NULL)
    { if (it->count >= it->argc)
        return (0);
      it->name = it->argv[it->count++];
    }
  else
    { char *eol;

      if (fgets(nbuffer,MAX_NAME+8,it->input) == NULL)
        { if (feof(it->input))
            return (0);
          fprintf(stderr,"%s: IO error reading line %d of -f file of names\n",Prog_Name,it->count);
          it->name = NULL;
          return (1);
        }
      if ((eol = index(nbuffer,'\n')) == NULL)
        { fprintf(stderr,"%s: Line %d in file list is longer than %d chars!\n",
                         Prog_Name,it->count,MAX_NAME+7);
          it->name = NULL;
          return (1);
        }
      *eol = '\0';
      it->count += 1;
      it->name  = nbuffer;
    }
  return (1);
}

static char Number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };


int main(int argc, char *argv[])
{ FILE  *istub, *ostub;
  char  *dbname;
  char  *root, *pwd;

  FILE  *bases, *indx, *quiva, *arrow;
  int64  boff, ioff, coff;

  int    ifiles, ofiles, ocells;
  char **flist;

  DAZZ_DB db;
  int     ureads;
  int64   offset;

  FILE   *IFILE;
  int     VERBOSE;
  int     LOSSY;
  int     ARROW;
  int     QUIVER;
  Filter *EXPR;

  //   Process command line

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("dex2DB")

    IFILE = NULL;
    EXPR  = NULL;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vlaq")
            break;
          case 'f':
            IFILE = fopen(argv[i]+2,"r");
            if (IFILE == NULL)
              { fprintf(stderr,"%s: Cannot open file of inputs '%s'\n",Prog_Name,argv[i]+2);
                exit (1);
              }
            break;
          case 'e':
            EXPR = parse_filter(argv[i]+2);
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    LOSSY   = flags['l'];
    ARROW   = flags['a'];
    QUIVER  = flags['q'];

    if (EXPR == NULL)
      EXPR = parse_filter("ln>=500 && rq>=750");
     
    if ( (IFILE == NULL && argc <= 2) || (IFILE != NULL && argc != 2) )
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -f: build or add to DB the files listed in the -f file.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -a: Build or add to an arrow DB.\n");
        fprintf(stderr,"      -q: Build or add to a quiva DB.\n");
        fprintf(stderr,"      -l: Use lossy compression (with -q option only).\n");
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
    if (ARROW && QUIVER)
      { fprintf(stderr,"%s: Cannot set both -a(rrow) and -q(uiver)\n",Prog_Name);
        exit (1);
      }
  }

  //  Try to open DB file, if present then adding to DB, otherwise creating new DB.  Set up
  //  variables as follows:
  //    dbname = full name of db = <pwd>/<root>.db
  //    istub  = open db file (if adding) or NULL (if creating)
  //    ostub  = new image of db file (will overwrite old image at end)
  //    bases  = .bps file positioned for appending
  //    indx   = .idx file positioned for appending
  //    quiva  = .qvs file positioned for appending (if QUIVER)
  //    arrow  = .arw file positioned for appending (if ARROW)
  //    ureads = # of reads currently in db
  //    offset = offset in .bps at which to place next sequence/arrow
  //    ioff   = offset in .idx file to truncate to if command fails
  //    boff   = offset in .bps file to truncate to if command fails
  //    coff   = offset in .qvs file to truncate to if command fails (if QUIVER)
  //    ifiles = # of .bam files to add
  //    ofiles = # of .bam files added so far
  //    ocells = # of SMRT cells already in db
  //    flist  = [0..ifiles+ocells] list of file names (root only) added to db so far

  { int i;

    root   = Root(argv[1],".db");
    pwd    = PathTo(argv[1]);
    dbname = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db name");
    if (dbname == NULL)
      exit (1);

    if (IFILE == NULL)
      ifiles = argc-2;
    else
      { File_Iterator *ng;

        ifiles = 0;
        ng = init_file_iterator(argc,argv,IFILE,2);
        if (ng == NULL)
          exit (1);
        while (next_file(ng))
          { if (ng->name == NULL)
              exit (1);
            ifiles += 1;
          }
        free(ng);
      }

    bases = NULL;
    indx  = NULL;
    quiva = NULL;
    arrow = NULL;
    ostub = NULL;
    ioff  = 0;
    boff  = 0;
    coff  = 0;

    istub = fopen(dbname,"r");
    if (istub == NULL)
      { ocells = 0;

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"w+");
        indx  = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"w+");
        if (bases == NULL || indx == NULL)
          goto error;

        if (QUIVER)
          { quiva = Fopen(Catenate(pwd,PATHSEP,root,".qvs"),"w+");
            if (quiva == NULL)
              goto error;
          }
        if (ARROW)
          { arrow = Fopen(Catenate(pwd,PATHSEP,root,".arw"),"w+");
            if (arrow == NULL)
              goto error;
          }

        fwrite(&db,sizeof(DAZZ_DB),1,indx);

        ureads  = 0;
        offset  = 0;
      }
    else
      { DAZZ_READ rec;

        if (fscanf(istub,DB_NFILE,&ocells) != 1)
          { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
            exit (1);
          }

        indx = Fopen(Catenate(pwd,PATHSEP,root,".idx"),"r+");
        if (indx == NULL)
          exit (1);

        if (fread(&db,sizeof(DAZZ_DB),1,indx) != 1)
          { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
            exit (1);
          }

        fseeko(indx, -sizeof(DAZZ_READ), SEEK_END);
        fread(&rec,sizeof(DAZZ_READ),1,indx);
        if (rec.coff < 0)
          { if (ARROW || QUIVER)
              { fprintf(stderr,"%s: Sequence DB but you set either the -a or -q flag?\n",
                               Prog_Name);
                exit (1);
              }
          }
        else if ((db.allarr & DB_ARROW) != 0)
          { if (QUIVER)
              { fprintf(stderr,"%s: Arrow DB but you set the -q flag?\n",Prog_Name);
                exit (1);
              }
            ARROW = 1;
          }
        else
          { if (ARROW)
              { fprintf(stderr,"%s: Quiver DB but you set the -a flag?\n",Prog_Name);
                exit (1);
              }
            QUIVER = 1;
          }

        bases = Fopen(Catenate(pwd,PATHSEP,root,".bps"),"r+");
        if (bases == NULL)
          exit (1);
        if (QUIVER)
          { quiva = Fopen(Catenate(pwd,PATHSEP,root,".qvs"),"r+");
            if (quiva == NULL)
              exit (1);
            fseeko(quiva,0,SEEK_END);
            coff = ftello(quiva);
          }
        if (ARROW)
          { arrow = Fopen(Catenate(pwd,PATHSEP,root,".arw"),"r+");
            if (arrow == NULL)
              exit (1);
            fseeko(arrow,0,SEEK_END);
          }

        fseeko(bases,0,SEEK_END);
        fseeko(indx, 0,SEEK_END);

        ureads = db.ureads;
        offset = ftello(bases);
        boff   = offset;
        ioff   = ftello(indx);
      }

    if (!QUIVER && LOSSY)
      fprintf(stderr,"%s: Warning: Option -l set but not adding Quiver information?\n",Prog_Name);

    flist  = (char **) Malloc(sizeof(char *)*(ocells+ifiles),"Allocating file list");
    ostub  = Fopen(Catenate(pwd,"/",root,".dbx"),"w+");
    if (ostub == NULL || flist == NULL)
      goto error;

    fprintf(ostub,DB_NFILE,ocells+ifiles);   //  Will write again with correct value at end
    ofiles = 0;
    for (i = 0; i < ocells; i++)
      { int  last;
        char prolog[MAX_NAME], fname[MAX_NAME];

        if (fscanf(istub,DB_FDATA,&last,fname,prolog) != 3)
          { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
            goto error;
          }
        if (ofiles == 0 || strcmp(flist[ofiles-1],fname) != 0)
          if ((flist[ofiles++] = Strdup(fname,"Adding to file list")) == NULL)
            goto error;
        fprintf(ostub,DB_FDATA,last,fname,prolog);
      }
  }

  { int            maxlen;
    int64          totlen, count[4];
    int            pmax;
    DAZZ_READ     *prec;
    int            c;
    File_Iterator *ng = NULL;
    BaxData       _bax, *bax = &_bax;
    samFile       *input;

    //  Buffer for reads all in the same well

    pmax = 100;
    prec = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*pmax,"Allocating record buffer");
    if (prec == NULL)
      goto error;

    totlen = 0;              //  total # of bases in new .bam files
    maxlen = 0;              //  longest read in new .bam files
    for (c = 0; c < 4; c++)  //  count of acgt in new .bam files
      count[c] = 0;

    //  For each new input source do

    ng = init_file_iterator(argc,argv,IFILE,2);  //  Setup to read .bam's
    if (ng == NULL)                              //    from command line or file
      goto error;

    initBaxData(bax,0,QUIVER,ARROW);

    while (next_file(ng))
      { FILE    *file;
        char    *path, *core;
        int      status, empty, intype;

        if (ng->name == NULL) goto error;

        //  Determine file type

#define IS_BAX 0
#define IS_BAM 1
#define IS_SAM 2

        path  = PathTo(ng->name);
        core  = Root(ng->name,".subreads.bam");
        if ((file = fopen(Catenate(path,"/",core,".subreads.bam"),"r")) == NULL)
          { free(core);
            core  = Root(ng->name,".subreads.sam");
            if ((file = fopen(Catenate(path,"/",core,".subreads.sam"),"r")) == NULL)
              { free(core);
                core  = Root(ng->name,".bax.h5");
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
        empty = (fgetc(file) == EOF);
        fclose(file);

        if (empty)
          { free(path);
            free(core);
            fprintf(stderr,"Skipping '%s', file is empty!\n",core);
            continue;
          }

        //  Check that core is not too long and name is unique

        if (strlen(core) >= MAX_NAME)
          { fprintf(stderr,"%s: File name over %d chars: '%.200s'\n",
                           Prog_Name,MAX_NAME,core);
            goto error;
          }

        { int j;

          for (j = 0; j < ofiles; j++)
            if (strcmp(core,flist[j]) == 0)
              { fprintf(stderr,"%s: File %s is already in database %s.db\n",
                               Prog_Name,core,Root(argv[1],".db"));
                goto error;
              }
        }

        //   Add the file name to flist

        if (VERBOSE)
          { fprintf(stderr,"Adding '%s' ...\n",core);
            fflush(stderr);
          }
        flist[ofiles++] = core;

        //  Get all the data from the file

#define LOWER_OFFSET 32
#define PHRED_OFFSET 33

        if (intype == IS_BAX)
          { SubRead  *s;
            QVcoding *coding = NULL;
            int       pwell, pcnt;
            int       i, x;
            int64     qpos = 0;

            if (VERBOSE)
              { fprintf(stderr, "  Extracting subreads ...\n"); fflush(stderr); }

            if ((status = getBaxData(bax,Catenate(path,"/",core,".bax.h5"))) != 0)
             { fprintf(stderr, "%s: ", Prog_Name);
                printBaxError(status);
                goto error;
              }

            //  If QUIVER then in a first pass accumulate all the QV statistics and produce
            //    coding tables

            if (QUIVER)
              { if (VERBOSE)
                  { fprintf(stderr, "  Compressing streams ...\n"); fflush(stderr); }

                nextSubread(bax,1);
                QVcoding_Scan1(0,NULL,NULL,NULL,NULL,NULL);
                while ((s = nextSubread(bax,0)) != NULL)
                  { int   rlen;
                    char *delQV, *delTag, *insQV, *mergeQV, *subQV;

                    if ( ! evaluate_bax_filter(EXPR,bax,s))
                      continue;

                    rlen = s->lpulse - s->fpulse;
                    delQV   = bax->delQV   + s->fpulse + s->data_off;
                    delTag  = bax->delTag  + s->fpulse + s->data_off;
                    insQV   = bax->insQV   + s->fpulse + s->data_off;
                    mergeQV = bax->mergeQV + s->fpulse + s->data_off;
                    subQV   = bax->subQV   + s->fpulse + s->data_off;

                    if (isupper(delTag[0]))
                      for (i = 0; i < rlen; i++)
                        delTag[i] += LOWER_OFFSET;
                    x = bax->delLimit;
                    if (isupper(x))
                      x += LOWER_OFFSET;

                    for (i = 0; i < rlen; i++)
                      { if (delQV[i] == x)
                          delTag[i] = 'n';
                        if (delQV[i] > 93)
                          delQV[i] = 126;
                        else
                          delQV[i] += PHRED_OFFSET;
                        if (insQV[i] > 93)
                          insQV[i] = 126;
                        else
                          insQV[i] += PHRED_OFFSET;
                        if (mergeQV[i] > 93)
                          mergeQV[i] = 126;
                        else
                          mergeQV[i] += PHRED_OFFSET;
                        if (subQV[i] > 93)
                          subQV[i] = 126;
                        else
                          subQV[i] += PHRED_OFFSET;
                      }

                    QVcoding_Scan1(rlen,delQV,delTag,insQV,mergeQV,subQV);
                  }

                coding = Create_QVcoding(LOSSY);
                if (coding == NULL)
                  goto error;

                coding->prefix = Strdup(".qvs","Allocating header prefix");
                if (coding->prefix == NULL)
                  goto error;

                qpos = ftello(quiva);
                Write_QVcoding(quiva,coding);
              }

            //  In the penultimate pass, read each entry and accumulate in DB

            if (VERBOSE)
              { fprintf(stderr, "  Transferring data ...\n"); fflush(stderr); }

            pcnt  = 0;
            pwell = -1;
            nextSubread(bax,1);
            while ((s = nextSubread(bax,0)) != NULL)
              { int    rlen, clen;
                char  *read;

                if ( ! evaluate_bax_filter(EXPR,bax,s))
                  continue;

                rlen    = s->lpulse - s->fpulse;
                read    = bax->baseCall + s->fpulse + s->data_off;

                for (i = 0; i < rlen; i++)
                  { x = Number[(int) read[i]];
                    count[x] += 1;
                    read[i] = x;
                  }
                ureads += 1;
                totlen += rlen;
                if (rlen > maxlen)
                  maxlen = rlen;

                prec[pcnt].origin = s->well;
                prec[pcnt].fpulse = s->fpulse;
                prec[pcnt].rlen   = rlen;
                prec[pcnt].boff   = offset;
                prec[pcnt].flags  = s->qv;
                prec[pcnt].coff   = -1;

                Compress_Read(rlen,read);
                clen = COMPRESSED_LEN(rlen);
                fwrite(read,1,clen,bases);

                if (QUIVER)
                  { char  *delQV, *delTag, *insQV, *mergeQV, *subQV;

                    delQV   = bax->delQV   + s->fpulse + s->data_off;
                    delTag  = bax->delTag  + s->fpulse + s->data_off;
                    insQV   = bax->insQV   + s->fpulse + s->data_off;
                    mergeQV = bax->mergeQV + s->fpulse + s->data_off;
                    subQV   = bax->subQV   + s->fpulse + s->data_off;

                    prec[pcnt].coff = qpos;

                    Compress_Next_QVentry1(rlen,delQV,delTag,insQV,
                                           mergeQV,subQV,quiva,coding,LOSSY);
                    qpos = ftello(quiva);
                  }
                if (ARROW)
                  { float  *snr;
                    uint16 *raw;
                    char   *pulse;
                    uint16  cnr[4];

                    raw   = bax->pulseW + s->fpulse + s->data_off;
                    pulse = (char *) raw;

                    for (i = 0; i < rlen; i++)
                      pulse[i] = raw[i]-1;

                    snr = bax->snrVec + 4*s->zmw_off;
                    for (i = 0; i < 4; i++)
                      cnr[i] = (uint32) (snr[bax->chan[i]] * 100.);
                    *((uint64  *) &(prec[pcnt].coff)) = ((uint64) cnr[0]) << 48 |
                                                        ((uint64) cnr[1]) << 32 |
                                                        ((uint64) cnr[2]) << 16 |
                                                        ((uint64) cnr[3]);

                    Compress_Read(rlen,pulse);
                    fwrite(pulse,1,clen,arrow);
                  }

                offset += clen;

                if (pwell == s->well)
                  { prec[pcnt].flags |= DB_CSS;
                    pcnt += 1;
                    if (pcnt >= pmax)
                      { pmax = ((int) (pcnt*1.2)) + 100;
                        prec = (DAZZ_READ *) realloc(prec,sizeof(DAZZ_READ)*pmax);
                        if (prec == NULL)
                          { fprintf(stderr,"%s: Out of memory",Prog_Name);
                            fprintf(stderr," (Allocating %d read records)\n",pmax);
                            goto error;
                          }
                      }
                  }
                else if (pcnt == 0)
                  pcnt += 1;
                else
                  { x = 0;
                    for (i = 1; i < pcnt; i++)
                      if (prec[i].rlen > prec[x].rlen)
                        x = i;
                    prec[x].flags |= DB_BEST;
                    fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);
                    prec[0] = prec[pcnt];
                    pcnt = 1;
                  }
                pwell = s->well;
              }

            //  Complete processing of current file: flush last well group, write file line
            //      in db image, and close file

            x = 0;
            for (i = 1; i < pcnt; i++)
              if (prec[i].rlen > prec[x].rlen)
                x = i;
            prec[x].flags |= DB_BEST;
            fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);
  
            fprintf(ostub,DB_FDATA,ureads,core,bax->movieName);
            ocells += 1;
          }

        else
          { samRecord *rec;
            QVcoding  *coding = NULL;
            int        pwell, pcnt;
            int        i, x;
            int        qpos = 0;
            char      *hdr = NULL;

            if (intype == IS_BAM)
              { if ((input = sam_open(Catenate(path,"/",core,".subreads.bam"))) == NULL)
                  { fprintf(stderr, "%s: can't open %s as a Bam file\n", Prog_Name, ng->name);
                    goto error;
                  }
                fflush(stderr);
              }
            else
              { if ((input = sam_open(Catenate(path,"/",core,".subreads.sam"))) == NULL)
                  { fprintf(stderr, "%s: can't open %s as a Sam file\n", Prog_Name, ng->name);
                    goto error;
                  }
              }

            status = sam_header_process(input,1);
            if (status < 0)
              goto error;
            else if ((status & HASPW) == 0 && ARROW)
              { fprintf(stderr, "%s: %s does not have Arrow information\n", Prog_Name, ng->name);
                goto error;
              }
            else if ((status & HASQV) == 0 && QUIVER)
              { fprintf(stderr, "%s: %s does not have Quiver information\n", Prog_Name, ng->name);
                goto error;
              }

            //  If QUIVER then in a first pass accumulate all the QV statistics and produce
            //    coding tables

            if (QUIVER)
              { if (VERBOSE)
                  { fprintf(stderr, "  Compressing streams ...\n"); fflush(stderr); }

                QVcoding_Scan1(0,NULL,NULL,NULL,NULL,NULL);
                while ((rec = sam_record_extract(input,status)) != SAM_EOF)
                  { char *delQV, *delTag, *insQV, *mergeQV, *subQV;

                    if (rec == NULL)
                      goto error;

                    if ( ! evaluate_bam_filter(EXPR,rec))
                      continue;

                    delQV   = rec->qv[0];
                    delTag  = rec->qv[1];
                    insQV   = rec->qv[2];
                    mergeQV = rec->qv[3];
                    subQV   = rec->qv[4];

                    QVcoding_Scan1(rec->len,delQV,delTag,insQV,mergeQV,subQV);
                  }

                coding = Create_QVcoding(LOSSY);
                if (coding == NULL)
                  goto error;

                coding->prefix = Strdup(".qvs","Allocating header prefix");
                if (coding->prefix == NULL)
                  goto error;

                qpos = ftello(quiva);
                Write_QVcoding(quiva,coding);

                //  Restart input

                sam_close(input);
                if (intype == IS_BAM)
                  input = sam_open(Catenate(path,"/",core,".subreads.bam"));
                else
                  input = sam_open(Catenate(path,"/",core,".subreads.sam"));
                sam_header_process(input,1);
              }

            //  In the penultimate pass, read each entry and accumulate in DB

            if (VERBOSE)
              { fprintf(stderr, "  Transferring data ...\n"); fflush(stderr); }

            pcnt  = 0;
            pwell = -1;
            while ((rec = sam_record_extract(input,status)) != SAM_EOF)
              { int    rlen, clen;
                char  *read;

                if (rec == NULL)
                  goto error;

                if ( ! evaluate_bam_filter(EXPR,rec))
                  continue;

                rlen = rec->len;
                read = rec->seq;

                for (i = 0; i < rlen; i++)
                  count[(int) read[i]] += 1;
                ureads += 1;
                totlen += rlen;
                if (rlen > maxlen)
                  maxlen = rlen;

                prec[pcnt].origin = rec->well;
                prec[pcnt].fpulse = rec->beg;
                prec[pcnt].rlen   = rlen;
                prec[pcnt].boff   = offset;
                prec[pcnt].flags  = (int) (1000.*rec->qual);
                prec[pcnt].coff   = -1;

                hdr = rec->header;

                Compress_Read(rlen,read);
                clen = COMPRESSED_LEN(rlen);
                fwrite(read,1,clen,bases);

                if (QUIVER)
                  { char  *delQV, *delTag, *insQV, *mergeQV, *subQV;

                    delQV   = rec->qv[0];
                    delTag  = rec->qv[1];
                    insQV   = rec->qv[2];
                    mergeQV = rec->qv[3];
                    subQV   = rec->qv[4];

                    prec[pcnt].coff = qpos;

                    Compress_Next_QVentry1(rlen,delQV,delTag,insQV,mergeQV,subQV,quiva,coding,LOSSY);
                    qpos = ftello(quiva);
                  }
                if (ARROW)
                  { char   *pulse;
                    uint16  cnr[4];

                    pulse = rec->arr;

                    for (i = 0; i < 4; i++)
                      cnr[i] = (uint32) (rec->snr[i] * 100.);
                    *((uint64  *) &(prec[pcnt].coff)) = ((uint64) cnr[0]) << 48 |
                                                        ((uint64) cnr[1]) << 32 |
                                                        ((uint64) cnr[2]) << 16 |
                                                        ((uint64) cnr[3]);

                    Compress_Read(rlen,pulse);
                    fwrite(pulse,1,clen,arrow);
                  }

                offset += clen;

                if (pwell == rec->well)
                  { prec[pcnt].flags |= DB_CSS;
                    pcnt += 1;
                    if (pcnt >= pmax)
                      { pmax = ((int) (pcnt*1.2)) + 100;
                        prec = (DAZZ_READ *) realloc(prec,sizeof(DAZZ_READ)*pmax);
                        if (prec == NULL)
                          { fprintf(stderr,"%s: Out of memory",Prog_Name);
                            fprintf(stderr," (Allocating %d read records)\n",pmax);
                            goto error;
                          }
                      }
                  }
                else if (pcnt == 0)
                  pcnt += 1;
               else
                  { x = 0;
                    for (i = 1; i < pcnt; i++)
                      if (prec[i].rlen > prec[x].rlen)
                        x = i;
                    prec[x].flags |= DB_BEST;
                    fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);
                    prec[0] = prec[pcnt];
                    pcnt = 1;
                  }
                pwell = rec->well;
              }

            //  Complete processing of current file: flush last well group, write file line
            //      in db image, and close file

            x = 0;
            for (i = 1; i < pcnt; i++)
              if (prec[i].rlen > prec[x].rlen)
                x = i;
            prec[x].flags |= DB_BEST;
            fwrite(prec,sizeof(DAZZ_READ),pcnt,indx);

            fprintf(ostub,DB_FDATA,ureads,core,hdr);
            ocells += 1;

            sam_close(input);
          }
    
        free(path);
        if (VERBOSE)
          { fprintf(stderr,   "Done\n"); fflush(stdout); }
      }

    //  Finished loading all sequences: update relevant fields in db record

    db.ureads = ureads;
    if (istub == NULL)
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((1.*count[c])/totlen);
        db.totlen = totlen;
        db.maxlen = maxlen;
        db.cutoff = -1;
        if (ARROW)
          db.allarr = DB_ARROW;
        else
          db.allarr = 0;
      }
    else
      { for (c = 0; c < 4; c++)
          db.freq[c] = (float) ((db.freq[c]*db.totlen + (1.*count[c]))/(db.totlen + totlen));
        db.totlen += totlen;
        if (maxlen > db.maxlen)
          db.maxlen = maxlen;
      }
  }

  //  If db has been previously partitioned then calculate additional partition points and
  //    write to new db file image

  if (db.cutoff >= 0)
    { int64      totlen, dbpos, size;
      int        nblock, ireads, tfirst, rlen;
      int        ufirst, cutoff, allflag;
      DAZZ_READ  record;
      int        i;

      if (VERBOSE)
        { fprintf(stderr,"Updating block partition ...\n");
          fflush(stderr);
        }

      //  Read the block portion of the existing db image getting the indices of the first
      //    read in the last block of the exisiting db as well as the partition parameters.
      //    Copy the old image block information to the new block information (except for
      //    the indices of the last partial block)

      if (fscanf(istub,DB_NBLOCK,&nblock) != 1)
        { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
          goto error;
        }
      dbpos = ftello(ostub);
      fprintf(ostub,DB_NBLOCK,0);
      if (fscanf(istub,DB_PARAMS,&size,&cutoff,&allflag) != 3)
        { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
          goto error;
        }
      fprintf(ostub,DB_PARAMS,size,cutoff,allflag);
      if (allflag)
        allflag = 0;
      else
        allflag = DB_BEST;

      nblock -= 1;
      for (i = 0; i <= nblock; i++)
        { if (fscanf(istub,DB_BDATA,&ufirst,&tfirst) != 2)
            { fprintf(stderr,"%s: %s.db is corrupted, read failed\n",Prog_Name,root);
              goto error;
            }
          fprintf(ostub,DB_BDATA,ufirst,tfirst);
        }

      //  Seek the first record of the last block of the existing db in .idx, and then
      //    compute and record partition indices for the rest of the db from this point
      //    forward.

      fseeko(indx,sizeof(DAZZ_DB)+sizeof(DAZZ_READ)*ufirst,SEEK_SET);
      totlen = 0;
      ireads = 0;
      for (i = ufirst; i < ureads; i++)
        { if (fread(&record,sizeof(DAZZ_READ),1,indx) != 1)
            { fprintf(stderr,"%s: %s.idx is corrupted, read failed\n",Prog_Name,root);
              goto error;
            }
          rlen = record.rlen;
          if (rlen >= cutoff && (record.flags & DB_BEST) >= allflag)
            { ireads += 1;
              tfirst += 1;
              totlen += rlen;
              if (totlen >= size)
                { fprintf(ostub," %9d %9d\n",i+1,tfirst);
                  totlen = 0;
                  ireads = 0;
                  nblock += 1;
                }
            }
        }

      if (ireads > 0)
        { fprintf(ostub,DB_BDATA,ureads,tfirst);
          nblock += 1;
        }

      db.treads = tfirst;

      fseeko(ostub,dbpos,SEEK_SET);
      fprintf(ostub,DB_NBLOCK,nblock);    //  Rewind and record the new number of blocks
    }
  else
    db.treads = ureads;

  rewind(indx);
  fwrite(&db,sizeof(DAZZ_DB),1,indx);   //  Write the finalized db record into .idx

  rewind(ostub);                        //  Rewrite the number of files actually added
  fprintf(ostub,DB_NFILE,ocells);

  if (istub != NULL)
    fclose(istub);
  fclose(ostub);
  fclose(indx);
  fclose(bases);
  if (arrow != NULL)
    fclose(arrow);
  if (quiva != NULL)
    fclose(quiva);

  rename(Catenate(pwd,"/",root,".dbx"),dbname);   //  New image replaces old image

  exit (0);

  //  Error exit:  Either truncate or remove the .idx and .bps files as appropriate.
  //               Remove the new image file <pwd>/<root>.dbx

error:
  if (ioff != 0)
    { fseeko(indx,0,SEEK_SET);
      if (ftruncate(fileno(indx),ioff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.idx after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (boff != 0)
    { fseeko(bases,0,SEEK_SET);
      if (ftruncate(fileno(bases),boff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.bps after error, truncate failed\n",
                       Prog_Name,root);
      if (ARROW)
        { fseeko(arrow,0,SEEK_SET);
          if (ftruncate(fileno(arrow),boff) < 0)
            fprintf(stderr,"%s: Fatal: could not restore %s.arw after error, truncate failed\n",
                           Prog_Name,root);
        }
    }
  if (QUIVER && coff != 0)
    { fseeko(quiva,0,SEEK_SET);
      if (ftruncate(fileno(quiva),coff) < 0)
        fprintf(stderr,"%s: Fatal: could not restore %s.qvs after error, truncate failed\n",
                       Prog_Name,root);
    }
  if (indx != NULL)
    { fclose(indx);
      if (ioff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".idx"));
    }
  if (bases != NULL)
    { fclose(bases);
      if (boff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".bps"));
    }
  if (arrow != NULL)
    { fclose(arrow);
      if (boff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".arw"));
    }
  if (quiva != NULL)
    { fclose(quiva);
      if (coff == 0)
        unlink(Catenate(pwd,PATHSEP,root,".qvs"));
    }
  if (ostub != NULL)
    { fclose(ostub);
      unlink(Catenate(pwd,"/",root,".dbx"));
    }
  if (istub != NULL)
    fclose(istub);

  exit (1);
}
