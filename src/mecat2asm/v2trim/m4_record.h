#ifndef M4_RECORD_H
#define M4_RECORD_H

#include <stdio.h>

#include "../klib/kvec.h"
#include "../common/ontcns_defs.h"

typedef struct {
	int qid;
    int qdir;
    idx qoff;
    idx qend;
    idx qsize;
    int sid;
    int sdir;
    idx soff;
    idx send;
    idx ssize;
    double ident_perc;
    int vscore;
} M4Record;

typedef kvec_t(M4Record) vec_m4;

#define m4_ident_gt(a, b) ((a).ident_perc > (b).ident_perc)
void ks_introsort_m4_ident_gt(size_t n, M4Record* m4v);

#define m4_sid_lt(a, b) ((a).sid < (b).sid)
void ks_introsort_m4_sid_lt(size_t n, M4Record* m4v);

#define DUMP_M4_RECORD(output_func, out, m) \
	output_func(out, \
			   "%d\t" \
			   "%d\t" \
			   "%lf\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\t" \
			   "%d\n", \
			   (m).qid, \
			   (m).sid, \
			   (m).ident_perc, \
			   (m).vscore, \
			   (m).qdir, \
			   (m).qoff, \
			   (m).qend, \
			   (m).qsize, \
			   (m).sdir, \
			   (m).soff, \
			   (m).send, \
			   (m).ssize)

#define LOAD_M4_RECORD(input_func, in, m) \
	input_func(in, \
			   "%d" \
			   "%d" \
			   "%lf" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d" \
			   "%d", \
			   &(m).qid, \
			   &(m).sid, \
			   &(m).ident_perc, \
			   &(m).vscore, \
			   &(m).qdir, \
			   &(m).qoff, \
			   &(m).qend, \
			   &(m).qsize, \
			   &(m).sdir, \
			   &(m).soff, \
			   &(m).send, \
			   &(m).ssize)

#endif // M4_RECORD_H
