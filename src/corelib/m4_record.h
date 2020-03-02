#ifndef __M4_RECORD_H
#define __M4_RECORD_H

#include "hbn_aux.h"
#include "kvec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int qid;
    int qdir;
    int qoff;
    int qend;
    int qsize;
    int sid;
    int sdir;
    int soff;
    int send;
    int ssize;
    double ident_perc;
    int score;
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
			   (m).score, \
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
			   &(m).score, \
			   &(m).qdir, \
			   &(m).qoff, \
			   &(m).qend, \
			   &(m).qsize, \
			   &(m).sdir, \
			   &(m).soff, \
			   &(m).send, \
			   &(m).ssize)

#ifdef __cplusplus
}
#endif

#endif // __M4_RECORD_H