#ifndef __HBN_HIT_H
#define __HBN_HIT_H

#include "gapped_candidate.h"
#include "seqdb.h"
#include "small_object_alloc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    size_t qoff;
    size_t qend;
    size_t qsize;
    int sid;
    int sdir;
    size_t soff;
    size_t send;
    size_t ssize;
    int score;
    double ident_perc;
    int use_normalised_offset;
    size_t qaln_offset;
    size_t saln_offset;
} HbnHSP;

typedef kvec_t(HbnHSP) vec_hsp;

#define dump_hbn_hsp(output_func, stream, hsp) \
    output_func(stream, \
                "%d\t" \
				"%d\t" \
				"%.2lf\t" \
				"%d\t" \
				"%d\t" \
				"%zu\t" \
                "%zu\t" \
                "%zu\t" \
				"%d\t" \
				"%zu\t" \
                "%zu\t" \
                "%zu\n", \
                (hsp).qid, \
                (hsp).sid, \
                (hsp).ident_perc, \
                (hsp).score, \
                (hsp).qdir, \
                (hsp).qoff, \
                (hsp).qend, \
                (hsp).qsize, \
                (hsp).sdir, \
                (hsp).soff, \
                (hsp).send, \
                (hsp).ssize \
    )

#define dump_hbn_hsp_std(hsp) dump_hbn_hsp(fprintf, stderr, hsp)

#define dump_hbn_hsp_gi(output_func, stream, hsp) \
    output_func(stream, \
                "%s\t" \
				"%s\t" \
				"%.2f\t" \
				"%d\t" \
				"%d\t" \
				"%zu\t" \
                "%zu\t" \
                "%zu\t" \
				"%d\t" \
				"%zu\t" \
                "%zu\t" \
                "%zu\n", \
                qhdr, \
                shdr, \
                (hsp).ident_perc, \
                (hsp).score, \
                (hsp).qdir, \
                (hsp).qoff, \
                (hsp).qend, \
                (hsp).qsize, \
                (hsp).sdir, \
                (hsp).soff, \
                (hsp).send, \
                (hsp).ssize \
    )

#define load_hbn_hsp(input_func, stream, hsp) \
    input_func(stream, \
                "%d" \
				"%d" \
				"%lf" \
				"%d" \
				"%d" \
				"%zu" \
                "%zu" \
                "%zu" \
				"%d" \
				"%zu" \
                "%zu" \
                "%zu", \
                &(hsp).qid, \
                &(hsp).sid, \
                &(hsp).ident_perc, \
                &(hsp).score, \
                &(hsp).qdir, \
                &(hsp).qoff, \
                &(hsp).qend, \
                &(hsp).qsize, \
                &(hsp).sdir, \
                &(hsp).soff, \
                &(hsp).send, \
                &(hsp).ssize \
    )

void ks_introsort_hbn_hsp_score_gt(size_t n, HbnHSP* a);

void ks_introsort_hbn_hsp_ident_gt(size_t n, HbnHSP* a);

void ks_introsort_hbn_hsp_soff_lt(size_t n, HbnHSP* a);

void ks_introsort_hbn_hsp_sid_lt(size_t n, HbnHSP* a);

void normalise_hbn_hsp_offsets(HbnHSP* hsp);

void normalise_hbn_hsp_sdir(HbnHSP* hsp);

typedef struct {
    HbnHSP* hsp_array;
    int oid;
    int hsp_count;
    int best_score;
} HbnHSPList;

void ks_introsort_hbn_hsp_list_score_gt(size_t n, HbnHSPList* a);

typedef struct {
    HbnHSPList* hsplist_array;
    int hsplist_count;
} HbnHitList;

typedef struct {
    HbnHitList* hitlist_array;
    int num_queries;
    int hitlist_max;
    SmallObjectAlloc* hsp_alloc;
    SmallObjectAlloc* hsplist_alloc;
    kstring_t aligned_strings;
    vec_cns_hit cns_hit_list;
    kstring_t output_buf;
} HbnAlignResults;

void
HbnAlignResultsClear(HbnAlignResults* results, int num_queries);

HbnAlignResults*
HbnAlignResultsFree(HbnAlignResults* results);

HbnAlignResults*
HbnAlignResultsNew(const int hitlist_max);

#ifdef __cplusplus
}
#endif

#endif // 