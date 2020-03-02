#ifndef KSW2_WRAPPER_H
#define KSW2_WRAPPER_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void* km;
    int reward, penalty, ambi_penalty;
    int go, ge;
    int go1, ge1;
    int zdrop;
    int band_width;
    int end_bonus;
    int8_t mat[25];
    int score_param_is_set;
    vec_u8 qfrag;
    vec_u8 tfrag;
} Ksw2Data;

Ksw2Data*
Ksw2DataNew();

Ksw2Data*
Ksw2DataFree(Ksw2Data* data);

void
ksw2_set_params(Ksw2Data* data,
    int reward,
    int penalty,
    int ambi,
    int go,
    int ge,
    int zdrop,
    int band_width);

int
ksw2_align(Ksw2Data* data,
    const u8* query,
    int qoff,
    int qsize,
    const u8* target,
    int toff,
    int tsize,
    const int min_align_size,
    const double min_ident_perc,
	int* qbeg,
	int* qend,
	int* tbeg,
	int* tend,
	double* ident_perc,
	kstring_t* query_align,
	kstring_t* target_align);

void
ksw2_extd2_set_params(Ksw2Data* data);

int nw_ksw2_extd2(Ksw2Data* data,
        const int qid,
        const u8* query,
        const int qfrom,
        const int qto,
        const int qsize,
        const int sid,
        const u8* subject,
        const int sfrom,
        const int sto,
        const int ssize,
        const int min_align_size,
        const double min_ident_perc,
        int max_distance,
        int* qoff,
        int* qend,
        int* soff,
        int* send,
        double* ident_perc,
        kstring_t* qaln,
        kstring_t* saln);

#ifdef __cplusplus
}
#endif
#endif // KSW2_WRAPPER_H