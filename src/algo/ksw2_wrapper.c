#include "ksw2_wrapper.h"

#include "ksw2.h"
#include "hbn_traceback_aux.h"

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	sc_ambi = sc_ambi > 0? -sc_ambi : sc_ambi;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = sc_ambi;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = sc_ambi;
}

Ksw2Data*
Ksw2DataNew()
{
    Ksw2Data* data = (Ksw2Data*)calloc(1, sizeof(Ksw2Data));
#ifdef HAVE_KALLOC
    data->km = km_init();
#else
    data->km = NULL;
#endif
    kv_init(data->qfrag);
    kv_init(data->tfrag);
    return data;
}

void
ksw2_set_params(Ksw2Data* data,
    int reward,
    int penalty,
    int ambi,
    int go,
    int ge,
    int zdrop,
    int band_width)
{
    data->reward = reward;
    data->penalty = penalty;
    data->ambi_penalty = ambi;
    data->go = go;
    data->ge = ge;
    data->go1 = 0;
    data->ge1 = 0;
    data->zdrop = zdrop; // 400;
    data->band_width = band_width; // 2048;
    data->end_bonus = -1;
    ksw_gen_simple_mat(5, data->mat, reward, penalty, ambi);
    data->score_param_is_set = 1;
}

Ksw2Data*
Ksw2DataFree(Ksw2Data* data)
{
#ifdef HAVE_KALLOC
    km_destroy(data->km);
#endif
    kv_destroy(data->qfrag);
    kv_destroy(data->tfrag);
    free(data);
    return NULL;
}

void
ksw2_extd2_set_params(Ksw2Data* data)
{
    data->reward = 2;
    data->penalty = 4;
    data->go = 4;
    data->ge = 2;
    data->go1 = 24;
    data->ge1 = 1;
    data->ambi_penalty = 1;
    data->zdrop = -1;
    data->end_bonus = -1;
    data->band_width = -1;
    ksw_gen_simple_mat(5, data->mat, data->reward, data->penalty, data->ambi_penalty);
}

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
        kstring_t* saln)
{
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = 0;
    if (max_distance == 0) max_distance = hbn_max(qto - qfrom, sto - sfrom) * 0.1;
    int max_band_width = hbn_min(max_distance, 8192);
    const u8* qsubseq = query + qfrom;
    const int qsubseq_size = qto - qfrom;
    const u8* ssubseq = subject + sfrom;
    const int ssubseq_size = sto - sfrom;
    ksw_extd2_sse(data->km, qsubseq_size, qsubseq, ssubseq_size, ssubseq, 5, data->mat,
        data->go, data->ge, data->go1, data->ge1, max_band_width, data->zdrop, data->end_bonus, flag, &ez);
    if (ez.n_cigar == 0) return 0;

    int qi = 0;
    int si = 0;
    ks_clear(*qaln);
    ks_clear(*saln);
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = qsubseq[qi];
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                c = ssubseq[si];
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = qsubseq[qi];
                c = DECODE_RESIDUE(c);
                kputc(c, qaln);
                kputc(GAP_CHAR, saln);
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                kputc(GAP_CHAR, qaln);
                c = ssubseq[si];
                c = DECODE_RESIDUE(c);
                kputc(c, saln);
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }
    if (qi < qsubseq_size || si < ssubseq_size) {
        HBN_LOG("[%d, %d, %d, %d] x [%d, %d, %d, %d] %g%% sw fail, dist = %d", qid, *qoff, *qend, qsize, sid, *soff, *send, ssize, *ident_perc, max_distance);
        return 0;
    }
    hbn_assert(qi == qsubseq_size);
    hbn_assert(si == ssubseq_size);
    kfree(data->km, ez.cigar);

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, qsubseq, 0, qsubseq_size, ks_s(*qaln),
        0, ssubseq, 0, ssubseq_size, ks_s(*saln),
        ks_size(*qaln), TRUE);
    *ident_perc = calc_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), NULL, NULL);
    if (*ident_perc < min_ident_perc) return 0;
    *qoff = qfrom;
    *qend = qto;
    *soff = sfrom;
    *send = sto;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, *qoff, *qend, ks_s(*qaln),
        0, subject, *soff, *send, ks_s(*saln),
        ks_size(*qaln), TRUE);

    return 1;    
}