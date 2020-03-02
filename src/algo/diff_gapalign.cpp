#include "diff_gapalign.h"
#include "hbn_traceback_aux.h"
#include "ksw2.h"

#include <algorithm>

using namespace std;

#define SMALL_DIFF_SEGMENT      500
#define LARGE_DIFF_SEGMENT      1000
#define DIFF_ROW_SIZE           4096
#define DIFF_COLUMN_SIZE        4096
#define SEGMENT_ALIGN_SIZE      4096
#define DIFF_D_PATH_SIZE        5000000
#define DIFF_ALN_PATH_SIZE      5000000

static const int kMatLen = 8;
static const int kMaxOverHang = 1000;

void
DiffAlignParamsInit(DiffAlignParams* params)
{
    params->segment_size = SMALL_DIFF_SEGMENT;
    params->row_size = DIFF_ROW_SIZE;
    params->column_size = DIFF_COLUMN_SIZE;
    params->segment_align_size = SEGMENT_ALIGN_SIZE;
    params->d_path_size = DIFF_D_PATH_SIZE;
    params->aln_path_size = DIFF_ALN_PATH_SIZE;
}

void
DiffAlignmentClear(DiffAlignment* align)
{
    align->aln_str_size = 0;
    align->aln_q_s = 0;
    align->aln_q_e = 0;
    align->aln_t_s = 0;
    align->aln_t_e = 0;
}

DiffAlignment*
DiffAlignmentFree(DiffAlignment* align)
{
    free(align->q_aln_str);
    free(align->t_aln_str);
    free(align);
    return NULL;
}

DiffAlignment*
DiffAlignmentNew(DiffAlignParams* params)
{
    DiffAlignment* align = (DiffAlignment*)calloc(1, sizeof(DiffAlignment));
    align->q_aln_str = (char*)malloc(params->segment_align_size);
    align->t_aln_str = (char*)malloc(params->segment_align_size);
    return align;
}

void
DiffGapAlignDataInit(DiffGapAlignData* data,
    int qoff,
    int qsize,
    int soff,
    int ssize)
{
    int wrk_l = 2 * hbn_max(qsize, ssize);
    ks_reserve(&data->qabuf, wrk_l);
    ks_reserve(&data->sabuf, wrk_l);
    int wrk_ll = 2 * hbn_max(qoff, soff);
    data->qas = data->qae = ks_s(data->qabuf) + wrk_ll;
    data->sas = data->sae = ks_s(data->sabuf) + wrk_ll;
    *data->qae = '\0';
    *data->sae = '\0';
    data->qoff = data->qend = 0;
    data->soff = data->send = 0;
    data->dist = 0;
    data->ident_perc = 0;
}

DiffGapAlignData*
DiffGapAlignDataNew()
{
    DiffGapAlignData* data = (DiffGapAlignData*)calloc(1, sizeof(DiffGapAlignData));
    ks_init(data->qabuf);
    ks_init(data->sabuf);
    DiffAlignParamsInit(&data->params);
    data->dynq = (int*)malloc(sizeof(int) * data->params.row_size);
    data->dynt = (int*)malloc(sizeof(int) * data->params.column_size);
    data->align = DiffAlignmentNew(&data->params);
    data->d_path = (DPathData2*)malloc(sizeof(DPathData2) * data->params.d_path_size);
    data->aln_path = (PathPoint*)malloc(sizeof(PathPoint) * data->params.aln_path_size);
    data->ksw = Ksw2DataNew();
    ksw2_set_params(data->ksw, MATCH_REWARD, MISMATCH_PENALTY,
        AMB_PENALTY, GAP_OPEN, GAP_EXTEND, KSW_ZDROP, KSW_BAND_WIDTH);
    return data;
}

DiffGapAlignData*
DiffGapAlignDataFree(DiffGapAlignData* data)
{
    ks_destroy(data->qabuf);
    ks_destroy(data->sabuf);
    free(data->dynq);
    free(data->dynt);
    DiffAlignmentFree(data->align);
    free(data->d_path);
    free(data->aln_path);
    Ksw2DataFree(data->ksw);
    free(data);
    return NULL;
}

//////////////

struct SCompareDPathData2
{
    bool operator()(const DPathData2& a, const DPathData2& b) {
        return (a.d == b.d) ? (a.k < b.k) : (a.d < b.d);
    }
};

struct DPathDataExtractor
{
    DPathData2* operator()(const int d, const int k) {
        DPathData2 target;
        target.d = d;
        target.k = k;
        return lower_bound(d_path_list, d_path_list + d_path_list_size, target, scmp);
    }

    DPathDataExtractor(DPathData2* d_path_list_, const int d_path_list_size_)
        : d_path_list(d_path_list_),
          d_path_list_size(d_path_list_size_) {}

    SCompareDPathData2 scmp;
    DPathData2* d_path_list;
    const int d_path_list_size;
};

static inline char
extract_char(const u8* A, int i, bool forward)
{
    const int c = forward ? A[i] : A[-i];
    hbn_assert(c >= 0 && c < 4);
    return DECODE_RESIDUE(c);
}

static void
diff_align_get_align_string(const u8* query, const int q_len, 
	const u8* target, const int t_len,
	DPathData2* d_path, const int d_path_size,
	PathPoint* aln_path, DiffAlignment* align,
	const int qe, const int te, 
	int d, int k,
	const int right_extend)
{
	int cd = d;
	int ck = k;
	int aln_path_idx = 0;
	int i;
	DPathDataExtractor dpath_extrator(d_path, d_path_size);
	while (cd >= 0 && aln_path_idx < q_len + t_len + 1)
	{
		DPathData2* d_path_aux = dpath_extrator(cd, ck);
		aln_path[aln_path_idx].x = d_path_aux->x2;
		aln_path[aln_path_idx].y = d_path_aux->y2;
		++aln_path_idx;
		aln_path[aln_path_idx].x = d_path_aux->x1;
		aln_path[aln_path_idx].y = d_path_aux->y1;
		++aln_path_idx;
		ck = d_path_aux->pre_k;
		cd -= 1;
	}
	--aln_path_idx;
	int cx = aln_path[aln_path_idx].x;
	int cy = aln_path[aln_path_idx].y;
	align->aln_q_s = cx;
	align->aln_t_s = cy;
	int aln_pos = 0;
	while (aln_path_idx > 0)
	{
		--aln_path_idx;
		int nx = aln_path[aln_path_idx].x;
		int ny = aln_path[aln_path_idx].y;
		if (cx == nx && cy == ny) continue;
		
		if (cx == nx && cy != ny) {
			for (i = 0; i < ny - cy; ++i) {
				align->q_aln_str[aln_pos + i] = GAP_CHAR;
				align->t_aln_str[aln_pos + i] = extract_char(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		} else if (cx != nx && cy == ny) {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char(query, cx + i, right_extend);
				align->t_aln_str[aln_pos + i] = GAP_CHAR;
			}
			aln_pos += nx - cx;
		} else {
			for (i = 0; i < nx - cx; ++i) {
				align->q_aln_str[aln_pos + i] = extract_char(query, cx + i, right_extend);
			}
			for (i = 0; i < ny - cy; ++i) {
				align->t_aln_str[aln_pos + i] = extract_char(target, cy + i, right_extend);
			}
			aln_pos += ny - cy;
		}
		
		cx = nx;
		cy = ny;
	}
	align->aln_str_size = aln_pos;
}

static int
diff_align(const u8* query, const int q_len, const u8* target, const int t_len,
	const int band_tolerance, const int get_aln_str, DiffAlignment* align,
	int* V, int* U, DPathData2* d_path, PathPoint* aln_path, const int right_extend)
{
	int k_offset;
	int  d;
	int  k, k2;
	int best_m;
	int min_k, new_min_k, max_k, new_max_k, pre_k;
	int x = -1, y = -1;
	int max_d, band_size;
	unsigned long d_path_idx = 0, max_idx = 0;
	int aligned = 0;
	int best_x = -1, best_y = -1, best_d = q_len + t_len + 100, best_k = 0, best_d_path_idx = -1;

	max_d = (int)(.3 * (q_len + t_len));
	k_offset = max_d;
	band_size = band_tolerance * 2;
    DiffAlignmentClear(align);
	best_m = -1;
	min_k = 0;
	max_k = 0;
	d_path_idx = 0;
	max_idx = 0;

	for (d = 0; d < max_d; ++d)
	{
		if (max_k - min_k > band_size) break;
		
		for (k = min_k; k <= max_k; k += 2)
		{
			if( k == min_k || (k != max_k && V[k - 1 + k_offset] < V[k + 1 + k_offset]) )
			{ pre_k = k + 1; x = V[k + 1 + k_offset]; }
			else
			{ pre_k = k - 1; x = V[k - 1 + k_offset] + 1; }
			y = x - k;
			d_path[d_path_idx].d = d;
			d_path[d_path_idx].k = k;
			d_path[d_path_idx].x1 = x;
			d_path[d_path_idx].y1 = y;

			if (right_extend)
				while( x < q_len && y < t_len && query[x] == target[y]) { ++x; ++y; }
			else
				while( x < q_len && y < t_len && query[-x] == target[-y]) { ++x; ++y; }

			d_path[d_path_idx].x2 = x;
			d_path[d_path_idx].y2 = y;
			d_path[d_path_idx].pre_k = pre_k;
			++d_path_idx;

			V[k + k_offset] = x;
			U[k + k_offset] = x + y;
			if (x + y > best_m) {
				best_m = x + y;
				best_x = x;
				best_y = y;
				best_d = d;
				best_k = k;
				best_d_path_idx = d_path_idx;
			}
			if (x >= q_len || y >= t_len)
			{ aligned = 1; max_idx = d_path_idx; break; }
		}

		// for banding
		new_min_k = max_k;
		new_max_k = min_k;
		for (k2 = min_k; k2 <= max_k; k2 += 2)
			if (U[k2 + k_offset] >= best_m - band_tolerance)
			{ new_min_k = std::min(new_min_k, k2); new_max_k = std::max(new_max_k, k2); }
		max_k = new_max_k + 1;
		min_k = new_min_k - 1;

		if (aligned)
		{
			align->aln_q_e = x;
			align->aln_t_e = y;
			align->dist = d;
			align->aln_str_size = (x + y + d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;

			if (get_aln_str) {
				diff_align_get_align_string(query, q_len, 
                    target, t_len, 
                    d_path, max_idx, aln_path, align, x, y, d, k, right_extend);
			}
            //HBN_LOG("aligned, right_extend = %d", right_extend);
			break;
		} 
	}
	
	if (!aligned) {
        //HBN_LOG("not aligned, best_x = %d", best_x);
		if (best_x > 0) {
			align->aln_q_e = best_x;
			align->aln_t_e = best_y;
			align->dist = best_d;
			align->aln_str_size = (best_x + best_y + best_d) / 2;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
			if (get_aln_str) {
				diff_align_get_align_string(query, q_len, 
                    target, t_len, 
                    d_path, best_d_path_idx, aln_path, align, best_x, best_y, best_d, best_k, right_extend);
			}
		} else {
			align->aln_q_e = 0;
			align->aln_t_e = 0;
			align->dist = 0;
			align->aln_str_size = 0;
			align->aln_q_s = 0;
			align->aln_t_s = 0;
		}
	}
	
	return (align->aln_q_e == q_len || align->aln_t_e == t_len);
}

static bool
retrieve_next_aln_block(const u8* query,
						int qidx, 
						const int qsize, 
						const u8* target,
						int tidx, 
						const int tsize, 
						const int desired_block_size, 
						const bool forward,
						const u8*& Q,
						const u8*& T,
						int& qblk, 
						int& tblk)
{
	bool last_block;
	int qleft = qsize - qidx;
	int tleft = tsize - tidx;
	if (qleft < desired_block_size + 100 || tleft < desired_block_size + 100) {
		qblk = min(qleft, static_cast<int>(tleft + tleft * 0.2));
		tblk = min(tleft, static_cast<int>(qleft + qleft * 0.2));
		last_block = true;
	} else {
		qblk = desired_block_size;
		tblk = desired_block_size;
		last_block = false;
	}
	
	if (forward) {
		Q = query + qidx;
		T = target + tidx;
	} else {
		Q = query - qidx;
		T = target - tidx;
	}
	
	return last_block;
}

static void
dump_subseq(const u8* seq, int len)
{
    for (int i = 0; i < len; ++i) {
        char c = DECODE_RESIDUE(seq[i]);
        fprintf(stderr, "%c", c);
    }
    fprintf(stderr, "\n");
}

static void
dw_extend(const u8* query, const int query_size,
    const u8* target, const int target_size,
    int* U, int* V,
    DiffAlignment* align,
    DPathData2* d_path,
    PathPoint* aln_path,
    DiffAlignParams* params,
    const int right_extend,
    int* qend,
    int* tend,
    char** qas,
    char** qae,
    char** sas,
    char** sae)
{
    *qend = 0;
    *tend = 0;
    const int kBlkSize = params->segment_size;
	int qidx = 0, tidx = 0;
	int qblk, tblk;
	const u8* seq1;
	const u8* seq2;
	while (1) {
		fill(U, U + params->row_size, 0);
		fill(V, V + params->column_size, 0);
		bool last_block = retrieve_next_aln_block(query,
												  qidx,
												  query_size,
												  target,
												  tidx,
												  target_size,
												  kBlkSize,
												  right_extend,
												  seq1,
												  seq2,
												  qblk,
												  tblk);	
        //dump_subseq(seq1, qblk);
        //dump_subseq(seq2, tblk);  
		diff_align(seq1,
			qblk,
			seq2,
			tblk,
			0.3 * max(qblk, tblk),
			400,
			align,
			U, 
			V,
			d_path,
			aln_path,
			right_extend);
        //HBN_LOG("align_size = %d", align->aln_str_size);

#if 1
        validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            0, seq1, 0, align->aln_q_e, align->q_aln_str,
            0, seq2, 0, align->aln_t_e, align->t_aln_str,
            align->aln_str_size, right_extend);
#endif 
		//HBN_LOG("last_block = %d", last_block);
        int done = last_block;
        int acnt = 0, qcnt = 0, tcnt = 0;
        //HBN_LOG("qblk = %d, qe = %d, tblk = %d, te = %d", qblk, align->aln_q_e, tblk, align->aln_t_e);
        if (qblk - align->aln_q_e > 30 && tblk - align->aln_t_e > 30) done = 1;
        int align_size = align->aln_str_size;
        int k = align_size - 1, m = 0;
        while (k >= 0) {
            const char qc = align->q_aln_str[k];
            const char tc = align->t_aln_str[k];
            if (qc != GAP_CHAR) ++qcnt;
            if (tc != GAP_CHAR) ++tcnt;
            m = (qc == tc) ? (m+1) : 0;
            ++acnt;
            if (m == kMatLen) break;
            --k;
        }

//HBN_LOG("qidx = %d, tidx = %d, align_size = %d, qcnt = %d, tcnt = %d, m = %d, k = %d, done = %d",
//    qidx, tidx, align_size, qcnt, tcnt, m, k, done);

        if (m != kMatLen || k < 1) {
            align_size = 0;
            for (int i = 0; i < qblk && i < tblk; ++i) {
                const char qc = extract_char(seq1, i, right_extend);
                const char tc = extract_char(seq2, i, right_extend);
                if (qc != tc) break;
                align->q_aln_str[align_size] = qc;
                align->t_aln_str[align_size] = tc;
                ++align_size;
            }
            qidx += align_size;
            tidx += align_size;
            done = 1;
        } else {
            align_size -= acnt;
            qidx += (align->aln_q_e - qcnt);
            tidx += (align->aln_t_e - tcnt);
            if (done) {
                align_size += kMatLen;
                qidx += kMatLen;
                tidx += kMatLen;
            }
        }

//HBN_LOG("qidx = %d, tidx = %d, align_size = %d, qcnt = %d, tcnt = %d, m = %d, k = %d, done = %d",
//    qidx, tidx, align_size, qcnt, tcnt, m, k, done);

        if (right_extend) {
            for (int i = 0; i < align_size; ++i) {
                **qae = align->q_aln_str[i];
                ++(*qae);
                **sae = align->t_aln_str[i];
                ++(*sae);
            }
            **qae = '\0';
            **sae = '\0';
        } else {
            for (int i = 0; i < align_size; ++i) {
                --(*qas);
                **qas = align->q_aln_str[i];
                --(*sas);
                **sas = align->t_aln_str[i];
            }
        }
        if (done) break;
	}
    *qend = qidx;
    *tend = tidx;
}

static BOOL
diff_align_trim_left_end(DiffGapAlignData* data)
{
    //HBN_LOG("trim left end");
    if (data->qas == data->qae) {
        hbn_assert(data->sas == data->sae);
        hbn_assert(data->qoff == data->qend);
        hbn_assert(data->soff == data->send);
        return TRUE;
    }

    hbn_assert(data->qae - data->qas == data->sae - data->sas);
    int qcnt = 0, scnt = 0, acnt = 0, m = 0;
    const char* qaln = data->qas;
    const char* saln = data->sas;
    while (qaln < data->qae) {
        const char qc = *qaln;
        const char sc = *saln;
        ++qaln;
        ++saln;
        if (qc != GAP_CHAR) ++qcnt;
        if (sc != GAP_CHAR) ++scnt;
        ++acnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen) return FALSE;
    qcnt -= kMatLen;
    scnt -= kMatLen;
    acnt -= kMatLen;
    data->qoff += qcnt;
    data->soff += scnt;
    data->qas += acnt;
    data->sas += acnt;
    hbn_assert(data->qoff < data->qend);
    hbn_assert(data->soff < data->send);
    hbn_assert(data->qas < data->qae);
    hbn_assert(data->sas < data->sae);
    return TRUE;
}

static BOOL
diff_align_trim_right_end(DiffGapAlignData* data)
{
    if (data->qas == data->qae) {
        hbn_assert(data->sas == data->sae);
        hbn_assert(data->qoff == data->qend);
        hbn_assert(data->soff == data->send);
        return TRUE;
    }    
    hbn_assert(data->qae - data->qas == data->sae - data->sas);
    hbn_assert((*data->qae) == '\0');
    hbn_assert((*data->sae) == '\0');
    //DiffGapAlignDataDump(fprintf, stderr, data);

    int qcnt = 0, scnt = 0, acnt = 0, m = 0;
    const char* qaln = data->qae;
    const char* saln = data->sae;
    while (qaln > data->qas) {
        --qaln;
        --saln;
        const char qc = *qaln;
        const char sc = *saln;
        if (qc != GAP_CHAR) ++qcnt;
        if (sc != GAP_CHAR) ++scnt;
        ++acnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen) return FALSE;
    qcnt -= kMatLen;
    scnt -= kMatLen;
    acnt -= kMatLen;
    int org_qend = data->qend;
    int org_send = data->send;
    data->qend -= qcnt;
    data->send -= scnt;
    data->qae -= acnt;
    data->sae -= acnt;
    *data->qae = '\0';
    *data->sae = '\0';
    //DiffGapAlignDataDump(fprintf, stderr, data);
    hbn_assert(data->qoff < data->qend, 
        "qid = %d, sid = %d, qoff = %d, org_qend = %d, qend = %d, soff = %d, org_send = %d, send = %d, qcnt = %d, scnt = %d, acnt = %d, m = %d",
        data->qid, data->sid, data->qoff, org_qend, data->qend, data->soff, org_send, data->send,
        qcnt, scnt, acnt, m);
    hbn_assert(data->soff < data->send);
    hbn_assert(data->qas < data->qae);
    hbn_assert(data->sas < data->sae);
    return TRUE;
}

static int
left_extend(Ksw2Data* ksw,
    const u8* query,
    const u8* subject,
    vec_u8* qsbuf,
    vec_u8* ssbuf,
    int* qbeg,
    int* sbeg,
    char** qas,
    char** sas)
{
    int qls = *qbeg;
    int sls = *sbeg;
    int ls = hbn_min(qls, sls);
    if (ls > kMaxOverHang) return 0;
    if (ls == 0) return 0;

    qls += kMatLen;
    sls += kMatLen;

    kv_clear(*qsbuf);
    kv_clear(*ssbuf);
    for (ls = qls; ls; --ls) kv_push(u8, *qsbuf, query[ls-1]);
    for (ls = sls; ls; --ls) kv_push(u8, *ssbuf, subject[ls-1]);
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = KSW_EZ_RIGHT | KSW_EZ_EXTZ_ONLY;
    ksw_extz2_sse(ksw->km, qls, kv_data(*qsbuf), sls, kv_data(*ssbuf), 5, ksw->mat,
        ksw->go, ksw->ge, ksw->band_width, ksw->zdrop, ksw->end_bonus, flag, &ez);
    if (ez.n_cigar == 0) return 0;

    *qas += kMatLen;
    *sas += kMatLen;
    *qbeg += kMatLen;
    *sbeg += kMatLen;

    int qi = 0;
    int si = 0;
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = kv_A(*qsbuf, qi);
                c = DECODE_RESIDUE(c);
                --(*qas); **qas = c;
                c = kv_A(*ssbuf, si);
                c = DECODE_RESIDUE(c);
                --(*sas); **sas = c; 
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = kv_A(*qsbuf, qi);
                c = DECODE_RESIDUE(c);
                --(*qas); **qas = c;
                --(*sas); **sas = GAP_CHAR;
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                --(*qas); **qas = GAP_CHAR;
                c = kv_A(*ssbuf, si);
                c = DECODE_RESIDUE(c);
                --(*sas); **sas = c; 
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }
    
    kfree(ksw->km, ez.cigar);
    *qbeg -= qi;
    *sbeg -= si;
    return 1;
}

static int
right_extend(Ksw2Data* ksw,
    const u8* query,
    int* qend,
    const int query_length,
    const u8* subject,
    int* send,
    const int subject_length,
    char** qae,
    char** sae)
{
    int qrs = query_length - (*qend);
    int srs = subject_length - (*send);
    int rs = hbn_min(qrs, srs);
    if (rs > kMaxOverHang || rs == 0) return 0;

    qrs += kMatLen;
    srs += kMatLen;
    const u8* q = query + query_length - qrs;
    const u8* s = subject + subject_length - srs;
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = KSW_EZ_RIGHT | KSW_EZ_EXTZ_ONLY;
    ksw_extz2_sse(ksw->km, qrs, q, srs, s, 5, ksw->mat,
        ksw->go, ksw->ge, ksw->band_width, ksw->zdrop, ksw->end_bonus, flag, &ez);
    if (ez.n_cigar == 0) return 0;    

    *qend -= kMatLen;
    *send -= kMatLen;
    *qae -= kMatLen;
    *sae -= kMatLen;
    int qi = 0;
    int si = 0;
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = query[*qend + qi];
                c = DECODE_RESIDUE(c);
                **qae = c; ++(*qae);
                c = subject[*send + si];
                c = DECODE_RESIDUE(c);
                **sae = c; ++(*sae);
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = query[*qend + qi];
                c = DECODE_RESIDUE(c);
                **qae = c; ++(*qae);
                **sae = GAP_CHAR; ++(*sae);
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                **qae = GAP_CHAR; ++(*qae);
                c = subject[*send + si];
                c = DECODE_RESIDUE(c);
                **sae = c; ++(*sae);
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }    

    kfree(ksw->km, ez.cigar);
    **qae = '\0';
    **sae = '\0';
    *qend += qi;
    *send += si;
    return 1;
}

int
diff_align(DiffGapAlignData* data,
    const u8* query,
    int qoff,
    int qsize,
    const u8* target,
    int toff,
    int tsize,
    const int min_align_size,
    const double min_ident_perc,
    const BOOL process_over_hang)
{
    DiffGapAlignDataInit(data, qoff, qsize, toff, tsize);
    hbn_assert(data->qas == data->qae);
    hbn_assert(data->sas == data->sae);
 
    int qend;
    int tend;
    dw_extend(query + qoff - 1,
        qoff,
        target + toff - 1,
        toff,
        data->dynq,
        data->dynt,
        data->align,
        data->d_path,
        data->aln_path,
        &data->params,
        FALSE,
        &qend,
        &tend,
        &data->qas,
        &data->qae,
        &data->sas,
        &data->sae);
    //DiffGapAlignDataDump(fprintf, stderr, data);
    //HBN_LOG("qend = %d, tend = %d", qend, tend);
    hbn_assert(qend <= qoff);
    hbn_assert(tend <= toff);
    data->qoff = qoff - qend;
    data->qend = qoff;
    data->soff = toff - tend;
    data->send = toff;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, data->qoff, data->qend, data->qas,
        0, target, data->soff, data->send, data->sas,
        data->qae - data->qas, TRUE);
    int left_extend_success = 1;
    if (!diff_align_trim_right_end(data)) {
        DiffGapAlignDataInit(data, qoff, qsize, toff, tsize);
        left_extend_success = 0;
    }
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, data->qoff, data->qend, data->qas,
        0, target, data->soff, data->send, data->sas,
        data->qae - data->qas, TRUE);

    dw_extend(query + data->qend,
        qsize - data->qend,
        target + data->send,
        tsize - data->send,
        data->dynq,
        data->dynt,
        data->align,
        data->d_path,
        data->aln_path,
        &data->params,
        TRUE,
        &qend,
        &tend,
        &data->qas,
        &data->qae,
        &data->sas,
        &data->sae);
    data->qend += qend;
    data->send += tend;
    hbn_assert(data->qend <= qsize);
    hbn_assert(data->send <= tsize);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, data->qoff, data->qend, data->qas,
        0, target, data->soff, data->send, data->sas,
        data->qae - data->qas, TRUE);
    if ((!left_extend_success) && (!diff_align_trim_left_end(data))) {
        return FALSE;
    }
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, data->qoff, data->qend, data->qas,
        0, target, data->soff, data->send, data->sas,
        data->qae - data->qas, TRUE);

    data->ident_perc = calc_ident_perc(data->qas, data->sas, data->qae - data->qas, &data->dist, &data->score);
    if (data->ident_perc < min_ident_perc - 2.0) return FALSE;
    if (process_over_hang) {
        left_extend(data->ksw, query, target, &data->ksw->qfrag, &data->ksw->tfrag,
            &data->qoff, &data->soff, &data->qas, &data->sas);
        right_extend(data->ksw, query, &data->qend, qsize,
            target, &data->send, tsize, &data->qae, &data->sae);
    }   
    data->ident_perc = calc_ident_perc(data->qas, data->sas, data->qae - data->qas, &data->dist, &data->score);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, data->qoff, data->qend, data->qas,
        0, target, data->soff, data->send, data->sas,
        data->qae - data->qas, TRUE);

    int r = (data->qae - data->qas >= min_align_size) && (data->ident_perc >= min_ident_perc);
    return r;
}