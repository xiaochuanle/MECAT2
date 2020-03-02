#include "cns_one_read.h"

#include "../../algo/hbn_traceback_aux.h"

#include <algorithm>
#include <sstream>

using namespace std;

static BOOL
meap_consensus_one_segment(FCCnsData* cns_data,
    int *sfrom,
    int *sto,
    const int subject_size,
    const int min_size,
    kstring_t* cns_seq)
{
    build_backbone(kv_data(cns_data->tag_list),
        kv_size(cns_data->tag_list),
        subject_size,
        cns_data->dci_alloc,
        cns_data->li_alloc,
        &cns_data->item_list,
        &cns_data->cov_list);
    
    consensus_backbone_segment(kv_data(cns_data->item_list),
        *sfrom,
        *sto,
        kv_data(cns_data->cov_list),
        cns_seq,
        sfrom,
        sto);

    if (ks_size(*cns_seq) < min_size) return FALSE;
    for (size_t p = 0; p < ks_size(*cns_seq); ++p) {
        int c = ks_A(*cns_seq, p);
        c = DECODE_RESIDUE(c);
        ks_A(*cns_seq, p) = c;
    }    
    return TRUE;
}

struct CmpMappingRangeBySoff
{
	bool operator()(const MappingRange& a, const MappingRange& b)
	{
		return (a.start == b.start) ? (a.end > b.end) : (a.start < b.start);
	}
};

void 
get_effective_ranges(MappingRange* m_cov_array, 
    const int m_cov_count, 
    MappingRange* e_cov_array,
    int*  e_cov_count,
    const int read_size, 
    const int min_size)
{
   // HBN_LOG("cov_count = %d", m_cov_count);
    if (m_cov_count == 0) {
        *e_cov_count = 0;
        return;
    }
    for (int i = 0; i < m_cov_count; ++i) {
        if (m_cov_array[i].start <= 500 && read_size - m_cov_array[i].end <= 500) {
            e_cov_array[0].start = 0;
            e_cov_array[0].end = read_size;
            *e_cov_count = 1;
            return;
        }
    }

	std::sort(m_cov_array, m_cov_array + m_cov_count, CmpMappingRangeBySoff());
	int i = 0, j;
    int e_cov_idx = 0;
	int left = m_cov_array[i].start, right;
	while (i < m_cov_count) {
		j = i + 1;
		while (j < m_cov_count && m_cov_array[j].end <= m_cov_array[i].end) ++j;
		if (j == m_cov_count) {
			right = m_cov_array[i].end;
			if (right - left >= min_size * 0.95) {
                e_cov_array[e_cov_idx].start = left;
                e_cov_array[e_cov_idx].end = right;
                ++e_cov_idx;
            }
			break;
		}
		if (m_cov_array[i].end - m_cov_array[j].start < 1000) {
			right = std::min(m_cov_array[i].end, m_cov_array[j].start);
			if (right - left >= min_size * 0.95) {
                e_cov_array[e_cov_idx].start = left;
                e_cov_array[e_cov_idx].end = right;
                ++e_cov_idx;                
            }
			left = std::max(m_cov_array[i].end, m_cov_array[j].start);
		}
		i = j;
	}
    *e_cov_count = e_cov_idx;
}

#ifdef NDEBUG
static inline bool
#else
static bool
#endif
check_ovlp_mapping_range(const int qb, const int qe, const int qs,
						 const int sb, const int se, const int ss,
						 double ratio)
{
	const int oq = qe - qb;
	const int qqs = qs * ratio;
	const int os = se - sb;
	const int qss = ss * ratio;
	return oq >= qqs || os >= qss;
}

static bool
subject_subseq_cov_is_full(u8* cov_stats, int soff, int send)
{
    int n = 0;
    for (int i = soff; i < send; ++i) 
        if (cov_stats[i] >= MAX_CNS_COV) ++n;
    if (send - soff >= n + 200) return false;
    return true;
}

static bool
extend_cns_hit(DiffGapAlignData* diff_data,
    Ksw2Data* ksw,
    const HbnProgramOptions* opts,
    const HbnConsensusInitHit* hit,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_length,
    const u8* fwd_subject,
    const int subject_length,
    u8* cov_stats,
    kstring_t* qaln,
    kstring_t* saln,
    int* qoff_,
    int* qend_,
    int* soff_,
    int* send_,
    double* ident_perc_)
{
    const u8* read = (hit->strand == 1) ? fwd_read : rev_read;
    int read_gapped_start = (hit->strand == 1) ? (hit->qoff) : (read_length - 1 - hit->qoff);
    diff_data->qsize = read_length;
    diff_data->ssize = subject_length;
    int r = diff_align(diff_data,
                read,
                read_gapped_start,
                read_length,
                fwd_subject,
                hit->soff,
                subject_length,
                opts->ovlp_cov_res,
                opts->perc_identity,
                FALSE);
    //DiffGapAlignDataDump(fprintf, stderr, diff_data);
    if (!r) return false;
    hbn_assert(diff_data->qae - diff_data->qas == diff_data->sae - diff_data->sas);
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        hit->qid, read, diff_data->qoff, diff_data->qend, diff_data->qas,
        hit->sid, fwd_subject, diff_data->soff, diff_data->send, diff_data->sas,
        diff_data->qae - diff_data->qas, TRUE);
    int qoff = diff_data->qoff;
    int qend = diff_data->qend;
    int soff = diff_data->soff;
    int send = diff_data->send;
    if (subject_subseq_cov_is_full(cov_stats, soff, send)) return FALSE;
    if (!check_ovlp_mapping_range(qoff, qend, read_length,
            soff, send, subject_length, opts->ovlp_cov_perc / 100.0)) return false;
    normalize_gaps(diff_data->qas, diff_data->sas, diff_data->qae - diff_data->qas, qaln, saln, TRUE);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        hit->qid, read, qoff, qend, ks_s(*qaln),
        hit->sid, fwd_subject, soff, send, 
        ks_s(*saln), ks_size(*qaln), TRUE);
    *qoff_ = qoff;
    *qend_ = qend;
    *soff_ = soff;
    *send_ = send;
    *ident_perc_ = calc_ident_perc(ks_s(*qaln), ks_s(*saln), ks_size(*qaln), NULL, NULL);
    for (int i = soff; i < send; ++i) ++cov_stats[i];
    return true; 
}

extern "C"
void
consensus_one_read(CnsThreadData* data, const int cns_info_idx)
{
    hbn_assert(cns_info_idx < data->cns_info_count);
    RawReadCnsInfo* cns_info = data->cns_info_array + cns_info_idx;
    HbnConsensusInitHit* cns_hit_array = data->cns_hit_array + cns_info->can_from;
    int cns_hit_count = cns_info->can_to - cns_info->can_from;
    const HbnProgramOptions* opts = data->opts;
    if (cns_hit_count < opts->min_cov) return;
    ks_introsort_cns_hit_score_gt(cns_hit_count, cns_hit_array);
    if (cns_hit_count > MAX_CNS_OVLPS) cns_hit_count = MAX_CNS_OVLPS;
    const int subject_id = cns_hit_array[0].sid;
    const int subject_length = data->raw_reads->seqinfo_array[subject_id].seq_size;
    const char* subject_name = data->raw_reads->seq_names
                               + data->raw_reads->seqinfo_array[subject_id].hdr_offset;
    vec_u8* fwd_subject_v = &data->fwd_subject;
    RawReadsReaderExtractRead(data->raw_reads, subject_id, FWD, fwd_subject_v);
    const u8* fwd_subject = kv_data(*fwd_subject_v);
    vec_u8* rev_subject_v = &data->rev_subject;
    RawReadsReaderExtractRead(data->raw_reads, subject_id, REV, rev_subject_v);
    const u8* rev_subject = kv_data(*rev_subject_v);
    vec_u8* cov_stats_v = &data->cov_stats;
    kv_resize(u8, *cov_stats_v, subject_length);
    vec_u8* fwd_read_v = &data->fwd_read;
    vec_u8* rev_read_v = &data->rev_read;
    kstring_t* qaln = &data->qaux;
    kstring_t* saln = &data->saux;
    FCCnsData* cns_data = data->cns_data;
    int num_extended_can = 0;
    int num_added_aln = 0;
    CnsThreadDataInit(data, fwd_subject, rev_subject, subject_length);
    u8* cov_stats = kv_data(data->cov_stats);
    MappingRange m_ovlp_cov_array[MAX_CNS_OVLPS];
    int m_ovlp_cov_count = 0;

    //const char* target = "m130605_231954_42210_c100515112550000001823076608221304_s1_p0/25145/0_17258";
    //if (strcmp(target, subject_name)) return;
    //HBN_LOG("correcting %d:%s:%d:%d", subject_id, subject_name, subject_length, cns_hit_count);

    for (int i = 0; i < cns_hit_count; ++i) {
        const HbnConsensusInitHit* hit = cns_hit_array + i;
        //if (hit->qid != 4924) continue;
        ++num_extended_can;
        //dump_cns_hit(fprintf, stderr, *hit);
        hbn_assert(hit->sid == subject_id);
        RawReadsReaderExtractRead(data->raw_reads, hit->qid, FWD, fwd_read_v);
        const u8* fwd_read = kv_data(*fwd_read_v);
        RawReadsReaderExtractRead(data->raw_reads, hit->qid, REV, rev_read_v);
        const u8* rev_read = kv_data(*rev_read_v);
        const int read_length = data->raw_reads->seqinfo_array[hit->qid].seq_size;
        hbn_assert(read_length == kv_size(*fwd_read_v));
        //HBN_LOG("read_length = %d", read_length);
        int qoff, qend, soff, send;
        double ident_perc;
        if (!extend_cns_hit(data->diff_data,
                data->ksw,
                opts,
                hit,
                fwd_read,
                rev_read,
                read_length,
                fwd_subject,
                subject_length,
                cov_stats,
                qaln,
                saln,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc)) continue;
        //HBN_LOG("add %d [%d, %d, %d] x [%d, %d, %d], %g", 
        //    i, qoff, qend, read_length, soff, send, subject_length, ident_perc);
        ++num_added_aln;
        make_align_tags_from_ovlp(ks_s(*qaln),
            ks_s(*saln),
            ks_size(*qaln),
            qoff,
            qend,
            soff,
            send,
            DEFAULT_CNS_WEIGHT,
            &cns_data->tag_list);
        m_ovlp_cov_array[m_ovlp_cov_count].start = soff;
        m_ovlp_cov_array[m_ovlp_cov_count].end = send;
        m_ovlp_cov_count++;
        if (num_added_aln >= MAX_CNS_COV && subject_subseq_cov_is_full(cov_stats, 0, subject_length)) break;
    }

    MappingRange e_ovlp_cov_array[MAX_CNS_OVLPS];
    int e_ovlp_cov_count = 0;
    get_effective_ranges(m_ovlp_cov_array,
        m_ovlp_cov_count,
        e_ovlp_cov_array,
        &e_ovlp_cov_count,
        subject_length,
        data->opts->min_size);
    if (e_ovlp_cov_count == 0) return;
    int max_size = 0, max_i = -1;
    for (int i = 0; i < e_ovlp_cov_count; ++i) {
        MappingRange cov = e_ovlp_cov_array[i];
        if (cov.end - cov.start > max_size) {
            max_size = cov.end - cov.start;
            max_i = i;
        }
    }
    //HBN_LOG("max_from = %d, max_to = %d", e_ovlp_cov_array[max_i].start, e_ovlp_cov_array[max_i].end);
    hbn_assert(max_size > 0);
    int from = 0, to = 0;
    int i = e_ovlp_cov_array[max_i].start;
    max_size = 0;
    while (i < e_ovlp_cov_array[max_i].end) {
        while (i < e_ovlp_cov_array[max_i].end 
               && 
               cov_stats[i] < opts->min_cov) {
            ++i;
        }
        if (i >= e_ovlp_cov_array[max_i].end) break;
        int j = i + 1;
        while (j < e_ovlp_cov_array[max_i].end
               &&
              cov_stats[j] >= opts->min_cov) {
            ++j;
        }
        if (j - i >= opts->min_size) {
            if (j - i > max_size) {
                max_size = j - i;
                from = i;
                to = j;
            }
        }
        i = j;
    }
    if (max_size < opts->min_size) return;

    kstring_t* cns_subseq = qaln;
    if (!meap_consensus_one_segment(cns_data,
        &from,
        &to,
        subject_length,
        opts->min_size,
        cns_subseq)) return;

//HBN_LOG("from = %d, to = %d", from, to);
    cns_info->cns_from = from;
    cns_info->cns_to = to;
    cns_info->cns_read_size = ks_size(*cns_subseq);
    cns_info->raw_read_size = subject_length;
    kputc('\0', cns_subseq);
    ostringstream os;
    os << ">" << subject_name << ' '
       << "[Seeds:Ovlps:From:To:RawReadLength:CnsReadLength]=["
       << num_extended_can << ':'
       << num_added_aln << ':'
       << from << ':'
       << to << ':'
       << subject_length << ':'
       << cns_info->cns_read_size << "]\n"
       << ks_s(*cns_subseq) << '\n';
    string cns_fasta = os.str();
    cns_info->cns_fasta_size = cns_fasta.size();
    pthread_mutex_lock(data->cns_out_lock);
    cns_info->cns_fasta_offset = ks_size(*data->cns_out);
    kputsn(cns_fasta.c_str(), cns_fasta.size(), data->cns_out);
    pthread_mutex_unlock(data->cns_out_lock);

   // exit(0);
}