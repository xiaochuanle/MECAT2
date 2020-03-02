#include "fccns_align_tag.h"

#include "../../corelib/ksort.h"

#define align_tag_lt(a, b) ( \
	((a).t_pos < (b).t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta < (b).delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base < (b).q_base)\
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos < (b).p_t_pos) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta < (b).p_delta) \
	|| \
	((a).t_pos == (b).t_pos && (a).delta == (b).delta && (a).q_base == (b).q_base && (a).p_t_pos == (b).p_t_pos && (a).p_delta == (b).p_delta && (a).p_q_base < (b).p_q_base) \
	)

KSORT_INIT(align_tag_lt, AlignTag, align_tag_lt);

void
make_align_tags_from_ovlp(const char* qaln,
    const char* taln,
    const size_t aln_size,
    const int qoff,
    const int qend,
    const int toff,
    const int tend,
    const double weight,
    vec_align_tag* align_tag_list)
{
    AlignTag tag;
    tag.weight = weight;
    int jj = 0;
    int i = qoff - 1;
    int j = toff - 1;
    int p_j = -1;
    int p_jj = 0;
    char p_q_base = GAP_CHAR;

    for (size_t p = 0; p < aln_size; ++p) {
        if (qaln[p] != GAP_CHAR) {
            ++i;
            ++jj;
        }
        if (taln[p] != GAP_CHAR) {
            ++j;
            jj = 0;
        }
        hbn_assert(i >= qoff);
        hbn_assert(i < qend);
        hbn_assert(j >= toff);
        hbn_assert(j < tend);

        if (jj >= ALIGN_TAG_MAX_DELTA || p_jj >= ALIGN_TAG_MAX_DELTA) continue;

        tag.t_pos = j;
        tag.p_t_pos = p_j;
        tag.delta = jj;
        tag.p_delta = p_jj;
        tag.q_base = qaln[p];
        tag.p_q_base = p_q_base;

        p_j = j;
        p_jj = jj;
        p_q_base = qaln[p];

        kv_push(AlignTag, *align_tag_list, tag);
    }
}