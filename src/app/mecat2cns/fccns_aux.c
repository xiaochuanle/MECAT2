#include "fccns_aux.h"


static inline u8
encode_dna_base(const char c) 
{
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        case '-': return 4;
    }
    HBN_ERR("invalid dna base: %c", c);
    return 4;
}

static void
build_base_links(AlignTag* tag_array, const int tag_count, BaseLinks* link, SmallObjectAlloc* alloc)
{
    int n_link = 0;
    int i = 0;
    while (i < tag_count) {
        int j = i + 1;
        while (j < tag_count && align_tag_plink_eq(tag_array[i], tag_array[j])) ++j;
        ++n_link;
        i = j;
    }

    link->plinks = (LinkInfo*)SmallObjectAllocAlloc(alloc, n_link);
    link->n_link = n_link;
    link->coverage = tag_count;

    n_link = 0;
    i = 0;
    while (i < tag_count) {
        int j = i + 1;
        while (j < tag_count && align_tag_plink_eq(tag_array[i], tag_array[j])) ++j;
        LinkInfo* lnk_info = link->plinks + n_link;
        lnk_info->p_t_pos = tag_array[i].p_t_pos;
        lnk_info->p_delta = tag_array[i].p_delta;
        lnk_info->p_q_base = tag_array[i].p_q_base;
        lnk_info->link_count = j - i;
        lnk_info->weight = .0;
        for (int k = i; k < j; ++k) lnk_info->weight += tag_array[k].weight;

        ++n_link;
        i = j;
    }
    hbn_assert(n_link == link->n_link);
}

static void
build_delta_links(AlignTag* tag_array, const int tag_count, DeltaCovInfo* dci, SmallObjectAlloc* alloc)
{
    for (int i = 0; i < 5; ++i) init_base_links(dci->links[i]);
    int i = 0;
    while (i < tag_count) {
        int j = i;
        while (j < tag_count && tag_array[i].q_base == tag_array[j].q_base) ++j;
        const u8 c = encode_dna_base(tag_array[i].q_base);
        build_base_links(tag_array + i, j - i, dci->links + c, alloc);
        i = j;
    }
}

static void
build_backbaone_item(AlignTag* tag_array,
    const int tag_count,
    BackboneItem* item,
    SmallObjectAlloc* dci_alloc,
    SmallObjectAlloc* li_alloc,
    int* cov_array)
{
    item->n_delta = tag_array[tag_count-1].delta + 1;
    item->delta = (DeltaCovInfo*)SmallObjectAllocAlloc(dci_alloc, item->n_delta);
    int i = 0;
    while (i < tag_count) {
        int j = i + 1;
        while (j < tag_count && tag_array[i].delta == tag_array[j].delta) ++j;
        build_delta_links(tag_array + i, j - i, item->delta + tag_array[i].delta, li_alloc);
        if (tag_array[i].delta == 0) cov_array[ tag_array[i].t_pos ] = j - i;
        i = j;
    }
}

void
build_backbone(AlignTag* tag_array,
    const int tag_count,
    const int template_size,
    SmallObjectAlloc* dci_alloc,
    SmallObjectAlloc* li_alloc,
    vec_backbone_item* item_list,
    vec_int* cov_list)
{
    kv_resize(BackboneItem, *item_list, template_size);
    BackboneItem* backbone = kv_data(*item_list);
    for (int i = 0; i < template_size; ++i) clear_backbone_item(backbone[i]);
    kv_resize(int, *cov_list, template_size);
    kv_zero(int, *cov_list);
    int* cov_array = kv_data(*cov_list);

    ks_introsort_align_tag_lt(tag_count, tag_array);
    int i = 0;
    while (i < tag_count) {
        int j = i + 1;
        while (j < tag_count && tag_array[i].t_pos == tag_array[j].t_pos) ++j;
        hbn_assert(tag_array[i].t_pos < template_size);
        build_backbaone_item(tag_array + i, j - i, backbone + tag_array[i].t_pos, dci_alloc, li_alloc, cov_array);
        i = j;
    }
}

static void
reverse_kstring(kstring_t* str)
{
    char* p = ks_s(*str);
    int s = 0, e = ks_size(*str);
    while (s < e) {
        char tmp = p[s];
        p[s] = p[e];
        p[e] = tmp;
        ++s;
        --e;
    }
}

void
consensus_backbone_segment(BackboneItem* backbone,
        int from,
        int to,
        int* coverage,
        kstring_t* cns_seq,
        int* cns_from,
        int* cns_to)
{
    int g_best_ck = 0;
    BaseLinks* g_best_aln_col = 0;
    int g_best_t_pos = 0;
    double g_best_score = -1.0;
    int g_best_q_base = -1;
    int kk;
    int ck;
    int best_i;
    int best_j;
    int best_b;
    int best_ck = -1;
    double score;
    double best_score;
    BaseLinks* aln_col;
    const double INDEL_FACTOR = 0.1;

    for (int i = from; i < to; ++i) {
        for (int j = 0; j < backbone[i].n_delta; ++j) {
            for (kk = 0; kk < 5; ++kk) {
                aln_col = backbone[i].delta[j].links + kk;
                if (aln_col->coverage) {
                    best_score = -1;
                    best_i = -1;
                    best_j = -1;
                    best_b = -1;

                    for (ck = 0; ck < aln_col->n_link; ++ck) {
                        int pi = aln_col->plinks[ck].p_t_pos;
                        int pj = aln_col->plinks[ck].p_delta;
                        int pkk = encode_dna_base(aln_col->plinks[ck].p_q_base);
                        score = aln_col->plinks[ck].weight - INDEL_FACTOR * coverage[i];
                        if (pi != -1) {
                            score += backbone[pi].delta[pj].links[pkk].score;
                        }

                        if (score > best_score) {
                            best_score = score;
                            aln_col->best_p_t_pos = best_i = pi;
                            aln_col->best_p_delta = best_j = pj;
                            aln_col->best_p_q_base = best_b = pkk;
                            best_ck = ck;
                        }
                    }

                    aln_col->score = best_score;
                    if (best_score > g_best_score) {
                        g_best_score = best_score;
                        g_best_aln_col = aln_col;
                        g_best_ck = best_ck;
                        g_best_t_pos = i;
                        g_best_q_base = kk;
                    }
                }
            }
        }
    }
    hbn_assert(g_best_score != -1, "g_best_score = %d", g_best_score);

    char bb = '$';
    ck = g_best_q_base;
    int i = g_best_t_pos;
    int _cns_to = i + 1, _cns_from = 0;
    int j;
    ks_clear(*cns_seq);
    while (1) {
        bb = ck;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = backbone[i].delta[j].links + ck;
        _cns_from = i;

        if (bb != 4) {
            hbn_assert(bb >= 0 && bb < 4, "i = %d, bb = %d", i, bb);
            kputc(bb, cns_seq);
        }
    }
	
	reverse_kstring(cns_seq);
	if (cns_from) *cns_from = _cns_from;
	if (cns_to) *cns_to = _cns_to;
}
