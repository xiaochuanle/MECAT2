#ifndef __SEQ_TAG_REPORT_H
#define __SEQ_TAG_REPORT_H

#include "kstring.h"
#include "name2id_map.h"
#include "seq_tag.h"

#ifdef __cplusplus
extern "C" {
#endif

void 
print_tag_report(int seq_cnt_array[], size_t res_cnt_array[], kstring_t* out);

const char*
seqtag2string(SeqTag* tag, 
    NameToIdMap* query_name2id_map,
    NameToIdMap* subject_name2id_map,
    kstring_t* tagstr);

void string2seqtag(const char* str, SeqTag* tag);

void dump_seq_tag(SeqTag* tag, FILE* out);

#ifdef __cplusplus
}
#endif

#endif // __SEQ_TAG_REPORT_H