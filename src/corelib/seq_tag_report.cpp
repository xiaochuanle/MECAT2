#include "seq_tag_report.h"

#include "../ncbi_blast/str_util/ncbistr.hpp"
#include "seq_tag.h"

#include <sstream>
#include <vector>

using namespace ncbi;
using namespace std;

static const size_t kItemWidth = 16;

static void
print_one_tag_report(const char* tagname, int seq, size_t res, double ratio, ostream& os)
{
    string item = tagname;
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';
    
    item = NStr::IntToString(seq, NStr::fWithCommas, 10);
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    item = NStr::UInt8ToString_DataSize(res);
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    item = NStr::DoubleToString(ratio);
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    os << endl;
}

static void
print_report_label(ostream& os)
{
    string item = "tag";
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    item = "sequences";
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    item = "residues";
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';

    item = "percentage";
    os << item;
    for (size_t i = item.size(); i < kItemWidth; ++i) os << ' ';   

    os << endl; 
}

extern "C"
void 
print_tag_report(int seq_cnt_array[], size_t res_cnt_array[], kstring_t* out)
{
    ks_clear(*out);

    int total_seq = 0;
    size_t total_res = 0;
    for (int i = 0; i < SEQ_TAG_MAX; ++i) {
        total_seq += seq_cnt_array[i];
        total_res += res_cnt_array[i];
    }
    if (total_seq == 0) return;

    ostringstream os;
    print_report_label(os);

    print_one_tag_report("all", total_seq, total_res, 100.0, os);
    for (int i = 0; i < SEQ_TAG_MAX; ++i) {
        if (seq_cnt_array[i] == 0) continue;
        double percentage = 100.0 * res_cnt_array[i] / total_res;
        print_one_tag_report(GetSeqTagName(i), seq_cnt_array[i], res_cnt_array[i], percentage, os);
    }

    string report = os.str();
    kputsn(report.c_str(), report.size(), out);
    kputc('\0', out);
}

extern "C"
const char*
seqtag2string(SeqTag* tag, 
    NameToIdMap* query_name2id_map,
    NameToIdMap* subject_name2id_map,
    kstring_t* tagstr)
{
    ostringstream os;
    const char kDelim = '\t';

    os << tag->qid;
    if (query_name2id_map) os << ':' << Name2IdMap_id2name(query_name2id_map, tag->qid);
    os << kDelim;

    os << tag->sid;
    if (subject_name2id_map && tag->sid != I32_MAX) os << ':' << Name2IdMap_id2name(subject_name2id_map, tag->sid);
    os << kDelim;

    os << GetSeqTagName(tag->tag) << kDelim;

    os << tag->ident_perc << kDelim;

    os << tag->qoff << kDelim
       << tag->qend << kDelim
       << tag->qsize << kDelim
       << tag->soff << kDelim
       << tag->send;
    
    string s = os.str();
    ks_clear(*tagstr);
    kputsn(s.c_str(), s.size(), tagstr);
    kputc('\0', tagstr);
    return ks_s(*tagstr);
}

extern "C"
void string2seqtag(const char* str, SeqTag* tag)
{
    vector<string> components;
    CTempString kDelim("\t");
    NStr::Split(str, kDelim, components);
    hbn_assert(components.size() == 9);
    
    tag->qid = NStr::StringToInt(components[0], NStr::fAllowTrailingSymbols, 10);
    tag->sid = NStr::StringToInt(components[1], NStr::fAllowTrailingSymbols, 10);
    tag->tag = SeqNameToTag(components[2].c_str());
    tag->ident_perc = NStr::StringToDouble(components[3].c_str());
    tag->qoff = NStr::StringToUInt8(components[4]);
    tag->qend = NStr::StringToUInt8(components[5]);
    tag->qsize = NStr::StringToUInt8(components[6]);
    tag->soff = NStr::StringToUInt8(components[7]);
    tag->send = NStr::StringToUInt8(components[8]);
}

extern "C"
void dump_seq_tag(SeqTag* tag, FILE* out)
{
    ks_dinit(tagstr);
    seqtag2string(tag, NULL, NULL, &tagstr);
    fprintf(out, "%s\n", ks_s(tagstr));
    ks_destroy(tagstr);
}