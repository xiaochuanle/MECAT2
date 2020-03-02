#include "hsp2string.h"

#include "../str_util/ncbistr.hpp"

#include <sstream>
#include <stdarg.h>
#include <string>

using namespace std;
using namespace NStr;

extern "C"
const char* blasthsp_to_string_names(const BlastHSP* hsp, ...)
{
    va_list arg_ptr;
    va_start(arg_ptr, hsp);
    kstring_t* hspstr = va_arg(arg_ptr, kstring_t*);
    const char* query_name = va_arg(arg_ptr, const char*);
    const char* subject_name = va_arg(arg_ptr, const char*);
    va_end(arg_ptr);

    ostringstream os;
    const char kDelim = '\t';
    os << query_name 
       << kDelim
       << subject_name
       << kDelim
       << hsp->hsp_info.perc_identity
       << kDelim
       << hsp->hsp_info.raw_score
       << kDelim
       << hsp->hbn_query.strand
       << kDelim
       << hsp->hbn_query.offset
       << kDelim
       << hsp->hbn_query.end
       << kDelim
       << hsp->hbn_query.seq_size
       << kDelim
       << hsp->hbn_subject.strand
       << kDelim
       << hsp->hbn_subject.offset
       << kDelim
       << hsp->hbn_subject.end
       << kDelim
       << hsp->hbn_subject.seq_size
       << endl;
    string s = os.str();
    ks_clear(*hspstr);
    kputsn(s.c_str(), s.size(), hspstr);
    //kputc('\0', hspstr);
    return ks_s(*hspstr);    
}

extern "C"
const char* blasthsp_to_string_name2idmap(const BlastHSP* hsp, ...)
{
    va_list arg_ptr;
    va_start(arg_ptr, hsp);
    kstring_t* hspstr = va_arg(arg_ptr, kstring_t*);
    NameToIdMap* name2id_map = va_arg(arg_ptr, NameToIdMap*);
    va_end(arg_ptr);

    const char* query_name = Name2IdMap_id2name(name2id_map, hsp->hbn_query.oid);
    const char* subject_name = Name2IdMap_id2name(name2id_map, hsp->hbn_subject.oid);
    return blasthsp_to_string_names(hsp, hspstr, query_name, subject_name);
}

extern "C"
const char* blasthsp_to_string_ids(const BlastHSP* hsp, ...)
{
    va_list arg_ptr;
    va_start(arg_ptr, hsp);
    kstring_t* hspstr = va_arg(arg_ptr, kstring_t*);
    va_end(arg_ptr);

    string query_name = NStr::IntToString(hsp->hbn_query.oid);
    string subject_name = NStr::IntToString(hsp->hbn_subject.oid);
    return blasthsp_to_string_names(hsp, hspstr, query_name.c_str(), subject_name.c_str());
}

////////////

extern "C"
BlastHSP* string_to_blasthsp_ids(const char* str, ...)
{
    va_list arg_ptr;
    va_start(arg_ptr, str);
    BlastHSP* hsp = va_arg(arg_ptr, BlastHSP*);
    va_end(arg_ptr);

    vector<string> components;
    CTempString tstr(str);
    CTempString delim("\t");
    Split(tstr, delim, components);
    hbn_assert(components.size() == 12);
    
    hsp->hbn_query.oid = StringToInt(components[0]);
    hsp->hbn_subject.oid = StringToInt(components[1]);
    hsp->hsp_info.perc_identity = StringToDouble(components[2]);
    hsp->hsp_info.raw_score = StringToInt(components[3]);
    
    hsp->hbn_query.strand = StringToInt(components[4]);
    hbn_assert(hsp->hbn_query.strand == FWD || hsp->hbn_query.strand == REV);
    hsp->hbn_query.offset = StringToUInt8(components[5]);
    hsp->hbn_query.end = StringToUInt8(components[6]);
    hsp->hbn_query.seq_size = StringToUInt8(components[7]);

    hsp->hbn_subject.strand = StringToInt(components[8]);
    hbn_assert(hsp->hbn_subject.strand == FWD || hsp->hbn_subject.strand == REV);
    hsp->hbn_subject.offset = StringToUInt8(components[9]);
    hsp->hbn_subject.end = StringToUInt8(components[10]);
    hsp->hbn_subject.seq_size = StringToUInt8(components[11]);

    hsp->hsp_info.query_align_offset = 0;
    hsp->hsp_info.subject_align_offset = 0;

    return hsp;    
}

extern "C"
BlastHSP* string_to_blasthsp_name2idmap(const char* str, ...)
{
    va_list arg_ptr;
    va_start(arg_ptr, str);
    BlastHSP* hsp = va_arg(arg_ptr, BlastHSP*);
    NameToIdMap* query_name2id_map = va_arg(arg_ptr, NameToIdMap*);
    NameToIdMap* subject_name2id_map = va_arg(arg_ptr, NameToIdMap*);
    va_end(arg_ptr);

    vector<string> components;
    CTempString tstr(str);
    CTempString delim("\t");
    Split(tstr, delim, components);
    hbn_assert(components.size() == 12);
    
    hsp->hbn_query.oid = Name2IdMap_name2id(query_name2id_map, components[0].c_str());
    hsp->hbn_subject.oid = Name2IdMap_name2id(subject_name2id_map, components[1].c_str());
    hsp->hsp_info.perc_identity = StringToDouble(components[2]);
    hsp->hsp_info.raw_score = StringToInt(components[3]);
    
    hsp->hbn_query.strand = StringToInt(components[4]);
    hbn_assert(hsp->hbn_query.strand == FWD || hsp->hbn_query.strand == REV);
    hsp->hbn_query.offset = StringToUInt8(components[5]);
    hsp->hbn_query.end = StringToUInt8(components[6]);
    hsp->hbn_query.seq_size = StringToUInt8(components[7]);

    hsp->hbn_subject.strand = StringToInt(components[8]);
    hbn_assert(hsp->hbn_subject.strand == FWD || hsp->hbn_subject.strand == REV);
    hsp->hbn_subject.offset = StringToUInt8(components[9]);
    hsp->hbn_subject.end = StringToUInt8(components[10]);
    hsp->hbn_subject.seq_size = StringToUInt8(components[11]);

    hsp->hsp_info.query_align_offset = 0;
    hsp->hsp_info.subject_align_offset = 0;

    return hsp;
}

BlastHSPReader*
BlastHSPReaderNew(const char* hsp_path, 
    const char* query_names,
    const int num_qureis,
    const char* subject_names,
    const int num_subjects)
{
    HBN_LOG("number of quries: %d, number of subjects: %d", num_qureis, num_subjects);
    BlastHSPReader* reader = (BlastHSPReader*)calloc(1, sizeof(BlastHSPReader));
    reader->unget = 0;
    reader->line_reader = HbnLineReaderNew(hsp_path);
    if (query_names) {
        hbn_assert(num_qureis > 0);
        hbn_assert(subject_names != NULL);
        hbn_assert(num_subjects > 0);
        HBN_LOG("build query map");
        reader->query_name2id_map = NameToIdMapNew();
        NameToIdMapSet(query_names, num_qureis, reader->query_name2id_map);
        if (query_names == subject_names) {
            reader->subject_name2id_map = reader->query_name2id_map;
        } else {
            HBN_LOG("build subject map");
            reader->subject_name2id_map = NameToIdMapNew();
            NameToIdMapSet(subject_names, num_subjects, reader->subject_name2id_map);
        }
    } else {
        hbn_assert(num_qureis == 0);
        hbn_assert(subject_names == NULL);
        hbn_assert(num_subjects == 0);
    }

    return reader;
}

BlastHSPReader*
BlastHSPReaderFree(BlastHSPReader* reader)
{
    if (reader->line_reader) reader->line_reader = HbnLineReaderFree(reader->line_reader);
    int free_subject_map = reader->query_name2id_map != reader->subject_name2id_map;
    if (reader->query_name2id_map) NameToIdMapFree(reader->query_name2id_map);
    reader->query_name2id_map = NULL;
    if (free_subject_map && reader->subject_name2id_map) NameToIdMapFree(reader->subject_name2id_map);
    reader->subject_name2id_map = NULL;
    free(reader);
    return NULL;
}

int
BlastHSPReaderGet(BlastHSPReader* reader, BlastHSP* hsp)
{
    if (reader->unget) {
        memcpy(hsp, &reader->hsp, sizeof(BlastHSP));
        reader->unget = 0;
        return 1;
    }

    if (HbnLineReaderAtEof(reader->line_reader)) return 0;

    HbnLineReaderReadOneLine(reader->line_reader);
    kstring_t* line = HbnLineReaderLine(reader->line_reader);
    kputc('\0', line);
    if (reader->query_name2id_map) {
        hbn_assert(reader->subject_name2id_map);
        string2hsp(ks_s(*line), hsp, reader->query_name2id_map, reader->subject_name2id_map);
    } else {
        string2hsp(ks_s(*line), hsp);
    }
    memcpy(&reader->hsp, hsp, sizeof(BlastHSP));
    return 1;
}

void
BlastHSPReaderUnget(BlastHSPReader* reader)
{
    reader->unget = 1;
}