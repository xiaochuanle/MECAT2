#include "string2hsp.h"

#include "../ncbi_blast/str_util/ncbistr.hpp"

#include <sstream>

using namespace std;
using namespace NStr;

extern "C"
HbnHSP* string2hsp_names(const char* str, 
            NameToIdMap* query_name2id_map, 
            NameToIdMap* subject_name2id_map,
            HbnHSP* hsp)
{
    vector<string> components;
    CTempString tstr(str);
    CTempString delim("\t");
    Split(tstr, delim, components);
    hbn_assert(components.size() == 12);
    
    hsp->qid = Name2IdMap_name2id(query_name2id_map, components[0].c_str());
    hsp->sid = Name2IdMap_name2id(subject_name2id_map, components[1].c_str());
    hsp->ident_perc = StringToDouble(components[2]);
    hsp->score = StringToInt(components[3]);
    
    hsp->qdir = StringToInt(components[4]);
    hbn_assert(hsp->qdir == FWD || hsp->qdir == REV);
    hsp->qoff = StringToUInt8(components[5]);
    hsp->qend = StringToUInt8(components[6]);
    hsp->qsize = StringToUInt8(components[7]);

    hsp->sdir = StringToInt(components[8]);
    hbn_assert(hsp->sdir == FWD || hsp->sdir == REV);
    hsp->soff = StringToUInt8(components[9]);
    hsp->send = StringToUInt8(components[10]);
    hsp->ssize = StringToUInt8(components[11]);

    hsp->use_normalised_offset = 0;
    hsp->qaln_offset = 0;
    hsp->saln_offset = 0;

    return hsp;
}

extern "C"
const char* hsp2string_names(const HbnHSP* hsp, 
    NameToIdMap* query_name2id_map, 
    NameToIdMap* subject_name2id_map,
    kstring_t* hspstr)
{
    ostringstream os;
    const char delim = '\t';

    os << Name2IdMap_id2name(query_name2id_map, hsp->qid) << delim
       << Name2IdMap_id2name(subject_name2id_map, hsp->sid) << delim
       << hsp->ident_perc << delim
       << hsp->score << delim
       << hsp->qdir << delim
       << hsp->qoff << delim
       << hsp->qend << delim
       << hsp->qsize << delim
       << hsp->sdir << delim
       << hsp->soff << delim
       << hsp->send << delim
       << hsp->ssize;

    string s = os.str();
    ks_clear(*hspstr);
    kputsn(s.c_str(), s.size(), hspstr);
    kputc('\0', hspstr);
    return ks_s(*hspstr);
}

extern "C"
HbnHSP* string2hsp(const char* str, HbnHSP* hsp)
{
    vector<string> components;
    CTempString tstr(str);
    CTempString delim("\t");
    Split(tstr, delim, components);
    hbn_assert(components.size() == 12);
    
    hsp->qid = StringToInt(components[0]);
    hsp->sid = StringToInt(components[1]);
    hsp->ident_perc = StringToDouble(components[2]);
    hsp->score = StringToInt(components[3]);
    
    hsp->qdir = StringToInt(components[4]);
    hbn_assert(hsp->qdir == FWD || hsp->qdir == REV);
    hsp->qoff = StringToUInt8(components[5]);
    hsp->qend = StringToUInt8(components[6]);
    hsp->qsize = StringToUInt8(components[7]);

    hsp->sdir = StringToInt(components[8]);
    hbn_assert(hsp->sdir == FWD || hsp->sdir == REV);
    hsp->soff = StringToUInt8(components[9]);
    hsp->send = StringToUInt8(components[10]);
    hsp->ssize = StringToUInt8(components[11]);

    hsp->use_normalised_offset = 0;
    hsp->qaln_offset = 0;
    hsp->saln_offset = 0;

    return hsp;    
}

extern "C"
const char* hsp2string(const HbnHSP* hsp, kstring_t* hspstr)
{
    ostringstream os;
    const char delim = '\t';

    os << hsp->qid << delim
       << hsp->sid << delim
       << hsp->ident_perc << delim
       << hsp->score << delim
       << hsp->qdir << delim
       << hsp->qoff << delim
       << hsp->qend << delim
       << hsp->qsize << delim
       << hsp->sdir << delim
       << hsp->soff << delim
       << hsp->send << delim
       << hsp->ssize;

    string s = os.str();
    ks_clear(*hspstr);
    kputsn(s.c_str(), s.size(), hspstr);
    kputc('\0', hspstr);
    return ks_s(*hspstr);    
}

extern "C" {

HbnHspReader*
HbnHspReaderNew(const char* hsp_path, 
    const char* query_names,
    const int num_qureis,
    const char* subject_names,
    const int num_subjects)
{
    HBN_LOG("number of quries: %d, number of subjects: %d", num_qureis, num_subjects);
    HbnHspReader* reader = (HbnHspReader*)calloc(1, sizeof(HbnHspReader));
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

HbnHspReader*
HbnHspReaderFree(HbnHspReader* reader)
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
HbnHspReaderGet(HbnHspReader* reader, HbnHSP* hsp)
{
    if (reader->unget) {
        memcpy(hsp, &reader->hsp, sizeof(HbnHSP));
        reader->unget = 0;
        return 1;
    }

    if (HbnLineReaderAtEof(reader->line_reader)) return 0;

    HbnLineReaderReadOneLine(reader->line_reader);
    kstring_t* line = HbnLineReaderLine(reader->line_reader);
    kputc('\0', line);
    if (reader->query_name2id_map) {
        hbn_assert(reader->subject_name2id_map);
        string2hsp_names(ks_s(*line), reader->query_name2id_map, reader->subject_name2id_map, hsp);
    } else {
        string2hsp(ks_s(*line), hsp);
    }
    memcpy(&reader->hsp, hsp, sizeof(HbnHSP));
    return 1;
}

void
HbnHspReaderUnget(HbnHspReader* reader)
{
    reader->unget = 1;
}

}