#include "read_store.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>

#include "fasta_reader.hpp"
#include "fastq_reader.hpp"
#include "logger.hpp"

int ReadStore::NameToId(const std::string &name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = names_to_ids_.find(name);
    if (it != names_to_ids_.end()) {
        return it->second;
    } else {
        
        names_.push_back(name);
        items_.push_back(Item());
        names_to_ids_[name] = (int)names_.size() - 1;
        return (int)names_.size() - 1;
    }
}


    std::string ReadStore::IdToName(int id) {
        assert(0 <= id && (size_t)id < names_.size());
        return names_[id];
    }

    void ReadStore::SetNameToId(const std::string &name, int id) {
        names_to_ids_[name] = id;

    assert(id >= 0);

    if (names_.size() <= (size_t)id) {
        names_.insert(names_.end(), id + 1 - names_.size(), std::string());
        names_.back() = name;
    }
    else {
        names_[id] = name;
    }

}

void ReadStore::Load(const std::string &fname, const std::string &type, int mode) {
    std::string t = type != "" ? type : DetectFileType(fname);
    if (t == "fasta") {
        LoadFasta(fname, mode);
    } else if (t == "fastq") {
        LoadFastq(fname, mode);
    } else if (t == "fofn") {
        LoadFofn(fname, mode);
    } else if (t == "txt") {
        LoadTxt(fname, mode);
    } else {
        LOG(ERROR)("Failed to recognize read files type: %s", t.c_str());
    }
}


void ReadStore::LoadFasta(const std::string &fname, int mode) {
    std::unordered_set<Seq::Id> ids;

    FastaReader* reader = new FastaReader(fname);
    readers_.push_back(reader);

    if (reader->IsValid()) {
        SeqReader::Item item;
        while (reader->Next(item)) {
            assert(!item.head.empty());
            ids.insert(Insert(item, reader, mode));
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(FATAL)("Failed to open file: %s", fname.c_str());
    }
    ids_in_file_[fname] = ids;
    LOG(INFO)("Load %zd reads from fasta file: %s", ids.size(), fname.c_str());
}


void ReadStore::LoadFastq(const std::string &fname, int mode) {
    std::unordered_set<Seq::Id> ids;

    FastqReader* reader = new FastqReader(fname);
    readers_.push_back(reader);
    
    if (reader->IsValid()) {
        SeqReader::Item item;

        while (reader->Next(item)) {
            assert(!item.head.empty());
            ids.insert(Insert(item, reader, mode));
        }
        if (!reader->IsFileEnd()) {
            LOG(WARNING)("No all reads in file are loaded: %s", fname.c_str());
        }
    } else {
        LOG(FATAL)("Failed to open file: %s", fname.c_str());
    }
    ids_in_file_[fname] = ids;
    LOG(INFO)("Load %zd reads from fastq file: %s", ids.size(), fname.c_str());
}

void ReadStore::LoadFofn(const std::string &fname, int mode) {
    std::ifstream in(fname);
    if (in.is_open()) {
        std::string line;
        while (std::getline(in, line)) {
            auto begin = std::find_if(line.begin(), line.end(), [](char a){return !::isspace(a); });
            if (begin != line.end()) {
                Load(line, "", mode);
            }
        }
    } else {
        LOG(FATAL)("Failed to open file: %s", fname.c_str());
    }
}

const std::unordered_set<Seq::Id>& ReadStore::IdsInFile(const std::string &fname) const {
    auto iter = ids_in_file_.find(fname);
    assert(iter != ids_in_file_.end());
    return iter->second;
}



std::string ReadStore::DetectFileType(const std::string &fname) {

    if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".fasta") {
        return "fasta";
    } else if (fname.size() >= 3 && fname.substr(fname.size()-3) == ".fa") {
        return "fasta";
    } else if (fname.size() >= 6 && fname.substr(fname.size()-6) == ".fastq") {
        return "fastq";
    } else if (fname.size() >= 5 && fname.substr(fname.size()-5) == ".fofn") {
        return "fofn";
    } else if (fname.size() >= 4 && fname.substr(fname.size()-4) == ".txt") {
        return "txt";
    } else {
        return "fasta";
    }
}

/*
Seq::Id ReadStore::Insert(std::string &&name, std::string &&seq) {
    Seq::Id id;
    auto it = names_to_ids_.find(name);
    if (it != names_to_ids_.end()) {
        sequences_[it->second] = std::move(seq);
        id = it->second;
    }
    else {
        names_.push_back(std::move(name));
        sequences_.push_back(std::move(seq));
        names_to_ids_[name] = (int)names_.size() - 1;
        id = (int)names_.size() - 1;
    }
    return id;
}
*/

Seq::Id ReadStore::Insert(const SeqReader::Item &item, SeqReader *reader, int mode) {
    Seq::Id id;
    auto it = names_to_ids_.find(item.head);
    
    if (it != names_to_ids_.end()) {
        items_[it->second].seq = mode == 0 ? item.seq : "";
        items_[it->second].id = item.id;
        items_[it->second].reader = reader;
        id = it->second;
    } else {
        names_.push_back(item.head);
        items_.push_back(Item(mode == 0 ? item.seq : "", item.id, reader));
        names_to_ids_[item.head] = (int)names_.size() - 1;
        id = (int)names_.size() - 1;
    }
    
    return id;

}

void ReadStore::SaveIdToName(const std::string &fname) const {
    std::ofstream out(fname);

    if (out.is_open()) {
        for (size_t i = 0; i < names_.size(); i++) {
            out << i << " " << names_[i] << "\n";
        }
    }
    else {
        LOG(ERROR)("Failed to open outfile: %s", fname.c_str());
    }
}


std::string ReadStore::GetSeq(const Seq::Area& sa) {
    const std::string &seq = GetSeq(sa.id);

    if (sa.strand == 0) {
        return seq.substr(sa.start, sa.end);
    }
    else {
        return Seq::ReverseComplement(seq.substr(sa.start, sa.end));
    }
}