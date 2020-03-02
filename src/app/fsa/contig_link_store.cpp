#include "contig_link_store.hpp"

#include <array>
#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <iostream>
#include <tuple>

#include "overlap.hpp"
#include "utility.hpp"



void ContigLinkStore::SetParameter(const std::string& name, int v) {
    if (name == "read2ctg_min_identity") read2ctg_min_identity_ = v;
    else if (name == "ctg2ctg_min_identity") ctg2ctg_min_identity_ = v;
    else if (name == "read_min_length") read_min_length_ = v;
    else if (name == "ctg_min_length") ctg_min_length_ = v;
    else if (name == "read2ctg_max_overhang") read2ctg_max_overhang_ = v;
    else if (name == "ctg2ctg_max_overhang") ctg2ctg_max_overhang_ = v;
    else if (name == "read2ctg_min_aligned_length") read2ctg_min_aligned_length_ = v;
    else if (name == "ctg2ctg_min_aligned_length") ctg2ctg_min_aligned_length_ = v;
    else if (name == "read2ctg_min_coverage") read2ctg_min_coverage_ = v;
    else if (name == "thread_size") thread_size_ = v;
    else;
}


void ContigLinkStore::LoadR2cFile(const std::string &fname) {

    auto filter_simple = [&](const Overlap& o)->bool {
        return o.identity_ >= read2ctg_min_identity_ && o.a_.id != o.b_.id &&
               o.a_.len >= read_min_length_ && o.b_.len >= ctg_min_length_ && 
               // o.AlignedLength() >= (size_t)read2ctg_min_aligned_length_/3 &&   // TODO removing condition is for short contigs. 
               o.Location(read2ctg_max_overhang_) != Overlap::Loc::Abnormal;
    };

    read2ctg_.Load(fname, "", thread_size_, filter_simple);

    auto better = [](const Overlap* a, const Overlap *b) { return a->AlignedLength() > b->AlignedLength(); };
    read2ctg_group_ = read2ctg_.GroupTarget(better);

    for (auto &ctg0 : read2ctg_group_) {
        for (auto &ctg1 : read2ctg_group_) {

            if (ctg0.first < ctg1.first) {

                auto keys = FindIntersectKeys(ctg0.second, ctg1.second);

                for (auto k : keys) {
                    auto &o0 = *ctg0.second[k];
                    auto &o1 = *ctg1.second[k];
                    
                    if (ContigLink::SimpleValid(o0, o1, read2ctg_max_overhang_)) {
                        links_[ctg0.first][ctg1.first].Add(o0, o1);
                    }
                }
            }
        }

    }
}


void ContigLinkStore::LoadC2cFile(const std::string &fname) {

    auto filter_simple = [&](const Overlap& o) -> bool {
        return o.identity_ >= ctg2ctg_min_identity_ && o.a_.id != o.b_.id &&
               o.a_.len >= ctg_min_length_ && o.b_.len >= ctg_min_length_ && 
               o.AlignedLength() >= (size_t)ctg2ctg_min_aligned_length_ && 
               o.Location(ctg2ctg_max_overhang_) != Overlap::Loc::Abnormal;
    };

    ctg2ctg_.Load(fname, "", thread_size_, filter_simple);

    for (auto &o : ctg2ctg_.Get()) {
        if (ContigLink::SimpleValid(o, ctg2ctg_max_overhang_)) {
            assert(o.a_.id != o.b_.id);
            if (o.a_.id < o.b_.id) {
                links_[o.a_.id][o.b_.id].Add(o, o.a_, o.b_, ctg2ctg_max_overhang_);
            }
            else {
                links_[o.b_.id][o.a_.id].Add(o, o.b_, o.a_, ctg2ctg_max_overhang_);
            }
        }
    }
}

void ContigLinkStore::AnalyzeSupport() {
    for (auto &i0 : links_) {
        for (auto &i1 : i0.second) {
            assert (i0.first < i1.first);
            i1.second.AnalyzeLinks(read2ctg_max_overhang_, ctg2ctg_max_overhang_, read2ctg_min_coverage_, read2ctg_min_aligned_length_);
        }
    }
}

ContigLink::Loc ContigLinkStore::Location(const ContigLink& bunch) const {
    if (bunch.best_c2c != nullptr) {
        return bunch.best_c2c->Location(ctg2ctg_max_overhang_);
    } else if (bunch.BestC2r2c() != nullptr) {
        return bunch.BestC2r2c()->Location(read2ctg_max_overhang_);
    } else {
        return ContigLink::Loc::Abnormal;
    }
}

void ContigLinkStore::Dump(const std::string &fname) {
    std::ofstream of(fname);
    if (of.is_open()) {
        for (auto &ctg0 : links_) {
            for (auto &ctg1 : ctg0.second) {
                auto& link = ctg1.second;

                of << ctg0.first << " " << read_store_.IdToName(ctg0.first) << " " << ctg1.first << " " << read_store_.IdToName(ctg1.first) << "\n";
                of << "Contig_Contig\n";
                for (auto &i : link.c2c_links) {
                    for (auto& o : i.ols) {
                        of << ctg2ctg_.ToM4aLine(*o);
                    }
                    of << "ol_expect:" << i.ol_expect[0] << " " << i.ol_expect[1] << "\n";
                }
                of << "Read_Contig\n";
                for (auto &i : link.c2r2c_links) {
                    for (auto& o : i.ols) {
                        of << read2ctg_.ToM4aLine(*o);
                    }
                    of << "s2t:" << i.pos_s2t[0] << " " << i.pos_s2t[1] << " " << i.pos_s2t[2] << " " << i.pos_s2t[3] << "\n";
                    of << "ol_expect:" << i.ol_expect[0] << " " << i.ol_expect[1] << "\n";
                }

                if (link.Best() != nullptr) {
                    of << "Best\n";
                    for (auto o : link.Best()->ols) {
                        of << read2ctg_.ToM4aLine(*o);

                    }
                }
            }
        }
    }
    else {
        LOG(ERROR)("Failed to write file %s", fname.c_str());
    }

}

