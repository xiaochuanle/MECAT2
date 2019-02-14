#ifndef FSA_CONTIG_LINK_STORE_HPP
#define FSA_CONTIG_LINK_STORE_HPP

#include <array>
#include <cassert>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "overlap_store.hpp"
#include "sequence.hpp"
#include "read_store.hpp"
#include "overlap.hpp"
#include "contig_link.hpp"

class ContigLinkStore {
public:


    ContigLinkStore(ReadStore &rs):read_store_(rs), ctg2ctg_(rs), read2ctg_(rs) {}

    void SetParameter(const std::string& name, int v); 

    void LoadR2cFile(const std::string &fname);
    void LoadC2cFile(const std::string &fname);
    void AnalyzeSupport();
    ContigLink::Loc Location(const ContigLink& link) const;
    std::unordered_map<int, std::unordered_map<int, ContigLink>>& Get() { return links_; }
  
    void Dump(const std::string &fname);
protected:

    int read2ctg_min_identity_{ 80 };
    int ctg2ctg_min_identity_{ 90 };

    int read2ctg_max_overhang_{ 100 };
    int ctg2ctg_max_overhang_{ 100 };

    int read_min_length_{ 10000 };
    int ctg_min_length_{ 50000 };

    int read2ctg_min_aligned_length_{ 5000 };
    int ctg2ctg_min_aligned_length_{ 5000 };

    int read2ctg_min_coverage_ { 3 };

    int thread_size_ {1};

    std::unordered_map<int, std::unordered_map<int, Overlap*>> read2ctg_group_;
    std::unordered_map<int, std::unordered_map<int, ContigLink>> links_;

    ReadStore &read_store_;
    OverlapStore ctg2ctg_;
    OverlapStore read2ctg_;
    
};
#endif // FSA_CONTIG_LINK_STORE_HPP  
