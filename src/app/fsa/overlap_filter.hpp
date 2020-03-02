#ifndef FSA_OVERLAP_FILTER_HPP
#define FSA_OVERLAP_FILTER_HPP


#include<climits>
#include <sstream>
#include <array>

#include "overlap_store.hpp"
#include "argument_parser.hpp"

class OverlapFilter {
public:
    OverlapFilter();
    virtual ~OverlapFilter();

    bool ParseArgument(int argc, char *const argv[]);

    void Run();
    void Usage();

protected:
    ArgumentParser GetArgumentParser();

protected:
    struct RdReason {
        enum Type {
            RS_OK = 0,
            RS_CONTAINED,
            RS_COVERAGE,
            RS_NO_LONGEST,
            RS_UNKNOWN
        };
        static RdReason Contained(int id) { return RdReason(Type::RS_CONTAINED, id); }
        static RdReason Coverage(int low, int high) { return RdReason(Type::RS_COVERAGE, low, high); }
        static RdReason Coverage(const std::array<int, 2>& a) { return RdReason(Type::RS_COVERAGE, a[0], a[1]); }
        static RdReason NoLongest() { return RdReason(Type::RS_NO_LONGEST); }
        Type type;
        std::array<int, 2> sub;

protected:
        RdReason(Type t=RS_UNKNOWN, int s0=0, int s1=0) { type = t; sub[0] = s0;  sub[1] = s1;}
    };

    struct OlReason {
        enum Type {
            RS_OK = 0,
            RS_SIMPLE,
            RS_DUPLICATE,
            RS_BESTN,
            RS_FILTERED_READ,
            RS_LACK_OF_SUPPORT,
            RS_LOCAL,
            RS_UNKNOWN
        };
        static OlReason Simple() { return OlReason(RS_SIMPLE); }
        static OlReason Duplicate() { return OlReason(RS_DUPLICATE); }
        static OlReason BestN() { return OlReason(RS_BESTN); }
        static OlReason FilteredRead(int id) { return OlReason(RS_FILTERED_READ, id); }
        static OlReason LackOfSupport() { return OlReason(RS_LACK_OF_SUPPORT); }
        static OlReason Local(int t, int v) { return OlReason(RS_LOCAL, t*10000+v); }
        
        OlReason(Type t=RS_UNKNOWN, int p0=0, int p1=0) : type(t), sub{p0, p1} { }

        Type type;
        std::array<int, 2> sub;
    };

    struct ReadStatInfo {
        double identity{0.0};
        int overhang {0}; 
        int score{0};
        int len {-1};
        int aligned { -1 };
        int count {0};  
        int oh_count {0};
        double all_score {0};
        double all_identity{0.0};
    };

    void LoadOverlaps(const std::string &fname);
    void SaveOverlaps(const std::string &fname);

    void FilterSimple();
    void FilterSimpleMt();
    void GroupAndFilterDuplicate();
    void GroupAndFilterDuplicateMt();
    void FilterContained();
    void FilterContainedMt();
    void FilterCoverage();
    void FilterCoverageMt();
    void FilterLackOfSupport();
    void FilterLackOfSupportMt();
    void FilterLocal();
    void FilterLocalMt();
    void FilterLocalStep(Seq::Id id, const std::unordered_map<Seq::Id, const Overlap*>& g);
    void FilterLongest();
    void FilterBestN();
    void FilterBestNMt();

    void AutoSelectParams();
    void AutoSelectMinLength(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectMinIdentity(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectMaxOverhang(const std::unordered_map<Seq::Id, ReadStatInfo> &info);
    void AutoSelectMinAlignedLength(const std::unordered_map<Seq::Id, ReadStatInfo> &info);

    void CheckSimple(const Overlap& o);

    bool CheckEnd(const Overlap &o) const;
    void ModifyEnd(const Overlap &o);
    bool CheckIdentity(Overlap& o) {
        return o.identity_ >= min_identity_;
    }


    bool CheckLength(Overlap &o) {
        return o.a_.len >= min_length_ && o.b_.len >= min_length_ && 
            o.a_.len <= max_length_ && o.b_.len <= max_length_;
    }

    bool CheckAlignedLength(Overlap &o) {
        return o.AlignedLength() >= (size_t)min_aligned_length_;
    }

    std::pair<int, int> CalcMinMaxCoverage(int id, const std::unordered_map<int, const Overlap*>& group);

    std::string OutputPath(const std::string &fname) const { return output_directory_+"/"+fname; }
    void PrintArguments();

    /** Calc coverage param: min_coverage, max_coverage and max_diff_coverage */
    std::array<int, 3> CoverageParam() const;
    std::array<int, 3> CoverageParam1() const;

    size_t FindLongestXHeap(std::vector<std::array<int,2>> &lengths, long long goal);
    size_t FindLongestXSort(std::vector<std::array<int,2>> &lengths, long long goal);   
    size_t FindLongestXHeap(std::vector<int> &lengths, long long goal); 

    void FilterCoverage(int min_coverage, int max_coverage, int max_diff_coverge);

    bool HasSupport(const Overlap &o, int count) const;
    bool HasAlignment(int a, int b, int end, int count, bool exceeding) const;
    std::unordered_set<const Overlap*> FindBestN(const std::pair<int, std::unordered_map<int, const Overlap*>> &groud) const;
  
    bool IsContained(const Overlap& o, std::array<int, 2> &rel);
    bool IsReserved(const Overlap &o) const {
        return GetOlReason(o).type == OlReason::RS_OK;
    }
    void UpdateFilteredRead(const std::unordered_map<Seq::Id, RdReason> &ignored);


    /** Record internal state and variables */
    void Dump() const;      
    void DumpCoverage(const std::string &fname) const;
    void DumpFilteredReads(const std::string &fname) const;
    void DumpFilteredOverlaps(const std::string &fname) const;

    static bool BetterAlignedLength(const Overlap &o0, const Overlap &o1) { return o0.AlignedLength() > o1.AlignedLength(); }
    static void SetOlReason(const Overlap &o, OlReason rs);
    static OlReason GetOlReason(const Overlap &o);

public:
    static bool ParamToGenomeSize(const std::string& str, long long *v);
    static int Percentile(const std::vector<int>& data, double percent);
    static int FirstTrough(const std::vector<int>& data, size_t last, size_t k);
protected:
    double min_identity_{ 90 };         //!< 
    int min_length_{ 2500 };            //!< 
    int max_length_{ INT_MAX };         //!< 
    int min_aligned_length_{ 2500 };    //!< 
    int max_overhang_{ 10 };               
    int min_coverage_{ -1 };             //!< 
    int max_coverage_{ -1 };           //!< 
    int max_diff_coverage_{ -1 };      //!< 
    double coverage_discard_ { 0.01 };
    int bestn_{ 10 };                   //!< 
    std::string overlap_file_type_{ "" };
    int thread_size_{ 4 };           //!< 
    std::string output_directory_ {"."};
    std::string coverage_fname_ { "coverage.txt" };  //!< variable this->coverages_

    double min_identity_median_ {-1};
    double max_overhang_median_ {-1};

    std::string filtered_read_fname { "filtered_reads.txt"};
    std::string filtered_overlap_fname {"filtered_overlaps.txt"};
    std::string ifname_;                       
    std::string ofname_;                        
    long long genome_size_ {0};
    int coverage_ {0};
    std::array<int, 3> coverage_params_;

    OverlapStore ol_store_;
    std::unordered_map<int, std::unordered_map<int, const Overlap*>> groups_;

    std::unordered_map<Seq::Id, std::array<int, 2>> coverages_;                         //!< record min and max base coverages of the reads
    std::unordered_map<Seq::Id, RdReason> filtered_reads_;                              //!< record filtered reads and reason for filtering
};

#endif // FSA_OVERLAP_FILTER_HPP
