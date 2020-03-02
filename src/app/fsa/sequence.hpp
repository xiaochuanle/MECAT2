#ifndef FSA_SEQUENCE_HPP
#define FSA_SEQUENCE_HPP

#include <cassert>
#include <string>

#include <fstream>

struct Seq {
    typedef int Id;
    typedef int EndId;
    static EndId IdToEndId(Id id, int end) {
        assert(id >= 0 && (end ==0 || end == 1));
        return end == 0 ? id + 1 : -id - 1;
    }


    static Id EndIdToId(EndId id) {
        assert(id != 0);
        return id > 0 ? id - 1 : (-id) - 1;
    }

    static EndId ReverseEndId(EndId id) { return -id; }

    struct Area {
        Seq::Id id;
        int strand;
        int start;
        int end;
    };
    
    static std::string ReverseComplement(const std::string &seq);
};


class SeqReader {
public:
    using ItemId = int64_t;
    struct Item {
        std::string head;
        std::string sub_head;
        std::string seq;
        std::string quality;
        ItemId id;            // location in file
    };
public:
    virtual bool Next(Item &item) = 0;
    virtual bool Get(ItemId id, Item &item) = 0;

    static std::string NextNonEmptyLine(std::ifstream &in);
};


#endif // FSA_SEQUENCE_HPP
