#ifndef FSA_FASTQ_READER_HPP
#define FSA_FASTQ_READER_HPP

#include <string>
#include <unordered_map>
#include <fstream>
#include <tuple>

#include "sequence.hpp"

class FastqReader : public SeqReader {
public:
    FastqReader(const std::string &fname);
    virtual ~FastqReader() {}

    bool IsValid() const { return in_.is_open(); }
    bool IsFileEnd() { return Tell() == -1; }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item); }

    
protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq) { return (bool)std::getline(in_, seq); }
    bool GetHead1() { std::string line = NextNonEmptyLine(); return line[0] == '+'; }
    bool GetQuality(std::string &qua) { return (bool)std::getline(in_, qua); }

    std::string NextNonEmptyLine() { return SeqReader::NextNonEmptyLine(in_); }

    ItemId Tell() { return in_.tellg(); }
    void Seek(ItemId id) { in_.clear(); in_.seekg(id, std::ios::beg); }
protected:
    std::ifstream in_;
};

#endif // FSA_FASTQ_READER_HPP 
