#ifndef FSA_FASTA_READER_HPP
#define FSA_FASTA_READER_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>

#include "sequence.hpp"

class FastaReader : public SeqReader {
public:
    FastaReader(const std::string &fname) : in_(fname) {}
    virtual ~FastaReader() {}

    bool IsValid() const { return in_.is_open(); }
    bool IsFileEnd() { return Tell() == -1; }

    virtual bool Next(Item &item);
    virtual bool Get(ItemId id, Item &item) { Seek(id); return Next(item);}

protected:
    bool GetHead(std::string &head, std::string &sub_head);
    bool GetSeq(std::string &seq);

    std::string NextNonEmptyLine();
    void ConsumeNextNonEmptyLine() { next_line_ = ""; pos_ = -1; }

    ItemId Tell() { return pos_ == -1 ? (std::streamoff)in_.tellg() : pos_; }
    void Seek(ItemId id) { next_line_.clear(); pos_ = -1; in_.clear(); in_.seekg(id, std::ios::beg); }

protected:
    std::ifstream in_;
    std::string next_line_ { "" };
    std::streamoff pos_ {-1} ;
};


#endif // FSA_FASTA_READER_HPP
