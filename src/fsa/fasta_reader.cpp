#include "fasta_reader.hpp"

#include <cassert>
#include <algorithm>
#include <cassert>



bool FastaReader::Next(Item &item) {
    item.id = Tell();
    item.quality = "";
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq);
}

bool FastaReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = NextNonEmptyLine();
    ConsumeNextNonEmptyLine();
    if (!line.empty() && line[0] == '>') {
        std::string::size_type s = 1;
        while (s < line.size() && ::isspace(line[s])) s++;

        std::string::size_type e = std::min(s+1, line.size());
        while (e < line.size() && !::isspace(line[e])) e++;

        if (e > s) {
            head = line.substr(s, e-s);
            while (e < line.size() && ::isspace(line[e])) e++;
            sub_head = line.substr(e);
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool FastaReader::GetSeq(std::string &seq) {
    seq = "";
    std::string line = NextNonEmptyLine();

    while (!line.empty() && line[0] != '>') {
        seq += line;
        ConsumeNextNonEmptyLine();
        line = NextNonEmptyLine();
    }

    return true;// !seq.empty(); allow empty sequence

}

std::string FastaReader::NextNonEmptyLine() {
    if (next_line_.empty()) {
        pos_ = in_.tellg();
        next_line_ = SeqReader::NextNonEmptyLine(in_);
    }
    return next_line_;
}