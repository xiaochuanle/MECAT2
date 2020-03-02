#include "fastq_reader.hpp"

#include <algorithm>
#include <cassert>

FastqReader::FastqReader(const std::string &fname) : in_(fname) {
}

bool FastqReader::Next(Item &item) {
    assert(IsValid());

    item.id = Tell();
    return GetHead(item.head, item.sub_head) && GetSeq(item.seq) &&
           GetHead1() && GetQuality(item.quality);
}

bool FastqReader::GetHead(std::string &head, std::string &sub_head) {
    std::string line = NextNonEmptyLine();

    if (!line.empty() && line[0] == '@') {
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
