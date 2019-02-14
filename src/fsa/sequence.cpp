#include "sequence.hpp"


std::string Seq::ReverseComplement(const std::string &seq) {
    std::string r(seq.rbegin(), seq.rend());

    for (size_t i = 0; i < r.length(); ++i) {
        switch (r[i]) {
        case 'A': r[i] = 'T'; break;
        case 'C': r[i] = 'G'; break;
        case 'G': r[i] = 'C'; break;
        case 'T': r[i] = 'A'; break;
        case 'a': r[i] = 't'; break;
        case 'c': r[i] = 'g'; break;
        case 'g': r[i] = 'c'; break;
        case 't': r[i] = 'a'; break;
        default: break;
        }
    }
    return r;
}

std::string SeqReader::NextNonEmptyLine(std::ifstream &in) {
    std::string line;
    bool r = (bool)std::getline(in, line);
    while (r) {
        std::string::size_type s = 0;
        while (s< line.size() && ::isspace(line[s])) s++;
        
        std::string::size_type e = line.size();
        while (e > s && ::isspace(line[e-1])) e--;
        
        if (e > s) return line.substr(s, e-s);
        
        r = (bool)std::getline(in, line);
    }
    return std::string();
}