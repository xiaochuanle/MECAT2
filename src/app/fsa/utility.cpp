#include "utility.hpp"

#include <algorithm>
#include <string>
#include <vector>

std::vector<std::string> SplitStringBySpace(const std::string &str) {
    std::vector<std::string> substrs;

    auto IsSpace = [](char c) { return c == ' ' || c == '\t'; };
    auto IsNotSpace = [](char c) { return c != ' ' && c != '\t'; };

    auto begin = std::find_if(str.begin(), str.end(), IsNotSpace);

    while (begin != str.end()) {
        auto end = std::find_if(begin, str.end(), IsSpace);
        substrs.push_back(std::string(begin, end));
        begin = std::find_if(end, str.end(), IsNotSpace);
    }

    return substrs;
}
