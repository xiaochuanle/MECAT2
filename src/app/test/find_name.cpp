#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <string>

using namespace std;


extern "C"
int find_chr21_name(const char* name)
{
    const char* chr21 = "chromosome 21";
    string s(name);
    return s.find(chr21) != string::npos;
}