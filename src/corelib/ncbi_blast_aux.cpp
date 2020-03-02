#include "ncbi_blast_aux.hpp"

#include <cctype>

BEGIN_NCBI_SCOPE

const string kEmptyStr = "";
const size_t kDfltLineLength = 60;
const string kArgNumDescriptions("num_descriptions");
const size_t kDfltArgNumDescriptions = 500;
const string kArgNumAlignments("num_alignments");
const size_t kDfltArgNumAlignments = 250;

inline bool s_IsArgNameChar(char c)
{
    return isalnum(c)  ||  c == '_'  ||  c == '-';
}

bool VerifyArgumentName(const string& name, bool extended)
{
    if ( name.empty() )
        return true;

    string::const_iterator it = name.begin();
    if (extended  &&  *it == '#') {
        for (++it;  it != name.end();  ++it) {
            if ( !isdigit((unsigned char)(*it)) ) {
                return false;
            }
        }
    } else {
        if (name[0] == '-') {
            // Prohibit names like '-' or '--foo'.
            // The second char must be present and may not be '-'.
            if (name.length() == 1  ||  name[1] == '-') {
                return false;
            }
        }
        for ( ;  it != name.end();  ++it) {
            if ( !s_IsArgNameChar((unsigned char)(*it)) )
                return false;
        }
    }

    return true;
}

END_NCBI_SCOPE