#include "str_util.hpp"

#include "../../corelib/hbn_aux.h"
#include "ncbistr_util.hpp"

#include <sstream>

BEGIN_HBNSTR_SCOPE

using namespace std;

template <class TStr>
TStr s_TruncateSpaces(const TStr& str, NStr::ETrunc where,
                      const TStr& empty_str)
{
    SIZE_TYPE length = str.length();
    if (length == 0) {
        return empty_str;
    }
    SIZE_TYPE beg = 0;
    if (where == eTrunc_Begin  ||  where == eTrunc_Both) {
        _ASSERT(beg < length);
        while ( isspace((unsigned char) str[beg]) ) {
            if (++beg == length) {
                return empty_str;
            }
        }
    }
    SIZE_TYPE end = length;
    if ( where == eTrunc_End  ||  where == eTrunc_Both ) {
        _ASSERT(beg < end);
        while (isspace((unsigned char) str[--end])) {
            if (beg == end) {
                return empty_str;
            }
        }
        _ASSERT(beg <= end  &&  !isspace((unsigned char) str[end]));
        ++end;
    }
    _ASSERT(beg < end  &&  end <= length);
    if ( beg | (end - length) ) { // if either beg != 0 or end != length
        return str.substr(beg, end - beg);
    }
    else {
        return str;
    }
}

string TruncateSpaces(const string& str, ETrunc where)
{
    return s_TruncateSpaces(str, where, kEmptyStr);
}

CTempString TruncateSpaces_Unsafe(const CTempString str, ETrunc where)
{
    return s_TruncateSpaces(str, where, CTempString());
}

void TruncateSpacesInPlace(CTempString& str, ETrunc where)
{
    str = s_TruncateSpaces(str, where, CTempString());
}

void TruncateSpacesInPlace(string& str, ETrunc where)
{
    SIZE_TYPE length = str.length();
    if (length == 0) {
        return;
    }
    SIZE_TYPE beg = 0;
    if ( where == ETrunc::eTrunc_Begin  ||  where == ETrunc::eTrunc_Both ) {
        // It's better to use str.data()[] to check string characters
        // to avoid implicit modification of the string by non-const operator[]
        _ASSERT(beg < length);
        while ( isspace((unsigned char) str.data()[beg]) ) {
            if (++beg == length) {
                str.erase();
                return;
            }
        }
    }

    SIZE_TYPE end = length;
    if ( where == ETrunc::eTrunc_End  ||  where == ETrunc::eTrunc_Both ) {
        // It's better to use str.data()[] to check string characters
        // to avoid implicit modification of the string by non-const operator[]
        _ASSERT(beg < end);
        while (isspace((unsigned char) str.data()[--end])) {
            if (beg == end) {
                str.erase();
                return;
            }
        }
        _ASSERT(beg <= end  &&  !isspace((unsigned char) str.data()[end]));
        ++end;
    }
    _ASSERT(beg < end  &&  end <= length);

#if defined(NCBI_COMPILER_GCC)  &&  (NCBI_COMPILER_VERSION == 304)
    // work around a library bug
    str.replace(end, length, kEmptyStr);
    str.replace(0, beg, kEmptyStr);
#else
    if ( beg | (end - length) ) { // if either beg != 0 or end != length
        str.replace(0, length, str, beg, end - beg);
    }
#endif
}

template<typename TString, typename TContainer>
TContainer& s_Split(const TString& str, const TString& delim,
                    TContainer& arr, NStr::TSplitFlags flags,
                    vector<SIZE_TYPE>* token_pos,
                    CTempString_Storage* storage = NULL)
{
    typedef CStrTokenPosAdapter<vector<SIZE_TYPE> >         TPosArray;
    typedef CStrDummyTargetReserve<TContainer, TPosArray>   TReserve;
    typedef CStrTokenize<TString, TContainer, TPosArray,
                         CStrDummyTokenCount, TReserve>     TSplitter;

    TPosArray token_pos_proxy(token_pos);
    TSplitter splitter(str, delim, flags, storage);
    splitter.Do(arr, token_pos_proxy, kEmptyStr);
    return arr;
}

#define CHECK_SPLIT_TEMPSTRING_FLAGS(where) \
    { \
        if ((flags & (NStr::fSplit_CanEscape | NStr::fSplit_CanQuote)) && !storage) { \
            string err = "NStr::" #where "(): the selected flags require non-NULL storage"; \
            HBN_ERR("%s", err.c_str()); \
    } \
}

list<string>& Split(const CTempString str, const CTempString delim,
                          list<string>& arr, TSplitFlags flags,
                          vector<SIZE_TYPE>* token_pos)
{
    return s_Split(str, delim, arr, flags, token_pos);
}

vector<string>& Split(const CTempString str, const CTempString delim,
                            vector<string>& arr, TSplitFlags flags,
                            vector<SIZE_TYPE>* token_pos)
{
    return s_Split(str, delim, arr, flags, token_pos);
}

list<CTempString>& Split(const CTempString str, const CTempString delim,
                               list<CTempString>& arr, TSplitFlags flags,
                               vector<SIZE_TYPE>* token_pos, CTempString_Storage* storage)
{
    CHECK_SPLIT_TEMPSTRING_FLAGS(Split);
    return s_Split(str, delim, arr, flags, token_pos, storage);
}

vector<CTempString>& Split(const CTempString str, const CTempString delim,
                                 vector<CTempString>& arr, TSplitFlags flags,
                                 vector<SIZE_TYPE>* token_pos, CTempString_Storage* storage)
{
    CHECK_SPLIT_TEMPSTRING_FLAGS(Split);
    return s_Split(str, delim, arr, flags, token_pos, storage);
}

list<CTempStringEx>& Split(const CTempString str, const CTempString delim,
                                 list<CTempStringEx>& arr, TSplitFlags flags,
                                 vector<SIZE_TYPE>* token_pos, CTempString_Storage* storage)
{
    CHECK_SPLIT_TEMPSTRING_FLAGS(Split);
    return s_Split(str, delim, arr, flags, token_pos, storage);
}

vector<CTempStringEx>& Split(const CTempString str, const CTempString delim,
                                   vector<CTempStringEx>& arr, TSplitFlags flags,
                                   vector<SIZE_TYPE>* token_pos, CTempString_Storage* storage)
{
    CHECK_SPLIT_TEMPSTRING_FLAGS(Split);
    return s_Split(str, delim, arr, flags, token_pos, storage);
}

bool SplitInTwo(const CTempString str, const CTempString delim,
                      string& str1, string& str2, TSplitFlags flags)
{
    CTempStringEx ts1, ts2;
    CTempString_Storage storage;
    bool result = SplitInTwo(str, delim, ts1, ts2, flags, &storage);
    str1 = ts1;
    str2 = ts2;
    return result;
}


bool SplitInTwo(const CTempString str, const CTempString delim,
                      CTempString& str1, CTempString& str2, TSplitFlags flags,
                      CTempString_Storage* storage)
{
    CTempStringEx ts1, ts2;
    bool result = SplitInTwo(str, delim, ts1, ts2, flags, storage);
    str1 = ts1;
    str2 = ts2;
    return result;
}


bool SplitInTwo(const CTempString str, const CTempString delim,
                      CTempStringEx& str1, CTempStringEx& str2,
                      TSplitFlags flags, CTempString_Storage* storage)
{
    CHECK_SPLIT_TEMPSTRING_FLAGS(SplitInTwo);
    typedef CStrTokenize<CTempString, int, CStrDummyTokenPos,
                         CStrDummyTokenCount,
                         CStrDummyTargetReserve<int, int> > TSplitter;

    CTempStringList part_collector(storage);
    TSplitter       splitter(str, delim, flags, storage);
    SIZE_TYPE       delim_pos = NPOS;

    // get first part
    splitter.Advance(&part_collector, NULL, &delim_pos);
    part_collector.Join(&str1);
    part_collector.Clear();

    // don't need further splitting, just quote and escape parsing
    splitter.SetDelim(kEmptyStr);
    splitter.Advance(&part_collector);
    part_collector.Join(&str2);

    return delim_pos != NPOS;
}

// Determines the end of an HTML <...> tag, accounting for attributes
// and comments (the latter allowed only within <!...>).
static SIZE_TYPE s_EndOfTag(const string& str, SIZE_TYPE start)
{
    _ASSERT(start < str.size()  &&  str[start] == '<');
    bool comments_ok = (start + 1 < str.size()  &&  str[start + 1] == '!');
    for (SIZE_TYPE pos = start + 1;  pos < str.size();  ++pos) {
        switch (str[pos]) {
        case '>': // found the end
            return pos;

        case '\"': // start of "string"; advance to end
            pos = str.find('\"', pos + 1);
            if (pos == NPOS) {
                HBN_ERR("Unclosed string in HTML tag at position %zu", start);
                // return pos;
            }
            break;

        case '-': // possible start of -- comment --; advance to end
            if (comments_ok  &&  pos + 1 < str.size()
                &&  str[pos + 1] == '-') {
                pos = str.find("--", pos + 2);
                if (pos == NPOS) {
                    HBN_ERR("Unclosed comment in HTML tag at position %zu", start);
                    // return pos;
                } else {
                    ++pos;
                }
            }
        }
    }
    HBN_ERR("Unclosed HTML tag at position %zu", start);
    // return NPOS;
}


// Determines the end of an HTML &foo; character/entity reference
// (which might not actually end with a semicolon  :-/ , but we ignore that case)
static SIZE_TYPE s_EndOfReference(const string& str, SIZE_TYPE start)
{
    _ASSERT(start < str.size()  &&  str[start] == '&');

    SIZE_TYPE pos = str.find_first_not_of
        ("#0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz",
         start + 1);
    if (pos != NPOS  &&  str[pos] == ';') {
        // found terminating semicolon, so it's valid, and we return that
        return pos;
    } else {
        // We consider it just a '&' by itself since it's invalid
        return start;
    }
}


static SIZE_TYPE s_VisibleHtmlWidth(const string& str)
{
    SIZE_TYPE width = 0, pos = 0;
    for (;;) {
        SIZE_TYPE pos2 = str.find_first_of("<&", pos);
        if (pos2 == NPOS) {
            width += str.size() - pos;
            break;
        } else {
            width += pos2 - pos;
            if (str[pos2] == '&') {
                ++width;
                pos = s_EndOfReference(str, pos);
            } else {
                pos = s_EndOfTag(str, pos);
            }
            if (pos == NPOS) {
                break;
            } else {
                ++pos;
            }
        }
    }
    return width;
}

static
inline bool _isspace(unsigned char c)
{
    return ((c>=0x09 && c<=0x0D) || (c==0x20));
}

template<typename _D>
void WrapIt(const string& str, SIZE_TYPE width,
    _D& dest, TWrapFlags flags,
    const string* prefix,
    const string* prefix1)
{
    if (prefix == 0) {
        prefix = &kEmptyStr;
    }

    if (prefix1 == 0)
        prefix1 = prefix;

    SIZE_TYPE     pos = 0, len = str.size(), nl_pos = 0;

    const bool          is_html = flags & fWrap_HTMLPre ? true : false;
    const bool          do_flat = (flags & fWrap_FlatFile) != 0;
    string temp_back; temp_back.reserve(width);

    enum EScore { // worst to best
        eForced,
        ePunct,
        eComma,
        eSpace,
        eNewline
    };

    // To avoid copying parts of str when we need to store a 
    // substr of str, we store the substr as a pair
    // representing start (inclusive) and end (exclusive).
    typedef pair<SIZE_TYPE, SIZE_TYPE> TWrapSubstr;

    // This variable is used for HTML links that cross line boundaries.
    // Since it's aesthetically displeasing for a link to cross a boundary, we 
    // close it at the end of each line and re-open it after the next line's 
    // prefix
    // (This is needed in, e.g. AE017351)
    TWrapSubstr best_link(0, 0); // last link found before current best_pos
    TWrapSubstr latest_link(0, 0); // last link found at all

    while (pos < len) {
        bool      hyphen = false; // "-" or empty
        SIZE_TYPE column = is_html ? s_VisibleHtmlWidth(*prefix1) : prefix1->size();
        SIZE_TYPE column0 = column;
        // the next line will start at best_pos
        SIZE_TYPE best_pos = NPOS;
        EScore    best_score = eForced;

        // certain logic can be skipped if this part has no backspace,
        // which is, by far, the most common case
        bool thisPartHasBackspace = false;

        temp_back = *prefix1;

        // append any still-open links from previous lines
        if (is_html && best_link.second != 0) {
            temp_back.append(
                str.begin() + best_link.first,
                str.begin() + best_link.second);
        }

        SIZE_TYPE pos0 = pos;

        // we can't do this in HTML mode because we might have to deal with
        // link tags that go across lines.
        if (!is_html) {
            if (nl_pos <= pos) {
                nl_pos = str.find('\n', pos);
                if (nl_pos == NPOS) {
                    nl_pos = len;
                }
            }
            if (column + (nl_pos - pos) <= width) {
                pos0 = nl_pos;
            }
        }

        for (SIZE_TYPE pos2 = pos0;  pos2 < len  &&  column <= width;
            ++pos2, ++column) {
            EScore    score = eForced;
            SIZE_TYPE score_pos = pos2;
            const char      c = str[pos2];

            if (c == '\n') {
                best_pos = pos2;
                best_score = eNewline;
                best_link = latest_link;
                break;
            }
            else if (_isspace((unsigned char)c)) {
                if (!do_flat  &&  pos2 > 0 &&
                    _isspace((unsigned char)str[pos2 - 1])) {
                    if (pos2 < len - 1 && str[pos2 + 1] == '\b') {
                        thisPartHasBackspace = true;
                    }
                    continue; // take the first space of a group
                }
                score = eSpace;
            }
            else if (is_html && c == '<') {
                // treat tags as zero-width...
                SIZE_TYPE start_of_tag = pos2;
                pos2 = s_EndOfTag(str, pos2);
                --column;
                if (pos2 == NPOS) {
                    break;
                }

                if ((pos2 - start_of_tag) >= 6 &&
                    str[start_of_tag + 1] == 'a' &&
                    str[start_of_tag + 2] == ' ' &&
                    str[start_of_tag + 3] == 'h' &&
                    str[start_of_tag + 4] == 'r' &&
                    str[start_of_tag + 5] == 'e' &&
                    str[start_of_tag + 6] == 'f')
                {
                    // remember current link in case of line wrap
                    latest_link.first = start_of_tag;
                    latest_link.second = pos2 + 1;
                }
                if ((pos2 - start_of_tag) >= 3 &&
                    str[start_of_tag + 1] == '/' &&
                    str[start_of_tag + 2] == 'a' &&
                    str[start_of_tag + 3] == '>')
                {
                    // link is closed
                    latest_link.first = 0;
                    latest_link.second = 0;
                }
            }
            else if (is_html && c == '&') {
                // ...and references as single characters
                pos2 = s_EndOfReference(str, pos2);
                if (pos2 == NPOS) {
                    break;
                }
            }
            else if (c == ','  && column < width && score_pos < len - 1) {
                score = eComma;
                ++score_pos;
            }
            else if (do_flat ? c == '-' : ispunct((unsigned char)c)) {
                // For flat files, only whitespace, hyphens and commas
                // are special.
                switch (c) {
                case '(': case '[': case '{': case '<': case '`':
                    score = ePunct;
                    break;
                default:
                    if (score_pos < len - 1 && column < width) {
                        score = ePunct;
                        ++score_pos;
                    }
                    break;
                }
            }

            if (score >= best_score  &&  score_pos > pos0) {
                best_pos = score_pos;
                best_score = score;
                best_link = latest_link;
            }

            while (pos2 < len - 1 && str[pos2 + 1] == '\b') {
                // Account for backspaces
                ++pos2;
                if (column > column0) {
                    --column;
                }
                thisPartHasBackspace = true;
            }
        }

        if (best_score != eNewline  &&  column <= width) {
            if (best_pos != len) {
                // If the whole remaining text can fit, don't split it...
                best_pos = len;
                best_link = latest_link;
                // Force backspace checking, to play it safe
                thisPartHasBackspace = true;
            }
        }
        else if (best_score == eForced && (flags & fWrap_Hyphenate)) {
            hyphen = true;
            --best_pos;
        }

        {{
            string::const_iterator begin = str.begin() + pos;
            string::const_iterator end = str.begin() + best_pos;
            if (thisPartHasBackspace) {
                // eat backspaces and the characters (if any) that precede them

                string::const_iterator bs; // position of next backspace
                while ((bs = find(begin, end, '\b')) != end) {
                    if (bs != begin) {
                        // add all except the last one
                        temp_back.append(begin, bs - 1);
                    }
                    else {
                        // The backspace is at the beginning of next substring,
                        // so we should remove previously added symbol if any.
                        SIZE_TYPE size = temp_back.size();
                        if (size > prefix1->size()) { // current size > prefix size
                            temp_back.resize(size - 1);
                        }
                    }
                    // skip over backspace
                    begin = bs + 1;
                }
            }
            if (begin != end) {
                // add remaining characters
                temp_back.append(begin, end);
            }
        }}

        // if we didn't close the link on this line, we 
        // close it here
        if (is_html && best_link.second != 0) {
            temp_back += "</a>";
        }

        if (hyphen) {
            temp_back += '-';
        }
        pos = best_pos;
        prefix1 = prefix;

        if (do_flat) {
            if (best_score == eSpace) {
                while (str[pos] == ' ') {
                    ++pos;
                }
                if (str[pos] == '\n') {
                    ++pos;
                }
            }
            if (best_score == eNewline) {
                ++pos;
            }
        }
        else {
            if (best_score == eSpace || best_score == eNewline) {
                ++pos;
            }
        }
        while (pos < len  &&  str[pos] == '\b') {
            ++pos;
        }

        dest.Append(temp_back);
    }
}

void Wrap(const string& str, SIZE_TYPE width,
                IWrapDest& dest, TWrapFlags flags,
                const string* prefix,
                const string* prefix1)
{
    WrapIt(str, width, dest, flags, prefix, prefix1);
}


list<string>& Wrap(const string& str, SIZE_TYPE width,
                         list<string>& arr2, NStr::TWrapFlags flags,
                         const string* prefix, const string* prefix1)
{
    CWrapDestStringList d(arr2);
    WrapIt(str, width, d, flags, prefix, prefix1);
    return arr2;
}

list<string>& WrapList(const list<string>& l, SIZE_TYPE width,
                             const string& delim, list<string>& arr,
                             NStr::TWrapFlags flags,
                             const string* prefix,
                             const string* prefix1)
{
    if (l.empty()) {
        return arr;
    }

    const string* pfx      = prefix1 ? prefix1 : prefix;
    string        s        = *pfx;
    bool          is_html  = flags & fWrap_HTMLPre ? true : false;
    SIZE_TYPE     column   = is_html? s_VisibleHtmlWidth(s)     : s.size();
    SIZE_TYPE     delwidth = is_html? s_VisibleHtmlWidth(delim) : delim.size();
    bool          at_start = true;

    ITERATE (list<string>, it, l) {
        SIZE_TYPE term_width = is_html ? s_VisibleHtmlWidth(*it) : it->size();
        if ( at_start ) {
            if (column + term_width <= width) {
                s += *it;
                column += term_width;
                at_start = false;
            } else {
                // Can't fit, even on its own line; break separately.
                Wrap(*it, width, arr, flags, prefix, pfx);
                pfx      = prefix;
                s        = *prefix;
                column   = is_html ? s_VisibleHtmlWidth(s) : s.size();
                at_start = true;
            }
        } else if (column + delwidth + term_width <= width) {
            s += delim;
            s += *it;
            column += delwidth + term_width;
            at_start = false;
        } else {
            // Can't fit on this line; break here and try again.
            arr.push_back(s);
            pfx      = prefix;
            s        = *prefix;
            column   = is_html ? s_VisibleHtmlWidth(s) : s.size();
            at_start = true;
            --it;
        }
    }
    arr.push_back(s);
    return arr;
}

enum ELanguage {
    eLanguage_C,
    eLanguage_Javascript
};


static string s_PrintableString(const CTempString    str,
                                NStr::TPrintableMode mode,
                                ELanguage            lang)
{
    //unique_ptr<CNcbiOstrstream> out;
    ostringstream out;
    SIZE_TYPE i, j = 0;

    for (i = 0;  i < str.size();  i++) {
        bool octal = false;
        char c = str[i];
        switch (c) {
        case '\t':
            c = 't';
            break;
        case '\v':
            c = 'v';
            break;
        case '\b':
            c = 'b';
            break;
        case '\r':
            c = 'r';
            break;
        case '\f':
            c = 'f';
            break;
        case '\a':
            c = 'a';
            break;
        case '\n':
            if (!(mode & NStr::fNewLine_Passthru))
                c = 'n';
            /*FALLTHRU*/
        case '\\':
        case '\'':
        case '"':
            break;
        case '&':
            if (lang != eLanguage_Javascript)
                continue;
            break;
        default:
            if (!isascii((unsigned char) c)) {
                if (mode & NStr::fNonAscii_Quote) {
                    octal = true;
                    break;
                }
            }
            if (!isprint((unsigned char) c)) {
                octal = true;
                break;
            }
            continue;
        }
        //if (!out.get()) {
        //    out.reset(new CNcbiOstrstream);
        //}
        if (i > j) {
            //out->write(str.data() + j, i - j);
            out.write(str.data() + j, i - j);
        }
        //out->put('\\');
        out.put('\\');
        if (c == '\n') {
            //out->write("n\\\n", 3);
            out.write("n\\\n", 3);
        } else if (octal) {
            bool reduce;
            if (!(mode & NStr::fPrintable_Full)) {
                reduce = (i == str.size() - 1  ||
                          str[i + 1] < '0' || str[i + 1] > '7' ? true : false);
            } else {
                reduce = false;
            }
            unsigned char v;
            char val[3];
            int k = 0;
            v = (unsigned char)((unsigned char)c >> 6);
            if (v  ||  !reduce) {
                val[k++] = char('0' + v);
                reduce = false;
            }
            v = ((unsigned char)c >> 3) & 7;
            if (v  ||  !reduce) {
                val[k++] = char('0' + v);
            }
            v = (unsigned char)c & 7;
            val[k++] = char('0' + v);
            //out->write(val, k);
            out.write(val, k);
        } else {
            //out->put(c);
            out.put(c);
        }
        j = i + 1;
    }
    if (j  &&  i > j) {
        //_ASSERT(out.get());
        //out->write(str.data() + j, i - j);
        _ASSERT(!out.str().empty());
        out.write(str.data() + j, i - j);
    }
    //if (out.get()) {
    //    // Return encoded string
    //    return CNcbiOstrstreamToString(*out);
    //}
    if (!out.str().empty()) return out.str();

    // All characters are good - return (a copy of) the original string
    return str;
}

string PrintableString(const CTempString str, NStr::TPrintableMode mode)
{
    return s_PrintableString(str, mode, eLanguage_C);
}

string& Replace(const string& src,
                      const string& search, const string& replace,
                      string& dst, SIZE_TYPE start_pos, SIZE_TYPE max_replace,
                      SIZE_TYPE* num_replace)
{
    // source and destination should not be the same
    if (&src == &dst) {
        //NCBI_THROW2(CStringException, eBadArgs,
        HBN_ERR("NStr::Replace():  source and destination are the same");
    }
    if (num_replace)
        *num_replace = 0;
    if (start_pos + search.size() > src.size() || search == replace) {
        dst = src;
        return dst;
    }

    // Use different algorithms depending on size or 'search' and 'replace'
    // for better performance (and for big strings only! > 16KB).

    if (replace.size() > search.size()  &&  src.size() > 16*1024) {
        // Replacing string is longer -- worst case.
        // Try to avoid memory reallocations inside std::string.
        // Count replacing strings first
        SIZE_TYPE n = 0;
        SIZE_TYPE start_orig = start_pos;
        for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
            start_pos = src.find(search, start_pos);
            if (start_pos == NPOS)
                break;
            n++;
            start_pos += search.size();
        }
        // Reallocate memory for destination string
        dst.resize(src.size() - n*search.size() + n*replace.size());

        // Use copy() to create destination string
        start_pos = start_orig;
        string::const_iterator src_start = src.begin();
        string::const_iterator src_end   = src.begin();
        string::iterator       dst_pos   = dst.begin();

        for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
            start_pos = src.find(search, start_pos);
            if (start_pos == NPOS)
                break;
            // Copy from source string up to 'search'
            src_end = src.begin() + start_pos;
            copy(src_start, src_end, dst_pos); 
            dst_pos += (src_end - src_start);
            // Append 'replace'
            copy(replace.begin(), replace.end(), dst_pos); 
            dst_pos   += replace.size();
            start_pos += search.size();
            src_start = src.begin() + start_pos;
        }
        // Copy source's string tail to the place
        copy(src_start, src.end(), dst_pos); 
        if (num_replace)
            *num_replace = n;

    } else {
        // Replacing string is shorter or have the same length.
        // ReplaceInPlace() can be faster on some platform, but not much,
        // so we use regular algorithm even for equal lengths here.
        dst = src;
        for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
            start_pos = dst.find(search, start_pos);
            if (start_pos == NPOS)
                break;
            dst.replace(start_pos, search.size(), replace);
            start_pos += replace.size();
            if (num_replace)
                (*num_replace)++;
        }
    }
    return dst;
}


string Replace(const string& src,
                     const string& search, const string& replace,
                     SIZE_TYPE start_pos, SIZE_TYPE max_replace,
                     SIZE_TYPE* num_replace)
{
    string dst;
    Replace(src, search, replace, dst, start_pos, max_replace, num_replace);
    return dst;
}

END_HBNSTR_SCOPE