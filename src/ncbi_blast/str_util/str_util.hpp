#ifndef __STR_UTIL_HPP
#define __STR_UTIL_HPP

#include "../ncbi_blast_aux.hpp"
#include "tempstr.hpp"
//#include "ncbistr_util.hpp"

#include <vector>
#include <initializer_list>

BEGIN_HBNSTR_SCOPE

    using namespace ncbi;
    using std::vector;

    /// Which end to truncate a string.
    enum ETrunc {
        eTrunc_Begin,  ///< Truncate leading spaces only
        eTrunc_End,    ///< Truncate trailing spaces only
        eTrunc_Both    ///< Truncate spaces at both begin and end of string
    };

    /// Truncate spaces in a string.
    ///
    /// @param str
    ///   String to truncate spaces from.
    /// @param where
    ///   Which end of the string to truncate space from. Default is to
    ///   truncate space from both ends (eTrunc_Both).
    /// @sa
    ///   TruncateSpaces_Unsafe
    string TruncateSpaces(const string& str,
                                 ETrunc        where = ETrunc::eTrunc_Both);

    /// Truncate spaces in a string.
    /// It can be faster but it is also more dangerous than TruncateSpaces()
    ///
    /// @param str
    ///   String to truncate spaces from.
    /// @param where
    ///   Which end of the string to truncate space from. Default is to
    ///   truncate space from both ends (eTrunc_Both).
    /// @attention
    ///   The lifespan of the result string is the same as one of the source.
    ///   So, for example, if the source is temporary string, or it changes somehow,
    ///   then the result will be invalid right away (will point to already released
    ///   or wrong range in the memory).
    /// @sa
    ///   TruncateSpaces
    CTempString TruncateSpaces_Unsafe(const CTempString str,
                                             ETrunc where = ETrunc::eTrunc_Both);

    /// @deprecated  Use TruncateSpaces_Unsafe() instead -- AND, do make sure
    ///              that you indeed use that in a safe manner!
    inline
    NCBI_DEPRECATED
    CTempString TruncateSpaces(const CTempString str,
                                      ETrunc where = ETrunc::eTrunc_Both) {
        return TruncateSpaces_Unsafe(str, where);
    }      

    /// @deprecated  Use TruncateSpaces_Unsafe() instead -- AND, do make sure
    ///              that you indeed use that in a safe manner!
    inline
    NCBI_DEPRECATED
    CTempString TruncateSpaces(const char* str,
                                      ETrunc where = ETrunc::eTrunc_Both) {
        return TruncateSpaces_Unsafe(str, where);
    }      

    /// Truncate spaces in a string (in-place)
    ///
    /// @param str
    ///   String to truncate spaces from.
    /// @param where
    ///   Which end of the string to truncate space from. Default is to
    ///   truncate space from both ends (eTrunc_Both).
    void TruncateSpacesInPlace(string& str,  ETrunc where = ETrunc::eTrunc_Both);
    void TruncateSpacesInPlace(CTempString&, ETrunc where = ETrunc::eTrunc_Both);     

    /// Whether to merge adjacent delimiters.
    /// Used by some methods that don't need full functionality of ESplitFlags.
    enum EMergeDelims {
        eMergeDelims   = fSplit_MergeDelimiters | fSplit_Truncate,
        eNoMergeDelims = 0
    };

    /// Split a string using specified delimiters.
    ///
    /// @param str
    ///   String to be split.
    /// @param delim
    ///   Delimiter(s) used to split string "str". The interpretation of
    ///   multi-character values depends on flags: by default, any of those
    ///   characters marks a split point (when unquoted), but with
    ///   fSplit_ByPattern, the entire string must occur.  (Meanwhile,
    ///   an empty value disables splitting.)
    /// @param arr
    ///   The split tokens are added to the list "arr" and also returned
    ///   by the function.
    /// @param flags
    ///   Flags directing splitting, characterized under ESplitFlags.
    /// @param token_pos
    ///   Optional array for the tokens' positions in "str".
    /// @attention
    ///   Modifying source CTempString object or destroying it,
    ///   will invalidate results.
    /// @return 
    ///   The list "arr" is also returned.
    /// @sa
    ///   ESplitFlags, SplitInTwo, SplitByPattern
    list<string>& Split( const CTempString    str,
                                const CTempString    delim,
                                list<string>&        arr,
                                TSplitFlags          flags = 0,
                                vector<SIZE_TYPE>*   token_pos = NULL);  

    vector<string>& Split(
                                const CTempString    str,
                                const CTempString    delim,
                                vector<string>&      arr,
                                TSplitFlags          flags = 0,
                                vector<SIZE_TYPE>*   token_pos = NULL); 

    list<CTempString>& Split(
                                const CTempString    str,
                                const CTempString    delim,
                                list<CTempString>&   arr,
                                TSplitFlags          flags = 0,
                                vector<SIZE_TYPE>*   token_pos = NULL,
                                CTempString_Storage* storage = NULL);  

    vector<CTempString>& Split(
                                const CTempString    str,
                                const CTempString    delim,
                                vector<CTempString>& arr,
                                TSplitFlags          flags = 0,
                                vector<SIZE_TYPE>*   token_pos = NULL,
                                CTempString_Storage* storage = NULL);

    list<CTempStringEx>& Split(
                                const CTempString    str,
                                const CTempString    delim,
                                list<CTempStringEx>& arr,
                                TSplitFlags          flags = 0,
                                vector<SIZE_TYPE>*   token_pos = NULL,
                                CTempString_Storage* storage = NULL);

    vector<CTempStringEx>& Split(
                                const CTempString      str,
                                const CTempString      delim,
                                vector<CTempStringEx>& arr,
                                TSplitFlags            flags = 0,
                                vector<SIZE_TYPE>*     token_pos = NULL,
                                CTempString_Storage*   storage = NULL);  

    /// Split a string into two pieces using the specified delimiters
    ///
    /// @param str 
    ///   String to be split.
    /// @param delim
    ///   Delimiters used to split string "str".
    /// @param str1
    ///   The sub-string of "str" before the first character of "delim".
    ///   It will not contain any characters in "delim".
    ///   Will be empty if "str" begin with a delimiter.
    /// @param str2
    ///   The sub-string of "str" after the first character of "delim" found.
    ///   May contain "delim" characters.
    ///   Will be empty if "str" had no "delim" characters or ended
    ///   with the "delim" character.
    /// @param flags
    ///   Flags directing splitting, characterized under ESplitFlags.
    ///   Note, that fSplit_Truncate_End don't have any effect due nature
    ///   of this method.
    /// @attention
    ///   Modifying source CTempString object or destroying it,
    ///   will invalidate results.
    /// @return
    ///   true if a symbol from "delim" was found in "str", false if not.
    ///   This lets you distinguish when there were no delimiters and when
    ///   the very last character was the first delimiter.
    /// @sa
    ///   ESplitFlags, Split
    bool SplitInTwo(const CTempString  str, 
                           const CTempString  delim,
                           string&            str1,
                           string&            str2,
                           TSplitFlags        flags = 0);

    bool SplitInTwo(const CTempString  str, 
                           const CTempString  delim,
                           CTempString&       str1,
                           CTempString&       str2,
                           TSplitFlags        flags = 0,
                           CTempString_Storage* storage = NULL);

    bool SplitInTwo(const CTempString  str, 
                           const CTempString  delim,
                           CTempStringEx&     str1,
                           CTempStringEx&     str2,
                           TSplitFlags        flags = 0,
                           CTempString_Storage* storage = NULL);

    /// Join strings using the specified delimiter.
    ///
    /// @param arr
    ///   Array of strings to be joined.
    /// @param delim
    ///   Delimiter used to join the string.
    /// @return 
    ///   The strings in "arr" are joined into a single string, separated
    ///   with "delim".
    /// @sa Split

    template<typename TIterator>
    static string xx_Join( TIterator from, TIterator to, const CTempString& delim);

    using std::enable_if;
    using std::is_same;
    using std::input_iterator_tag;
    using std::is_convertible;
    using std::forward_iterator_tag;
    using std::is_arithmetic;

    template<typename TIterator>
    static typename enable_if<is_same<typename TIterator::iterator_category, input_iterator_tag>::value &&
                              is_convertible<typename TIterator::value_type, string>::value, string>::type
    x_Join( TIterator from, TIterator to, const CTempString& delim)
    {
        return TransformJoin(from, to, delim, [](const typename TIterator::value_type& i){ return i;});
    }

    template<typename TIterator>
    static typename enable_if<is_convertible<typename TIterator::iterator_category, forward_iterator_tag>::value &&
                              is_convertible<typename TIterator::value_type, string>::value, string>::type
    x_Join( TIterator from, TIterator to, const CTempString& delim)
    {
        return xx_Join(from, to, delim);
    }

    template<typename TValue>
    static typename enable_if<is_convertible<TValue, string>::value, string>::type
    x_Join( TValue* from, TValue* to, const CTempString& delim)
    {
        return xx_Join(from, to, delim);
    }

    template<typename TIterator>
    static typename enable_if<is_convertible<typename TIterator::iterator_category, std::input_iterator_tag>::value &&
                              is_arithmetic< typename TIterator::value_type>::value, string>::type
    x_Join( TIterator from, TIterator to, const CTempString& delim)
    {
        return TransformJoin( from, to, delim, [](const typename TIterator::value_type& i){ return NumericToString(i);});
    }

    template<typename TValue>
    static typename enable_if<is_arithmetic<TValue>::value, string>::type
    x_Join( TValue* from, TValue* to, const CTempString& delim)
    {
        return TransformJoin( from, to, delim, [](const TValue& i){ return NumericToString(i);});
    }

    template<typename TContainer>
    string
    Join(const TContainer& arr, const CTempString& delim)
    {
        return x_Join(begin(arr), end(arr), delim);
    }
    template<typename TValue>
    string
    Join(const std::initializer_list<TValue>& arr, const CTempString& delim)
    {
        return x_Join(begin(arr), end(arr), delim);
    }
    template<typename TInputIterator>
    string
    Join( TInputIterator from, TInputIterator to, const CTempString& delim)
    {
        return x_Join(from, to, delim);
    }
    template<typename TInputIterator>
    string
    JoinNumeric( TInputIterator from, TInputIterator to, const CTempString& delim)
    {
        return x_Join( from, to, delim);
    }
    template<typename TIterator, typename FTransform>
    string
    TransformJoin( TIterator from, TIterator to, const CTempString& delim, FTransform fnTransform);



    /// How to wrap the words in a string to a new line.
    enum EWrapFlags {
        fWrap_Hyphenate  = 0x1, ///< Add a hyphen when breaking words?
        fWrap_HTMLPre    = 0x2, ///< Wrap as pre-formatted HTML?
        fWrap_FlatFile   = 0x4  ///< Wrap for flat file use.
    };
    typedef int TWrapFlags;     ///< Bitwise OR of "EWrapFlags"

    /// Wrap the specified string into lines of a specified width.
    ///
    /// Split string "str" into lines of width "width" and add the
    /// resulting lines to the list "arr".  Normally, all
    /// lines will begin with "prefix" (counted against "width"),
    /// but the first line will instead begin with "prefix1" if
    /// you supply it.
    ///
    /// @param str
    ///   String to be split into wrapped lines.
    /// @param width
    ///   Width of each wrapped line.
    /// @param arr
    ///   List of strings containing wrapped lines.
    /// @param flags
    ///   How to wrap the words to a new line. See EWrapFlags documentation.
    /// @param prefix
    ///   The prefix string added to each wrapped line, except the first line,
    ///   unless "prefix1" is set.
    ///   If "prefix" is set to 0(default), do not add a prefix string to the
    ///   wrapped lines.
    /// @param prefix1
    ///   The prefix string for the first line. Use this for the first line
    ///   instead of "prefix".
    ///   If "prefix1" is set to 0(default), do not add a prefix string to the
    ///   first line.
    /// @return
    ///   Return "arr", the list of wrapped lines.
    template<typename _D>
    static void WrapIt(const string& str, SIZE_TYPE width,
        _D& dest, TWrapFlags flags = 0,
        const string* prefix = 0,
        const string* prefix1 = 0);

    class IWrapDest
    {
    public:
        virtual ~IWrapDest() {}
        virtual void Append(const string& s) = 0;
        virtual void Append(const CTempString& s) = 0;
    };

    class CWrapDestStringList : public IWrapDest
    {
    protected:
        list<string>& m_list;
    public:
        CWrapDestStringList(list<string>& l) : m_list(l) {};
        virtual void Append(const string& s)
        {
            m_list.push_back(s);
        }
        virtual void Append(const CTempString& s)
        {
            m_list.push_back(kEmptyStr);
            m_list.back().assign(s.data(), s.length());
        }
    };

    void Wrap(const string& str, SIZE_TYPE width,
                              IWrapDest& dest, TWrapFlags flags,
                              const string* prefix,
                              const string* prefix1);

    list<string>& Wrap(const string& str, SIZE_TYPE width,
                              list<string>& arr, TWrapFlags flags = 0,
                              const string* prefix = 0,
                              const string* prefix1 = 0);

    list<string>& Wrap(const string& str, SIZE_TYPE width,
                              list<string>& arr, TWrapFlags flags,
                              const string& prefix,
                              const string* prefix1 = 0);

    list<string>& Wrap(const string& str, SIZE_TYPE width,
                              list<string>& arr, TWrapFlags flags,
                              const string& prefix,
                              const string& prefix1);

    /// Wrap the list using the specified criteria.
    ///
    /// WrapList() is similar to Wrap(), but tries to avoid splitting any
    /// elements of the list to be wrapped.  Also, the "delim" only applies
    /// between elements on the same line; if you want everything to end with
    /// commas or such, you should add them first.
    ///
    /// @param l
    ///   The list to be wrapped.
    /// @param width
    ///   Width of each wrapped line.
    /// @param delim
    ///   Delimiters used to split elements on the same line.
    /// @param arr
    ///   List containing the wrapped list result.
    /// @param flags
    ///   How to wrap the words to a new line. See EWrapFlags documentation.
    /// @param prefix
    ///   The prefix string added to each wrapped line, except the first line,
    ///   unless "prefix1" is set.
    ///   If "prefix" is set to 0(default), do not add a prefix string to the
    ///   wrapped lines.
    /// @param prefix1
    ///   The prefix string for the first line. Use this for the first line
    ///   instead of "prefix".
    ///   If "prefix1" is set to 0(default), do not add a prefix string to the
    ///   first line.
    /// @return
    ///   Return "arr", the wrapped list.
    list<string>& WrapList(const list<string>& l, SIZE_TYPE width,
                                  const string& delim, list<string>& arr,
                                  TWrapFlags    flags = 0,
                                  const string* prefix = 0,
                                  const string* prefix1 = 0);

    list<string>& WrapList(const list<string>& l, SIZE_TYPE width,
                                  const string& delim, list<string>& arr,
                                  TWrapFlags    flags,
                                  const string& prefix,
                                  const string* prefix1 = 0);
        
    list<string>& WrapList(const list<string>& l, SIZE_TYPE width,
                                  const string& delim, list<string>& arr,
                                  TWrapFlags    flags,
                                  const string& prefix,
                                  const string& prefix1);      

    /// How to display printable strings.
    ///
    /// Assists in making a printable version of "str".
    enum EPrintableMode {
        fNewLine_Quote     = 0,  ///< Display "\n" instead of actual linebreak
        eNewLine_Quote     = fNewLine_Quote,
        fNewLine_Passthru  = 1,  ///< Break the line at every "\n" occurrence
        eNewLine_Passthru  = fNewLine_Passthru,
        fNonAscii_Passthru = 0,  ///< Allow non-ASCII but printable characters
        fNonAscii_Quote    = 2,  ///< Octal for all non-ASCII characters
        fPrintable_Full    = 64  ///< Show all octal digits at all times
    };
    typedef int TPrintableMode;  ///< Bitwise OR of EPrintableMode flags

    /// Get a printable version of the specified string. 
    ///
    /// All non-printable characters will be represented as "\r", "\n", "\v",
    /// "\t", "\"", "\\\\", etc, or "\\ooo" where 'ooo' is an octal code of the
    /// character.  The resultant string is a well-formed C string literal,
    /// which, without alterations, can be compiled by a C/C++ compiler.
    /// In many instances, octal representations of non-printable characters
    /// can be reduced to take less than all 3 digits, if there is no
    /// ambiguity in the interpretation.  fPrintable_Full cancels the
    /// reduction, and forces to produce full 3-digit octal codes throughout.
    ///
    /// @param str
    ///   The string whose printable version is wanted.
    /// @param mode
    ///   How to display the string.  The default setting of fNewLine_Quote
    ///   displays the new lines as "\n", and uses the octal code reduction.
    ///   When set to fNewLine_Passthru, line breaks are actually
    ///   produced on output but preceded with trailing backslashes.
    /// @return
    ///   Return a printable version of "str".
    /// @sa
    ///   ParseEscapes, Escape, CEncode, CParse, Sanitize
    string PrintableString(const CTempString str,
                                  TPrintableMode    mode = fNewLine_Quote | fNonAscii_Passthru);  

    /// Replace occurrences of a substring within a string.
    ///
    /// @param src
    ///   Source string from which specified substring occurrences are
    ///   replaced.
    /// @param search
    ///   Substring value in "src" that is replaced.
    /// @param replace
    ///   Replace "search" substring with this value.
    /// @param dst
    ///   Result of replacing the "search" string with "replace" in "src".
    ///   This value is also returned by the function.
    /// @param start_pos
    ///   Position to start search from.
    /// @param max_replace
    ///   Replace no more than "max_replace" occurrences of substring "search"
    ///   If "max_replace" is zero(default), then replace all occurrences with
    ///   "replace".
    /// @param num_replace
    ///   Optional pointer to a value which receives number of replacements occurred.
    /// @return
    ///   Result of replacing the "search" string with "replace" in "src". This
    ///   value is placed in "dst" as well.
    /// @sa
    ///   Version of Replace() that returns a new string.
    string& Replace(const string& src,
                           const string& search,
                           const string& replace,
                           string&       dst,
                           SIZE_TYPE     start_pos = 0,
                           SIZE_TYPE     max_replace = 0,
                           SIZE_TYPE*    num_replace = 0);

    /// Replace occurrences of a substring within a string and returns the
    /// result as a new string.
    ///
    /// @param src
    ///   Source string from which specified substring occurrences are
    ///   replaced.
    /// @param search
    ///   Substring value in "src" that is replaced.
    /// @param replace
    ///   Replace "search" substring with this value.
    /// @param start_pos
    ///   Position to start search from.
    /// @param max_replace
    ///   Replace no more than "max_replace" occurrences of substring "search"
    ///   If "max_replace" is zero(default), then replace all occurrences with
    ///   "replace".
    /// @param num_replace
    ///   Optional pointer to a value which receives number of replacements occurred.
    /// @return
    ///   A new string containing the result of replacing the "search" string
    ///   with "replace" in "src"
    /// @sa
    ///   Version of Replace() that has a destination parameter to accept
    ///   result.
    string Replace(const string& src,
                          const string& search,
                          const string& replace,
                          SIZE_TYPE     start_pos = 0,
                          SIZE_TYPE     max_replace = 0,
                          SIZE_TYPE*    num_replace = 0);
                                                                                                 

END_HBNSTR_SCOPE

template<typename TIterator, typename FTransform>
std::string
NStr::TransformJoin( TIterator from, TIterator to, const ncbi::CTempString& delim, FTransform fnTransform)
{
    if (from == to) {
        return kEmptyStr;
    }
    string result(fnTransform(*from++));
    for ( ; from != to; ++from) {
        result.append(delim).append(fnTransform(*from));
    }
    return result;    
}

template<typename TIterator>
std::string
NStr::xx_Join( TIterator from, TIterator to, const ncbi::CTempString& delim)
{
    if (from == to) {
        return kEmptyStr;
    }
    string result(*from++);
    size_t sz_all = 0, sz_delim = delim.size();
    for ( TIterator f = from; f != to; ++f) {
        sz_all += string(*f).size() + sz_delim;
    }
    result.reserve(result.size() + sz_all);
    for ( ; from != to; ++from) {
        result.append(delim).append(string(*from));
    }
    return result;    
}

inline
std::list<std::string>& NStr::Wrap(const string& str, SIZE_TYPE width, list<string>& arr,
                         NStr::TWrapFlags flags, const string& prefix,
                         const string* prefix1)
{
    return Wrap(str, width, arr, flags, &prefix, prefix1);
}

inline
std::list<std::string>& NStr::Wrap(const string& str, SIZE_TYPE width, list<string>& arr,
                         NStr::TWrapFlags flags, const string& prefix,
                         const string& prefix1)
{
    return Wrap(str, width, arr, flags, &prefix, &prefix1);
}

inline
std::list<std::string>& NStr::WrapList(const list<string>& l, SIZE_TYPE width,
                             const string& delim, list<string>& arr,
                             NStr::TWrapFlags flags, const string& prefix,
                             const string* prefix1)
{
    return WrapList(l, width, delim, arr, flags, &prefix, prefix1);
}

inline
std::list<std::string>& NStr::WrapList(const list<string>& l, SIZE_TYPE width,
                             const string& delim, list<string>& arr,
                             NStr::TWrapFlags flags, const string& prefix,
                             const string& prefix1)
{
    return WrapList(l, width, delim, arr, flags, &prefix, &prefix1);
}

#endif // __STR_UTIL_HPP