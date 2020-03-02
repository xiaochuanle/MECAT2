#ifndef __NCBI_BLAST_AUX_HPP
#define __NCBI_BLAST_AUX_HPP

#include "c_ncbi_blast_aux.h"

#include <iostream>
#include <list>
#include <string>
#include <limits>
#include <memory>


/// Define the name for the NCBI namespace.
#define NCBI_NS_NCBI ncbi

/// Define a new scope.
#define BEGIN_SCOPE(ns) namespace ns {

/// End the previously defined scope.
#define END_SCOPE(ns) }

/// Use the specified namespace.
#define USING_SCOPE(ns) using namespace ns

/// Define ncbi namespace.
///
/// Place at beginning of file for NCBI related code.
#define BEGIN_NCBI_SCOPE BEGIN_SCOPE(NCBI_NS_NCBI)

/// End previously defined NCBI scope.
#define END_NCBI_SCOPE   END_SCOPE(NCBI_NS_NCBI)

#define BEGIN_HBNSTR_SCOPE BEGIN_SCOPE(ns_hbnstr)
#define END_HBNSTR_SCOPE END_SCOPE(ns_hbnstr)
#define NStr ns_hbnstr

BEGIN_NCBI_SCOPE

#define NCBI_XNCBI_EXPORT
#define NCBI_BLASTINPUT_EXPORT
#define NCBI_ALIGN_FORMAT_EXPORT

#if (SIZEOF_LONG_DOUBLE > SIZEOF_DOUBLE)
#  define NCBI_CONST_LONGDOUBLE(v)   v##L
#else
#  define NCBI_CONST_LONGDOUBLE(v)   v
#endif

#ifndef _ASSERT
#define _ASSERT hbn_assert
#endif

#define NPOS std::string::npos

#define _TROUBLE

#define ITERATE(container_type_, iter_, container_) \
    for (auto iter_ = (container_).begin(); iter_ != (container_).end(); ++iter_)

#define NON_CONST_ITERATE(container_type_, iter_, container_) ITERATE(container_type_, iter_, container_)

#define HBN_REMOVE_THIS 1

#define NCBI_DEPRECATED

#define AutoPtr unique_ptr
#define CRef shared_ptr
#define CConstRef shared_ptr

using std::list;
using std::min;
using std::max;
using std::string;

typedef std::istream CNcbiIstream;
typedef std::ostream CNcbiOstream;
typedef std::iostream CNcbiIostream;

typedef size_t SIZE_TYPE;
typedef int TIntId;
typedef int CDir;
typedef int CTime;

extern const string kEmptyStr;
#if 0
extern const int kMin_Int;
extern const int kMax_Int;
extern const unsigned int kMax_UInt;
extern const Int8 kMax_I8;
extern const Int8 kMin_I8;
extern const Uint8 kMax_UI8;
extern const double kMin_Double;
extern const double kMax_Double;
#else
#define kMin_Int (std::numeric_limits<int>::min())
#define kMax_Int (std::numeric_limits<int>::max())
#define kMax_UInt (std::numeric_limits<unsigned int>::max())
#define kMin_Long (std::numeric_limits<long>::min())
#define kMax_Long (std::numeric_limits<long>::max())
#define kMax_ULong (std::numeric_limits<unsigned long>::max())
#define kMin_I8 (std::numeric_limits<Int8>::min())
#define kMax_I8 (std::numeric_limits<Int8>::max())
#define kMax_UI8 (std::numeric_limits<Uint8>::max())
#define kMin_Double (std::numeric_limits<double>::min())
#define kMax_Double (std::numeric_limits<double>::max())
#endif

extern const size_t kDfltLineLength;
extern const size_t kDfltArgNumDescriptions;
extern const size_t kDfltArgNumAlignments;

/// Type for sequence locations and lengths.
///
/// Use this typedef rather than its expansion, which may change.
typedef unsigned int TSeqPos;

/// Define special value for invalid sequence position.
const TSeqPos kInvalidSeqPos = ((TSeqPos) (-1));

/// Low level macro for declaring safe bool operator.
#define DECLARE_SAFE_BOOL_METHOD(Expr)                                  \
    struct SSafeBoolTag {                                               \
        void SafeBoolTrue(SSafeBoolTag*) {}                             \
    };                                                                  \
    typedef void (SSafeBoolTag::*TBoolType)(SSafeBoolTag*);             \
    operator TBoolType() const {                                        \
        return (Expr)? &SSafeBoolTag::SafeBoolTrue: 0;                  \
    }                                                                   \
    private:                                                            \
    bool operator==(TBoolType) const;                                   \
    bool operator!=(TBoolType) const;                                   \
    public:                                                             \
    friend struct SSafeBoolTag


/// Declaration of safe bool operator from boolean expression.
/// Actual operator declaration will be:
///    operator TBoolType(void) const;
/// where TBoolType is a typedef convertible to bool (member pointer).
#define DECLARE_OPERATOR_BOOL(Expr)             \
    DECLARE_SAFE_BOOL_METHOD(Expr)


/// Declaration of safe bool operator from pointer expression.
/// Actual operator declaration will be:
///    operator bool(void) const;
#define DECLARE_OPERATOR_BOOL_PTR(Ptr)          \
    DECLARE_OPERATOR_BOOL((Ptr) != 0)

/// Template used to replace bool type arguments with some strict equivalent.
/// This allow to prevent compiler to do an implicit casts from other types
/// to bool. "TEnum" should be an enumerated type with two values:
///   - negative (FALSE/OFF) should have value 0
//    - positive (TRUE/ON) have value 1.

template <class TEnum>
class CBoolEnum
{
public:
    // Constructors
    CBoolEnum(bool  value) : m_Value(value) {}
    CBoolEnum(TEnum value) : m_Value(value ? true : false) {}

    /// Operator bool
    operator bool() const   { return m_Value; }
    /// Operator enum
    operator TEnum () const { return TEnum(m_Value); }

private:
    bool m_Value;

private:
    // Disable implicit conversions from/to other types
    CBoolEnum(char*);
    CBoolEnum(int);
    CBoolEnum(unsigned int);
    template <class T> CBoolEnum(T);

    operator int() const;
    operator unsigned int() const;
};

    /// Flags for Split*() methods. 
    /// 
    /// @note
    ///   With quote support enabled, doubling a quote character suppresses
    ///   its special meaning, as does escaping it if that's enabled too;
    ///   unescaped trailing backslashes and unbalanced quotes result in
    ///   exceptions.
    /// @note
    ///   All escape symbols, single or double quotes became removed
    ///   if a corresponding fSplit_Can* flag is used.
    enum ESplitFlags {
        fSplit_MergeDelimiters = 1 << 0,  ///< Merge adjacent delimiters
        fSplit_Truncate_Begin  = 1 << 1,  ///< Truncate leading delimiters
        fSplit_Truncate_End    = 1 << 2,  ///< Truncate trailing delimiters
        fSplit_Truncate        = fSplit_Truncate_Begin | fSplit_Truncate_End,
        fSplit_ByPattern       = 1 << 3,  ///< Require full delimiter strings
        fSplit_CanEscape       = 1 << 4,  ///< Allow \\... escaping
        fSplit_CanSingleQuote  = 1 << 5,  ///< Allow '...' quoting
        fSplit_CanDoubleQuote  = 1 << 6,  ///< Allow "..." quoting
        fSplit_CanQuote        = fSplit_CanSingleQuote | fSplit_CanDoubleQuote,
        /// All delimiters are merged and trimmed, to get non-empty tokens only
        fSplit_Tokenize        = fSplit_MergeDelimiters | fSplit_Truncate
    };
    typedef int TSplitFlags; ///< Bitwise OR of ESplitFlags

    /// Which type of string comparison.
    enum class ECase {
        eCase,      ///< Case sensitive compare
        eNocase     ///< Case insensitive compare
    };

bool VerifyArgumentName(const string& name, bool extended = false);

struct CObject {};

END_NCBI_SCOPE

#endif // __NCBI_BLAST_AUX_HPP
