#ifndef __NCBIARGS_ALLOW_HPP
#define __NCBIARGS_ALLOW_HPP

#include "ncbiargs_types.hpp"
#include "../str_util/ncbistr.hpp"

#include <list>
#include <map>
#include <set>
#include <sstream>

BEGIN_NCBI_SCOPE

using namespace NStr;

using std::list;
using std::map;
using std::multimap;
using std::set;
using std::pair;

/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow --
///
/// Abstract base class for defining argument constraints.
///
/// Other user defined constraints are defined by deriving from this abstract
/// base class:
///
///  - CArgAllow_Symbols  -- symbol from a set of allowed symbols
///  - CArgAllow_String   -- string to contain only allowed symbols 
///  - CArgAllow_Strings  -- string from a set of allowed strings
///  - CArgAllow_Int8s    -- Int8    value to fall within a given interval
///  - CArgAllow_Integers -- integer value to fall within a given interval
///  - CArgAllow_Doubles  -- floating-point value to fall in a given interval
///
/// @sa CArgAllow_Regexp

class NCBI_XNCBI_EXPORT CArgAllow
{
public:
    /// Verify if specified value is allowed.
    virtual bool Verify(const string &value) const = 0;

    /// Get usage information.
    virtual
    string GetUsage(void) const = 0;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    virtual ~CArgAllow(void);

    /// Create object's clone, moving it from stack memory into heap.
    /// The method is required only when using objects created on stack.
    ///
    /// NOTE: Default implementation returns NULL.
    /// Inherited classes must override this method.
    ///
    /// @sa CArgDescriptions::SetConstraint
    virtual CArgAllow* Clone(void) const;

protected:
    // In the absence of the following constructor, new compilers (as required
    // by the new C++ standard) may fill the object memory with zeros,
    // erasing flags set by CObject::operator new (see CXX-1808)
    CArgAllow(void) {}
};



/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_Symbols --
///
/// Define constraint to describe exactly one symbol.
///
/// Argument to be exactly one symbol from the specified set of symbols.
///
/// Examples:
/// - To allow only symbols 'a', 'b' and 'Z' for argument "MyArg":
///   SetConstraint("MyArg", new CArgAllow_Symbols("abZ"))
/// - To allow only printable symbols (according to "isprint()" from <ctype.h>):
///   SetConstraint("MyArg", new CArgAllow_Symbols(CArgAllow_Symbols::ePrint))

class NCBI_XNCBI_EXPORT CArgAllow_Symbols : public CArgAllow
{
public:
    /// Symbol class for defining sets of characters.
    ///
    /// Symbol character classes patterned after those defined in <ctype.h>.
    enum ESymbolClass {
        // Standard character class from <ctype.h>:  isalpha(), isdigit(), etc.
        eAlnum,  ///< Alphanumeric characters
        eAlpha,  ///< Alphabet characters
        eCntrl,  ///< Control characters
        eDigit,  ///< Digit characters
        eGraph,  ///< Graphical characters
        eLower,  ///< Lowercase characters
        ePrint,  ///< Printable characters
        ePunct,  ///< Punctuation characters
        eSpace,  ///< Space characters
        eUpper,  ///< Uppercase characters
        eXdigit, ///< Hexadecimal characters
        eUser    ///< User defined characters using constructor with string&
    };

    /// Constructor.
    CArgAllow_Symbols(ESymbolClass symbol_class);

    /// Constructor for user defined eUser class.
    CArgAllow_Symbols(const string& symbol_set);

    /// Add allowed symbols
    CArgAllow_Symbols& Allow(ESymbolClass  symbol_class);
    CArgAllow_Symbols& Allow(const string& symbol_set);

protected:
    /// Verify if specified value is allowed.
    virtual bool Verify(const string& value) const;

    /// Get usage information.
    virtual string GetUsage(void) const;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    CArgAllow_Symbols(void) {
    }
    virtual CArgAllow* Clone(void) const;

    typedef pair<ESymbolClass, string> TSymClass;
    set< TSymClass > m_SymClass;
};


/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_String --
///
/// Define constraint to describe string argument.
///
/// Argument to be a string containing only allowed symbols.
///
/// Examples:
/// - To allow string containing only symbols 'a', 'b' and 'Z' for arg MyArg:
///   SetConstraint("MyArg", new CArgAllow_String("abZ"))
/// - To allow only numeric symbols (according to "isdigit()" from <ctype.h>):
///   SetConstraint("MyArg", new CArgAllow_String(CArgAllow_String::eDigit))

class NCBI_XNCBI_EXPORT CArgAllow_String : public CArgAllow_Symbols
{
public:
    /// Constructor.
    CArgAllow_String(ESymbolClass symbol_class);

    /// Constructor for user defined eUser class.
    CArgAllow_String(const string& symbol_set);

protected:
    /// Verify if specified value is allowed.
    virtual bool Verify(const string& value) const;

    /// Get usage information.
    virtual string GetUsage(void) const;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    CArgAllow_String(void) {
    }
    virtual CArgAllow* Clone(void) const;
};



/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_Strings --
///
/// Define constraint to describe set of string values.
///
/// Argument to have only particular string values. Use the Allow() method to
/// add the allowed string values, which can be daisy-chained.
///
/// Examples:
/// - SetConstraint("a", (new CArgAllow_Strings)->
///                  Allow("foo")->Allow("bar")->Allow("etc"))
/// - You can use "operator,()" to shorten the notation:
///   SetConstraint("b", &(*new CArgAllow_Strings, "foo", "bar", "etc"))

class NCBI_XNCBI_EXPORT CArgAllow_Strings : public CArgAllow
{
public:
    /// Constructor
    /// @param use_case
    ///   If to ignore the case of the characters
    CArgAllow_Strings(NStr::ECase use_case = NStr::ECase::eCase);

    /// Add allowed string values
    /// @param value
    ///   String to add to the set of allowed string values
    CArgAllow_Strings* Allow(const string& value);

    /// Add allowed string values
    /// @param value
    ///   String to add to the set of allowed string values
    CArgAllow_Strings& AllowValue(const string& value);

    /// Short notation operator for adding allowed string values
    /// @param value
    ///   String to add to the set of allowed string values
    /// @sa
    ///   Allow()
    CArgAllow_Strings& operator,(const string& value) {
        return AllowValue(value);
    }

protected:
    /// Verify if specified value is allowed.
    virtual bool Verify(const string& value) const;

    /// Get usage information.
    virtual string GetUsage(void) const;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    virtual CArgAllow* Clone(void) const;

    /// Type of the container that contains the allowed string values
    /// @sa m_Strings
    typedef set<string, PNocase_Conditional> TStrings;

    TStrings     m_Strings;  ///< Set of allowed string values
};


/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_Int8s --
///
/// Define constraint to describe range of 8-byte integer values and TIntIds.
///
/// Argument to have only integer values falling within given interval.
///
/// Example:
/// - SetConstraint("a2", new CArgAllow_Int8s(-1001, 123456789012))

class NCBI_XNCBI_EXPORT CArgAllow_Int8s : public CArgAllow
{
public:
    /// Constructor specifying an allowed integer value.
    CArgAllow_Int8s(Int8 x_value);

    /// Constructor specifying range of allowed integer values.
    CArgAllow_Int8s(Int8 x_min, Int8 x_max);

    /// Add allow values
    CArgAllow_Int8s& AllowRange(Int8 from, Int8 to);
    CArgAllow_Int8s& Allow(Int8 value);

protected:
    /// Verify if specified value is allowed.
    virtual bool   Verify(const string& value) const;

    /// Get usage information.
    virtual string GetUsage(void) const;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    CArgAllow_Int8s(void) {
    }
    virtual CArgAllow* Clone(void) const;


    typedef pair<Int8, Int8> TInterval;
    set< TInterval >  m_MinMax;
};



/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_Integers --
///
/// Define constraint to describe range of integer id values.
///
/// Argument to have only integer values falling within given interval.
///
/// Example:
/// - SetConstraint("i", new CArgAllow_Integers(kMin_Int, 34))

class NCBI_XNCBI_EXPORT CArgAllow_Integers : public CArgAllow_Int8s
{
public:
    /// Constructor specifying an allowed integer value.
    CArgAllow_Integers(int x_value);
    /// Constructor specifying range of allowed integer values.
    CArgAllow_Integers(int x_min, int x_max);

protected:
    /// Get usage information.
    virtual string GetUsage(void) const;

    CArgAllow_Integers(void) {
    }
    virtual CArgAllow* Clone(void) const;
};



/////////////////////////////////////////////////////////////////////////////
///
/// CArgAllow_Doubles --
///
/// Define constraint to describe range of double values.
///
/// Argument to have only double values falling within given interval.
///
/// Example:
/// - SetConstraint("d", new CArgAllow_Doubles(0.01, 0.99))

class NCBI_XNCBI_EXPORT CArgAllow_Doubles : public CArgAllow
{
public:
    /// Constructor specifying an allowed double value.
    CArgAllow_Doubles(double x_value);

    /// Constructor specifying range of allowed double values.
    CArgAllow_Doubles(double x_min, double x_max);

    /// Add allowed values
    CArgAllow_Doubles& AllowRange(double from, double to);
    CArgAllow_Doubles& Allow(double value);

protected:
    /// Verify if specified value is allowed.
    virtual bool   Verify(const string& value) const;

    /// Get usage information.
    virtual string GetUsage(void) const;

    /// Print constraints in XML format
    virtual void PrintUsageXml(CNcbiOstream& out) const;

    CArgAllow_Doubles(void) {
    }
    virtual CArgAllow* Clone(void) const;

    typedef pair<double,double> TInterval;
    set< TInterval >  m_MinMax;
};

/// Class to constrain the values of an argument to those greater than or equal
/// to the value specified in the constructor
class NCBI_BLASTINPUT_EXPORT CArgAllowValuesGreaterThanOrEqual : public CArgAllow
{
public:
    /// Constructor taking an integer
    CArgAllowValuesGreaterThanOrEqual(int min) : m_MinValue(min) {}
    /// Constructor taking a double
    CArgAllowValuesGreaterThanOrEqual(double min) : m_MinValue(min) {}

    virtual CArgAllowValuesGreaterThanOrEqual* Clone() const {
        return new CArgAllowValuesGreaterThanOrEqual(m_MinValue);
    }

protected:
    /// Overloaded method from CArgAllow
    virtual bool Verify(const string& value) const {
        return NStr::StringToDouble(value) >= m_MinValue;
    }

    /// Overloaded method from CArgAllow
    virtual string GetUsage(void) const {
        return ">=" + NStr::DoubleToString(m_MinValue);
    }
    
private:
    double m_MinValue;  /**< Minimum value for this object */
};

/// Class to constrain the values of an argument to those less than or equal
/// to the value specified in the constructor
class NCBI_BLASTINPUT_EXPORT CArgAllowValuesLessThanOrEqual : public CArgAllow
{
public:
    /// Constructor taking an integer
    CArgAllowValuesLessThanOrEqual(int max) : m_MaxValue(max) {}
    /// Constructor taking a double
    CArgAllowValuesLessThanOrEqual(double max) : m_MaxValue(max) {}

protected:
    /// Overloaded method from CArgAllow
    virtual bool Verify(const string& value) const {
        return NStr::StringToDouble(value) <= m_MaxValue;
    }

    /// Overloaded method from CArgAllow
    virtual string GetUsage(void) const {
        return "<=" + NStr::DoubleToString(m_MaxValue);
    }
    
private:
    double m_MaxValue;  /**< Maximum value for this object */
};

/// Class to constrain the values of an argument to those in between the values
/// specified in the constructor
class NCBI_BLASTINPUT_EXPORT CArgAllowValuesBetween : public CArgAllow
{
public:
    /// Constructor taking an integer
    CArgAllowValuesBetween(int min, int max, bool inclusive = false) 
        : m_MinValue(min), m_MaxValue(max), m_Inclusive(inclusive) {}
    /// Constructor taking a double
    CArgAllowValuesBetween(double min, double max, bool inclusive = false)
        : m_MinValue(min), m_MaxValue(max), m_Inclusive(inclusive) {}

    virtual CArgAllowValuesBetween* Clone() const {
        return new CArgAllowValuesBetween(m_MinValue, m_MaxValue, m_Inclusive);
    }

protected:
    /// Overloaded method from CArgAllow
    virtual bool Verify(const string& value) const {
        double val = NStr::StringToDouble(value);
        bool retval = false;
        if ( !m_Inclusive ) {
            retval = (val > m_MinValue && val < m_MaxValue);
        } else {
            retval = (val >= m_MinValue && val <= m_MaxValue);
        }
        return retval;
    }

    /// Overloaded method from CArgAllow
    virtual string GetUsage(void) const {
        string retval;
        if ( !m_Inclusive ) {
            retval = "(>" + NStr::DoubleToString(m_MinValue) + " and <"
                + NStr::DoubleToString(m_MaxValue) + ")";
        } else {
            retval = "(>=" + NStr::DoubleToString(m_MinValue) + " and =<"
                + NStr::DoubleToString(m_MaxValue) + ")";
        }
        return retval;
    }
    
private:
    double m_MinValue;  /**< Minimum value for this object */
    double m_MaxValue;  /**< Maximum value for this object */
    bool m_Inclusive;   /**< Whether the values above should be included or not */
};

/** 
 * @brief Macro to create a subclass of CArgAllow that allows the specification
 * of sets of data
 * 
 * @param ClassName Name of the class to be created [in]
 * @param DataType data type of the allowed arguments [in]
 * @param String2DataTypeFn Conversion function from a string to DataType [in]
 */
#define DEFINE_CARGALLOW_SET_CLASS(ClassName, DataType, String2DataTypeFn)  \
class NCBI_BLASTINPUT_EXPORT ClassName : public CArgAllow                       \
{                                                                           \
public:                                                                     \
    ClassName(const set<DataType>& values)                                  \
        : m_AllowedValues(values)                                           \
    {                                                                       \
        if (values.empty()) {                                               \
            throw runtime_error("Allowed values set must not be empty");    \
        }                                                                   \
    }                                                                       \
                                                                            \
protected:                                                                  \
    virtual bool Verify(const string& value) const {                        \
        DataType value2check = String2DataTypeFn(value);                    \
        ITERATE(set<DataType>, itr, m_AllowedValues) {                      \
            if (*itr == value2check) {                                      \
                return true;                                                \
            }                                                               \
        }                                                                   \
        return false;                                                       \
    }                                                                       \
                                                                            \
    virtual string GetUsage(void) const {                                   \
        std::ostringstream os;                                                 \
        os << "Permissible values: ";                                       \
        ITERATE(set<DataType>, itr, m_AllowedValues) {                      \
            os << "'" << *itr << "' ";                                      \
        }                                                                   \
        return os.str();                                 \
    }                                                                       \
                                                                            \
private:                                                                    \
    /* Set containing the permissible values */                             \
    set<DataType> m_AllowedValues;                                          \
}

#ifndef SKIP_DOXYGEN_PROCESSING
DEFINE_CARGALLOW_SET_CLASS(CArgAllowIntegerSet, int, NStr::StringToInt);
DEFINE_CARGALLOW_SET_CLASS(CArgAllowStringSet, string, string);
#endif /* SKIP_DOXYGEN_PROCESSING */

END_NCBI_SCOPE

#endif // __NCBIARGS_ALLOW_HPP