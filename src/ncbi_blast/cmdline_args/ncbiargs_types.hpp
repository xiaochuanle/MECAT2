#ifndef __NCBIARG_TYPES_HPP
#define __NCBIARG_TYPES_HPP

/// @file ncbiargs.hpp
/// Defines command line argument related classes.
///
/// The CArgDescriptions and CArgDesc classes are used for describing
/// unparsed arguments; CArgs and CArgValue for parsed argument values;
/// CArgException and CArgHelpException for argument exceptions; and CArgAllow, 
/// CArgAllow_{Strings, ..., Integers, Doubles} for argument constraints.
///
/// The following description is included as applies to several classes in
/// this file:
///
/// Parsing and validation of command-line arguments are done according to
/// user-provided descriptions. The command line has the following syntax:
///
/// Command string:
///
///    progname  {arg_key, arg_key_opt, arg_key_dflt, arg_flag} [--]
///              {arg_pos} {arg_pos_opt, arg_pos_dflt}
///              {arg_extra} {arg_extra_opt}
///
/// where:
///
///   arg_key        :=  -<key> <value>    -- (mandatory)
///   arg_key_opt    := [-<key> <value>]   -- (optional, without default value)
///   arg_key_dflt   := [-<key> <value>]   -- (optional, with default value)
///   arg_flag       := -<flag>            -- (always optional)
///   "--" is an optional delimiter to indicate the beginning of pos. args
///   arg_pos        := <value>            -- (mandatory)
///   arg_pos_opt    := [<value>]          -- (optional, without default value)
///   arg_pos_dflt   := [<value>]          -- (optional, with default value)
///   arg_extra      := <value>            -- (dep. on the constraint policy)
///   arg_extra_opt  := [<value>]          -- (dep. on the constraint policy)
///
/// and:
///
///   <key> must be followed by <value>
///   <flag> and <key> are case-sensitive, and they can contain
///                    only alphanumeric characters
///   <value> is an arbitrary string (additional constraints can
///           be applied in the argument description, see "EType")
///
/// {arg_pos***} and {arg_extra***} -- position-dependent arguments, with
/// no tag preceding them.
/// {arg_pos***} -- have individual names and descriptions (see methods
/// AddPositional***).
/// {arg_extra***} have one description for all (see method AddExtra).
/// User can apply constraints on the number of mandatory and optional
/// {arg_extra***} arguments.

#include "../ncbi_blast_aux.hpp"
#include "../../corelib/hbn_aux.h"

#include <vector>

BEGIN_NCBI_SCOPE

using namespace std;

/////////////////////////////////////////////////////////////////////////////
///
/// CArgValue --
///
/// Generic abstract base class for argument values.

class NCBI_XNCBI_EXPORT CArgValue
{
public:
    /// Get argument name.
    const string& GetName(void) const { return m_Name; }

    /// Check if argument holds a value.
    ///
    /// Argument does not hold value if it was described as optional argument
    /// without default value, and if it was not passed a value in the command
    /// line.  On attempt to retrieve the value from such "no-value" argument,
    /// exception will be thrown.
    virtual bool HasValue(void) const = 0;

    virtual int ArgValueExpected(void) const = 0;

    DECLARE_OPERATOR_BOOL(HasValue());

    /// Get the argument's string value.
    ///
    /// If it is a value of a flag argument, then return either "true"
    /// or "false".
    /// @sa
    ///   AsInteger(), AsInt8(), AsDouble(), AsBoolean()
    virtual const string& AsString(void) const = 0;

    /// Get the argument's integer (8-byte long) value.
    ///
    /// If you request a wrong value type, such as a call to "AsInt8()"
    /// for a "boolean" argument, an exception is thrown.
    /// This will however work okay for "plain integer" argument.
    /// @sa
    ///   AsInteger(), AsString(), AsDouble, AsBoolean()
    virtual Int8 AsInt8(void) const = 0;

    /// Get the argument's integer value.
    ///
    /// If you request a wrong value type, such as a call to "AsInteger()"
    /// for a "boolean" or even "Int8" argument, an exception is thrown.
    /// @sa
    ///   AsInt8(), AsString(), AsDouble, AsBoolean()
    virtual int    AsInteger(void) const = 0;

    /// Get the argument's value as an integer id (TIntId). The actual value is
    /// Int4 or Int8 depending on the NCBI_INT8_GI definition.
    ///
    /// If you request a wrong value type, such as a call to "AsIntId()"
    /// for a "boolean", an exception is thrown. Calling AsIntId() on an
    /// integer argument is always allowed. For an Int8 argument it will
    /// throw an exception if NCBI_INT8_GI is not defined.
    /// @sa
    ///   AsInteger(), AsInt8()
    virtual TIntId AsIntId(void) const = 0;

    /// Get the argument's double value.
    ///
    /// If you request a wrong value type, such as a call to "AsDouble()"
    /// for a "boolean" argument, an exception is thrown.
    /// @sa
    ///   AsString(), AsInt8(), AsInteger, AsBoolean()
    virtual double AsDouble (void) const = 0;

    /// Get the argument's boolean value.
    ///
    /// If you request a wrong value type, such as a call to "AsBoolean()"
    /// for a "integer" argument, an exception is thrown.
    /// @sa
    ///   AsString(), AsInt8(), AsInteger, AsDouble()
    virtual bool   AsBoolean(void) const = 0;

    enum EFileFlags {
        fBinary   = (1 <<  1),  ///< Open file in binary mode.
        fText     =  0,         ///< Open file in text mode.
        fAppend   = (1 <<  2),  ///< Open file in append mode.
        fTruncate = (1 << 12),  ///< Open file in truncate mode.
        fNoCreate = (1 << 11),  ///< Open existing file, never create it
        fCreatePath = (1 << 8)  ///< If needed, create directory where the file is located
    };
    typedef unsigned int TFileFlags;   ///< Bitwise OR of "EFileFlags"

    /// Get the argument as an input file stream.
    virtual CNcbiIstream& AsInputFile (TFileFlags flags = 0) const = 0;

    /// Get the argument as an output file stream.
    virtual CNcbiOstream& AsOutputFile(TFileFlags flags = 0) const = 0;

    /// Get the argument as a file stream.
    virtual CNcbiIostream& AsIOFile(TFileFlags flags = 0) const = 0;

    /// Get the argument as a directory.
    virtual const CDir& AsDirectory(void) const = 0;

    /// Get the argument as a DateTime.
    virtual const CTime& AsDateTime(void) const = 0;

    /// Close the file.
    virtual void CloseFile (void) const = 0;



    /// Some values types can contain several value lists
    ///
    /// Example: CGIs pass list selections by repeating the same name
    typedef vector<string>  TStringArray;

    /// Get the value list
    virtual const TStringArray& GetStringList() const;

    /// Get reference on value list for further modification
    virtual TStringArray& SetStringList();
    
    /// Get ordinal position of the value.
    /// NOTE: this is not the position in command line, rather
    /// this reflects the order in which values were added to the list.
    size_t GetOrdinalPosition(void) const
    {
        return m_Ordinal;
    }

    /// Whether the argument:
    /// @sa GetDefault()
    enum EArgValueFlags {
        fArgValue_HasDefault  = (1 << 0),  ///< Has default value
        fArgValue_FromDefault = (1 << 1)   ///< Not provided in command line
    };
    typedef unsigned int TArgValueFlags;  ///< Bitwise OR of "EArgValueFlags"

    /// Get default value of the argument.
    ///
    /// @param flags
    ///   Indicate whether the argument has default value, and if the arg's
    ///   value was set from the command line or from the default value.
    /// @return
    ///   Default value, if specified for this argument.
    ///   If the argument doesn't have a default value: empty string.
    ///   If the argument is a flag: "false" or "true".
    const string& GetDefault(TArgValueFlags* flags = NULL) const;

    virtual ~CArgValue(void);

protected:
    friend class CArgs;
    friend class CArgDescDefault;
    friend class CArgDescMandatory;
    friend class CArgDesc_Flag;

    /// Protected constructor and destructor.
    ///
    /// Prohibit explicit instantiation of CArgValue with name.
    CArgValue(const string& name);
    
    void SetOrdinalPosition(size_t pos)
    {
        m_Ordinal = pos;
    }
    void x_SetDefault(const string& def_value, bool from_def);

    string m_Name;          ///< Argument name
    size_t m_Ordinal;
    string m_Default;
    TArgValueFlags m_Flags;
};

//  Overload the comparison operator -- to handle "CRef<CArgValue>" elements
//  in "CArgs::m_Args" stored as "set< CRef<CArgValue> >"
//
inline bool operator< (const CRef<CArgValue>& x, const CRef<CArgValue>& y)
{
    return x->GetName() < y->GetName();
}

struct CArgValueCRefLessThan
{
    bool operator()(const CRef<CArgValue>& x, const CRef<CArgValue>& y) const {
        return x < y;
    }
};

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//  CArg_***::   classes representing various types of argument value
//
//    CArg_NoValue     : CArgValue
//
//    CArg_String      : CArgValue
//
//       CArg_Alnum       : CArg_String
//       CArg_Int8        : CArg_String
//          CArg_Integer  : CArg_Int8
//          CArg_IntIds   : CArg_Int8
//       CArg_DataSize    : CArg_String
//       CArg_Double      : CArg_String
//       CArg_Boolean     : CArg_String
//       CArg_InputFile   : CArg_String
//       CArg_OutputFile  : CArg_String
//       CArg_IOFile      : CArg_String
//       CArg_Dir         : CArg_String
//       CArg_DateTime    : CArg_String
//    


class CArg_NoValue : public CArgValue
{
public:
    CArg_NoValue(const string& name);
    virtual bool HasValue(void) const;

    virtual int ArgValueExpected(void) const { return 0; }

    virtual const string&  AsString (void) const;
    virtual Int8           AsInt8   (void) const;
    virtual int            AsInteger(void) const;
    virtual TIntId         AsIntId  (void) const;
    virtual double         AsDouble (void) const;
    virtual bool           AsBoolean(void) const;
    virtual const CDir&    AsDirectory(void) const;
    virtual const CTime&   AsDateTime(void) const;

    virtual CNcbiIstream&  AsInputFile (TFileFlags flags = 0) const;
    virtual CNcbiOstream&  AsOutputFile(TFileFlags flags = 0) const;
    virtual CNcbiIostream& AsIOFile(TFileFlags flags = 0) const;
    virtual void           CloseFile   (void) const;
};


// Generates error (like CArg_NoValue) for excluded arguments.
class CArg_ExcludedValue : public CArgValue
{
public:
    CArg_ExcludedValue(const string& name);
    virtual bool HasValue(void) const;

    virtual int ArgValueExpected(void) const { return 0; }

    virtual const string&  AsString (void) const;
    virtual Int8           AsInt8   (void) const;
    virtual int            AsInteger(void) const;
    virtual TIntId         AsIntId  (void) const;
    virtual double         AsDouble (void) const;
    virtual bool           AsBoolean(void) const;
    virtual const CDir&    AsDirectory(void) const;
    virtual const CTime&   AsDateTime(void) const;

    virtual CNcbiIstream&  AsInputFile (TFileFlags flags = 0) const;
    virtual CNcbiOstream&  AsOutputFile(TFileFlags flags = 0) const;
    virtual CNcbiIostream& AsIOFile(TFileFlags flags = 0) const;
    virtual void           CloseFile   (void) const;
};


class CArg_String : public CArgValue
{
public:
    CArg_String(const string& name, const string& value);
    virtual bool HasValue(void) const;

    virtual int ArgValueExpected(void) const { return 1; }

    virtual const string&  AsString (void) const;
    virtual Int8           AsInt8   (void) const;
    virtual int            AsInteger(void) const;
    virtual TIntId         AsIntId  (void) const;
    virtual double         AsDouble (void) const;
    virtual bool           AsBoolean(void) const;
    virtual const CDir&    AsDirectory(void) const;
    virtual const CTime&   AsDateTime(void) const;

    virtual CNcbiIstream&  AsInputFile (TFileFlags flags = 0) const;
    virtual CNcbiOstream&  AsOutputFile(TFileFlags flags = 0) const;
    virtual CNcbiIostream& AsIOFile(TFileFlags flags = 0) const;
    virtual void           CloseFile   (void) const;

    virtual const TStringArray& GetStringList() const;
    virtual TStringArray& SetStringList();

private:
    /// Value of the argument as passed to the constructor ("value")
    /// becomes the first element in the value list
    /// AsString() and other methods then use it 
    TStringArray  m_StringList;
};



class CArg_Int8 : public CArg_String
{
public:
    CArg_Int8(const string& name, const string& value);
    virtual Int8 AsInt8(void) const;
    /// An Int8 argument can be used as an integer id only if NCBI_INT8_GI is defined.
    virtual TIntId AsIntId(void) const;
protected:
    Int8 m_Integer;
};



class CArg_Integer : public CArg_Int8
{
public:
    CArg_Integer(const string& name, const string& value);
    virtual int AsInteger(void) const;
    /// An integer argument can also be used as an integer id.
    virtual TIntId AsIntId(void) const;
};


class CArg_IntId : public CArg_Int8
{
public:
    CArg_IntId(const string& name, const string& value);
    /// An IntId argument can be used as an integer only if NCBI_INT8_GI is not defined.
    /// Otherwise an exception will be thrown.
    virtual int AsInteger(void) const;
    virtual TIntId AsIntId(void) const;
};


class CArg_DataSize : public CArg_String
{
public:
    CArg_DataSize(const string& name, const string& value);
    virtual Int8 AsInt8(void) const;
protected:
    Uint8 m_Integer;
};



class CArg_Double : public CArg_String
{
public:
    CArg_Double(const string& name, const string& value);
    virtual double AsDouble(void) const;
private:
    double m_Double;
};



class CArg_Boolean : public CArg_String
{
public:
    CArg_Boolean(const string& name, bool value);
    CArg_Boolean(const string& name, const string& value);
    virtual bool AsBoolean(void) const;
private:
    bool m_Boolean;
};



class CArg_Flag : public CArg_Boolean
{
public:
    CArg_Flag(const string& name, bool value);
    virtual bool HasValue(void) const;
    virtual int ArgValueExpected(void) const { return 0; }
};

class CArg_Grid : public CArg_String
{
public:
    CArg_Grid(const string& name, const string& node_id, const string& num_nodes);
    CArg_Grid(const string& name, const int node_id, const int num_nodes);
    virtual bool HasValue(void) const;

    virtual int ArgValueExpected(void) const { return 2; }

    int GetNodeId() const;
    void SetNodeId(const int node_id);
    int GetNumNodes() const;
    void SetNumNodes(const int num_nodes);

    static void Validate(const int node_id, const int num_nodes);

private:

    int m_NodeId;
    int m_NumNodes;
};

END_NCBI_SCOPE

#endif // __NCBIARG_TYPES_HPP