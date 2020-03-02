#include "ncbiargs_types.hpp"

#include "../str_util/ncbistr.hpp"

BEGIN_NCBI_SCOPE

// Allow autodetection among decimal and hex, but NOT octal, in case
// anyone's been relying on leading zeros being meaningless.
Int8 s_StringToInt8(const string& value)
{
    NStr::EHbnStrConvErrorMode orig_mode = NStr::GetStrConvErrorMode();
    NStr::SetStrConvErrorMode(NStr::EHbnStrConvErrorMode::eHbnStrConvErrorSetErrno);
    errno = 0;

    Int8 n = NStr::StringToInt8(value);
    if (n == 0 && errno != 0) {
        if (NStr::StartsWith(value, "0x", NStr::ECase::eNocase)) {
            n = 0;
            errno = 0;
            n = NStr::StringToInt8(value, 0, 16);
        }
    }
    if (n == 0 && errno != 0) HBN_ERR("Cannot convert string '%s' to Int8", value.c_str());
    NStr::SetStrConvErrorMode(orig_mode);
    return n;
}

/////////////////////////////////////////////////////////////////////////////
//  CArg_***::   classes representing various types of argument value
//
//    CArgValue
//       CArg_NoValue        : CArgValue
//       CArg_String         : CArgValue
//          CArg_Int8        : CArg_String
//             CArg_Integer  : CArg_Int8
//             CArg_IntId    : CArg_IntId
//          CArg_Double      : CArg_String
//          CArg_Boolean     : CArg_String
//          CArg_InputFile   : CArg_String
//          CArg_OutputFile  : CArg_String
//          CArg_IOFile      : CArg_String
//


///////////////////////////////////////////////////////
//  CArgValue::

CArgValue::CArgValue(const string& name)
    : m_Name(name), m_Ordinal(0), m_Flags(0)
{

    if (!VerifyArgumentName(m_Name)) HBN_ERR("Invalid argument name '%s'", m_Name.c_str());
}


CArgValue::~CArgValue(void)
{
    return;
}

const string& CArgValue::GetDefault(TArgValueFlags* has_default) const
{
    if (has_default) {
        *has_default = m_Flags;
    }
    return m_Default;
}

void CArgValue::x_SetDefault(const string& def_value, bool from_def)
{
    m_Default = def_value;
    m_Flags |= fArgValue_HasDefault;
    if (from_def) {
        m_Flags |= fArgValue_FromDefault;
    }
}

const CArgValue::TStringArray& CArgValue::GetStringList() const
{
    //NCBI_THROW(CArgException, eInvalidArg,
    string err_msg = (
        "Value lists not implemented for this argument: " + m_Name);
    HBN_ERR("%s", err_msg.c_str());
}


CArgValue::TStringArray& CArgValue::SetStringList()
{
    //NCBI_THROW(CArgException, eInvalidArg,
    string err_msg = (
        "Value lists not implemented for this argument: " + m_Name);
    HBN_ERR("%s", err_msg.c_str());
}

///////////////////////////////////////////////////////
//  CArg_NoValue::

CArg_NoValue::CArg_NoValue(const string& name)
    : CArgValue(name)
{
    return;
}


bool CArg_NoValue::HasValue(void) const
{
    return false;
}


#define THROW_CArg_NoValue HBN_ERR("The argument '%s' has no value", GetName().c_str())

const string& CArg_NoValue::AsString    (void) const { THROW_CArg_NoValue; }
Int8          CArg_NoValue::AsInt8      (void) const { THROW_CArg_NoValue; }
int           CArg_NoValue::AsInteger   (void) const { THROW_CArg_NoValue; }
TIntId        CArg_NoValue::AsIntId     (void) const { THROW_CArg_NoValue; }
double        CArg_NoValue::AsDouble    (void) const { THROW_CArg_NoValue; }
bool          CArg_NoValue::AsBoolean   (void) const { THROW_CArg_NoValue; }
const CDir&   CArg_NoValue::AsDirectory (void) const { THROW_CArg_NoValue; }
const CTime&  CArg_NoValue::AsDateTime  (void) const { THROW_CArg_NoValue; }
CNcbiIstream& CArg_NoValue::AsInputFile (CArgValue::TFileFlags) const { THROW_CArg_NoValue; }
CNcbiOstream& CArg_NoValue::AsOutputFile(CArgValue::TFileFlags) const { THROW_CArg_NoValue; }
CNcbiIostream& CArg_NoValue::AsIOFile(CArgValue::TFileFlags) const { THROW_CArg_NoValue; }
void          CArg_NoValue::CloseFile   (void) const { THROW_CArg_NoValue; }


///////////////////////////////////////////////////////
//  CArg_ExcludedValue::

CArg_ExcludedValue::CArg_ExcludedValue(const string& name)
    : CArgValue(name)
{
    return;
}


bool CArg_ExcludedValue::HasValue(void) const
{
    return false;
}


#define THROW_CArg_ExcludedValue \
    HBN_ERR("The value of argument '%s' is excluded by other arguments", GetName().c_str())

const string& CArg_ExcludedValue::AsString    (void) const { THROW_CArg_ExcludedValue; }
Int8          CArg_ExcludedValue::AsInt8      (void) const { THROW_CArg_ExcludedValue; }
int           CArg_ExcludedValue::AsInteger   (void) const { THROW_CArg_ExcludedValue; }
TIntId        CArg_ExcludedValue::AsIntId     (void) const { THROW_CArg_ExcludedValue; }
double        CArg_ExcludedValue::AsDouble    (void) const { THROW_CArg_ExcludedValue; }
bool          CArg_ExcludedValue::AsBoolean   (void) const { THROW_CArg_ExcludedValue; }
const CDir&   CArg_ExcludedValue::AsDirectory (void) const { THROW_CArg_ExcludedValue; }
const CTime&  CArg_ExcludedValue::AsDateTime  (void) const { THROW_CArg_ExcludedValue; }
CNcbiIstream& CArg_ExcludedValue::AsInputFile (CArgValue::TFileFlags) const { THROW_CArg_ExcludedValue; }
CNcbiOstream& CArg_ExcludedValue::AsOutputFile(CArgValue::TFileFlags) const { THROW_CArg_ExcludedValue; }
CNcbiIostream& CArg_ExcludedValue::AsIOFile(CArgValue::TFileFlags) const { THROW_CArg_ExcludedValue; }
void          CArg_ExcludedValue::CloseFile   (void) const { THROW_CArg_ExcludedValue; }


///////////////////////////////////////////////////////
//  CArg_String::

CArg_String::CArg_String(const string& name, const string& value)
    : CArgValue(name)
{
    m_StringList.push_back(value);
}


bool CArg_String::HasValue(void) const
{
    return !m_StringList.empty();
}


const string& CArg_String::AsString(void) const
{
    if (m_StringList.empty()) {
        return kEmptyStr;
    }
    return m_StringList[0];
}


const CArgValue::TStringArray& CArg_String::GetStringList() const
{
    return m_StringList;
}


CArgValue::TStringArray& CArg_String::SetStringList()
{
    return m_StringList;
}


Int8 CArg_String::AsInt8(void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (Int8) type");
    HBN_ERR("%s", err.c_str());
}

int CArg_String::AsInteger(void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (Integer) type");
    HBN_ERR("%s", err.c_str());
}

TIntId CArg_String::AsIntId(void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (TIntId) type");
    HBN_ERR("%s", err.c_str());
}

double CArg_String::AsDouble(void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (Double) type");
    HBN_ERR("%s", err.c_str());
}

bool CArg_String::AsBoolean(void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (Boolean) type");
    HBN_ERR("%s", err.c_str());
}

const CDir& CArg_String::AsDirectory (void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (CDir) type");
    HBN_ERR("%s", err.c_str());
}

const CTime& CArg_String::AsDateTime (void) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (CTime) type");
    HBN_ERR("%s", err.c_str());
}

CNcbiIstream& CArg_String::AsInputFile(CArgValue::TFileFlags) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (InputFile) type");
    HBN_ERR("%s", err.c_str());
}

CNcbiOstream& CArg_String::AsOutputFile(CArgValue::TFileFlags) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (OutputFile) type");
    HBN_ERR("%s", err.c_str());
}

CNcbiIostream& CArg_String::AsIOFile(CArgValue::TFileFlags) const
{
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (IOFile) type");
    HBN_ERR("%s", err.c_str());
}

void CArg_String::CloseFile(void) const
{
    string sval = AsString();
    string err = string("Attempt to close an argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') of non-file type");
    HBN_ERR("%s", err.c_str());
}



///////////////////////////////////////////////////////
//  CArg_Int8::

CArg_Int8::CArg_Int8(const string& name, const string& value)
    : CArg_String(name, value)
{
    m_Integer = s_StringToInt8(value);
}


Int8 CArg_Int8::AsInt8(void) const
{
    return m_Integer;
}


TIntId CArg_Int8::AsIntId(void) const
{
#ifdef NCBI_INT8_GI
    return m_Integer;
#else
    string sval = AsString();
    string err = string("Attempt to cast argument (name = '")
                 +
                 GetName();
    if (!sval.empty()) err += string("', value = '") + sval;
    err += string("') to a wrong (TIntId) type");
    HBN_ERR("%s", err.c_str());
#endif
}


///////////////////////////////////////////////////////
//  CArg_Integer::

CArg_Integer::CArg_Integer(const string& name, const string& value)
    : CArg_Int8(name, value)
{
    if (m_Integer < kMin_Int  ||  m_Integer > kMax_Int) {
        HBN_ERR("Integer value '%s' of argument '%s' is out of range", value.c_str(), GetName().c_str());
    }
}


int CArg_Integer::AsInteger(void) const
{
    return static_cast<int> (m_Integer);
}


TIntId CArg_Integer::AsIntId(void) const
{
    return AsInteger();
}


///////////////////////////////////////////////////////
//  CArg_IntId::

CArg_IntId::CArg_IntId(const string& name, const string& value)
    : CArg_Int8(name, value)
{
#ifndef NCBI_INT8_GI
    if (m_Integer < kMin_Int || m_Integer > kMax_Int) {
        HBN_ERR("Integer value '%s' of argument '%s' is out of range", value.c_str(), GetName().c_str());
    }
#endif
}


int CArg_IntId::AsInteger(void) const
{
#ifndef NCBI_INT8_GI
    return static_cast<int> (m_Integer);
#else
    NCBI_THROW(CArgException, eWrongCast, s_ArgExptMsg(GetName(),
        "Attempt to cast to a wrong (Integer) type", AsString()));
#endif
}


TIntId CArg_IntId::AsIntId(void) const
{
    return static_cast<TIntId> (m_Integer);
}


///////////////////////////////////////////////////////
//  CArg_DataSize::

CArg_DataSize::CArg_DataSize(const string& name, const string& value)
    : CArg_String(name, value)
{
    m_Integer = NStr::StringToUInt8_DataSize(value);
}

Int8 CArg_DataSize::AsInt8(void) const
{
    return m_Integer;
}


///////////////////////////////////////////////////////
//  CArg_Double::

CArg_Double::CArg_Double(const string& name, const string& value)
    : CArg_String(name, value)
{
    m_Double = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
}


double CArg_Double::AsDouble(void) const
{
    return m_Double;
}



///////////////////////////////////////////////////////
//  CArg_Boolean::

CArg_Boolean::CArg_Boolean(const string& name, bool value)
    : CArg_String(name, NStr::BoolToString(value))
{
    m_Boolean = value;
}


CArg_Boolean::CArg_Boolean(const string& name, const string& value)
    : CArg_String(name, value)
{
    m_Boolean = NStr::StringToBool(value);
}

bool CArg_Boolean::AsBoolean(void) const
{
    return m_Boolean;
}


///////////////////////////////////////////////////////
//  CArg_Flag

CArg_Flag::CArg_Flag(const string& name, bool value)
    : CArg_Boolean(name, value)
{
}
bool CArg_Flag::HasValue(void) const
{
    return AsBoolean();
}

///////////////////////////////////////////////////////
//  CArg_Grid

string s_GridInfoToString(int node_id, int num_nodes)
{
    string gs = NStr::IntToString(node_id) + " " + NStr::IntToString(num_nodes);
    return gs;
}

CArg_Grid::CArg_Grid(const string& name, const string& node_id, const string& num_nodes)
    : CArg_String(name, node_id + " " + num_nodes)
{
    m_NodeId = NStr::StringToInt(node_id);
    m_NumNodes = NStr::StringToInt(num_nodes);
    Validate(m_NodeId, m_NumNodes);
}

CArg_Grid::CArg_Grid(const string& name, const int node_id, const int num_nodes)
    : CArg_String(name, s_GridInfoToString(node_id, num_nodes))
{
    m_NodeId = node_id;
    m_NumNodes = num_nodes;
    Validate(m_NodeId, m_NumNodes);
}

bool CArg_Grid::HasValue() const 
{
    return CArg_String::HasValue();
}

void CArg_Grid::Validate(const int node_id, const int num_nodes)
{
    if (node_id < 0) HBN_ERR("node_id must be >= 0: %s", node_id);
    if (num_nodes <= 0) HBN_ERR("number of nodes must be > 0: %d", num_nodes);
    if (node_id >= num_nodes) HBN_ERR("node id (%d) must be less than number of nodes (%d)",
        node_id, num_nodes);
}

END_NCBI_SCOPE