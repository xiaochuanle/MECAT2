#include "ncbiargs_desc.hpp"

#include "cmdline_flags.hpp"

BEGIN_NCBI_SCOPE
USING_SCOPE(blast);

static const char* s_ExtraName    = "....";

/////////////////////////////////////////////////////////////////////////////
//  CArgs::
//


CArgs::CArgs(void)
{
    m_nExtra = 0;
}


CArgs::~CArgs(void)
{
    return;
}

CArgs& CArgs::Assign(const CArgs& other)
{
    if (this != &other) {
        m_Args = other.m_Args;
        m_nExtra = other.m_nExtra;
        m_Command = other.m_Command;
    }
    return *this;
}

static string s_ComposeNameExtra(size_t idx)
{
    return '#' + NStr::UInt8ToString(idx);
}


inline bool s_IsArgNameChar(char c)
{
    return isalnum(c)  ||  c == '_'  ||  c == '-';
}


CArgs::TArgsCI CArgs::x_Find(const string& name) const
{
    CArgs::TArgsCI arg = m_Args.find(CRef<CArgValue> (new CArg_NoValue(name)));
    if (arg != m_Args.end() || name.empty() || name[0] == '-'  ||
        !s_IsArgNameChar(name[0])) {
        return arg;
    }
    return m_Args.find(CRef<CArgValue> (new CArg_NoValue("-" + name)));
}

CArgs::TArgsI CArgs::x_Find(const string& name)
{
    CArgs::TArgsI arg = m_Args.find(CRef<CArgValue> (new CArg_NoValue(name)));
    if (arg != m_Args.end() || name.empty() || name[0] == '-'  ||
        !s_IsArgNameChar(name[0])) {
        return arg;
    }
    return m_Args.find(CRef<CArgValue> (new CArg_NoValue("-" + name)));
}


bool CArgs::Exist(const string& name) const
{
    return (x_Find(name) != m_Args.end());
}


const CArgValue& CArgs::operator[] (const string& name) const
{
    TArgsCI arg = x_Find(name);
    if (arg == m_Args.end()) {
        // Special diagnostics for "extra" args
        if (!name.empty()  &&  name[0] == '#') {
            size_t idx;
            try {
                idx = NStr::StringToUInt(name.c_str() + 1);
            } catch (...) {
                idx = kMax_UInt;
            }
            if (idx == kMax_UInt) {
                //NCBI_THROW(CArgException, eInvalidArg,
                string err_msg = (
                           "Asked for an argument with invalid name: \"" +
                           name + "\"");
                HBN_ERR("%s", err_msg.c_str());
            }
            if (m_nExtra == 0) {
                //NCBI_THROW(CArgException, eInvalidArg,
                string err_msg = (
                           "No \"extra\" (unnamed positional) arguments "
                           "provided, cannot Get: " + s_ComposeNameExtra(idx));
                HBN_ERR("%s", err_msg.c_str());
            }
            if (idx == 0  ||  idx >= m_nExtra) {
                //NCBI_THROW(CArgException, eInvalidArg,
                string err_msg = (
                           "\"Extra\" (unnamed positional) arg is "
                           "out-of-range (#1.." + s_ComposeNameExtra(m_nExtra)
                           + "): " + s_ComposeNameExtra(idx));
                HBN_ERR("%s", err_msg.c_str());
            }
        }

        // Diagnostics for all other argument classes
        //NCBI_THROW(CArgException, eInvalidArg,
        string err_msg = (
                   "Unknown argument requested: \"" + name + "\"");
        HBN_ERR("%s", err_msg.c_str());
    }

    // Found arg with name "name"
    return **arg;
}


const CArgValue& CArgs::operator[] (size_t idx) const
{
    return (*this)[s_ComposeNameExtra(idx)];
}

vector< CRef<CArgValue> > CArgs::GetAll(void) const
{
    vector< CRef<CArgValue> > res;
    ITERATE( TArgs, a, m_Args) {
        if ((**a).HasValue()) {
            res.push_back( *a );
        }
    }
    return res;
}


string& CArgs::Print(string& str) const
{
    for (TArgsCI arg = m_Args.begin();  arg != m_Args.end();  ++arg) {
        // Arg. name
        const string& arg_name = (*arg)->GetName();
        str += arg_name;

        // Arg. value, if any
        const CArgValue& arg_value = (*this)[arg_name];
        if ( arg_value ) {
            str += " = `";
            string tmp;
            try {
                tmp = NStr::Join( arg_value.GetStringList(), " "); 
            } catch (...) {
                tmp = arg_value.AsString();
            }
            str += tmp;
            str += "'\n";
        } else {
            str += ":  <not assigned>\n";
        }
    }
    return str;
}


void CArgs::Remove(const string& name)
{
    CArgs::TArgsI it =  m_Args.find(CRef<CArgValue> (new CArg_NoValue(name)));
    m_Args.erase(it);
}


void CArgs::Reset(void)
{
    m_nExtra = 0;
    m_Args.clear();
}


void CArgs::Add(CArgValue* arg, bool update, bool add_value)
{
    // special case:  add an "extra" arg (generate virtual name for it)
    bool is_extra = false;
    if ( arg->GetName().empty() ) {
        arg->m_Name = s_ComposeNameExtra(m_nExtra + 1);
        is_extra = true;
    }

    // check-up
    VerifyArgumentName(arg->GetName());
    CArgs::TArgsI arg_it = x_Find(arg->GetName());
    if ( arg_it !=  m_Args.end()) {
        if (update) {
            Remove(arg->GetName());
        } else {
            if (add_value) {
                const string& v = arg->AsString();
                CRef<CArgValue> av = *arg_it;
                av->SetStringList().push_back(v);
            } else {
                //NCBI_THROW(CArgException,eSynopsis,
                string err_msg = (
                   "Argument with this name is defined already: " 
                   + arg->GetName());
                HBN_ERR("%s", err_msg.c_str());
            }
        }
    }

    // add
    arg->SetOrdinalPosition(m_Args.size()+1);
    m_Args.insert(CRef<CArgValue>(arg));

    if ( is_extra ) {
        m_nExtra++;
    }
}


bool CArgs::IsEmpty(void) const
{
    return m_Args.empty();
}

/////////////////////////////////////////////////////////////////////////////
//  Aux.functions to figure out various arg. features
//
//    s_IsPositional(arg)
//    s_IsOptional(arg)
//    s_IsFlag(arg)
//

class CArgDesc;

inline bool s_IsKey(const CArgDesc& arg)
{
    return (dynamic_cast<const CArgDescSynopsis*> (&arg) != 0);
}


inline bool s_IsPositional(const CArgDesc& arg)
{
    return dynamic_cast<const CArgDesc_Pos*> (&arg) &&  !s_IsKey(arg);
}


inline bool s_IsOpening(const CArgDesc& arg)
{
    return dynamic_cast<const CArgDesc_Opening*> (&arg) != NULL;
}


inline bool s_IsOptional(const CArgDesc& arg)
{
    return (dynamic_cast<const CArgDescOptional*> (&arg) != 0);
}


inline bool s_IsFlag(const CArgDesc& arg)
{
    return (dynamic_cast<const CArgDesc_Flag*> (&arg) != 0);
}


inline bool s_IsAlias(const CArgDesc& arg)
{
    return (dynamic_cast<const CArgDesc_Alias*> (&arg) != 0);
}

bool ArgDescIsKey(const CArgDesc& arg)
{
    return s_IsKey(arg);
}

bool ArgDescIsPositional(const CArgDesc& arg)
{
    return s_IsPositional(arg);
}

bool ArgDescIsOpening(const CArgDesc& arg)
{
    return s_IsOpening(arg);
}

bool ArgDescIsOptional(const CArgDesc& arg)
{
    return s_IsOptional(arg);
}

bool ArgDescIsFlag(const CArgDesc& arg)
{
    return s_IsFlag(arg);
}

bool ArgDescIsAlias(const CArgDesc& arg)
{
    return s_IsAlias(arg);
}

/////////////////////////////////////////////////////////////////////////////
//  CArgDesc***::   abstract base classes for argument descriptions
//
//    CArgDesc
//
//    CArgDescMandatory  : CArgDesc
//    CArgDescOptional   : virtual CArgDescMandatory
//    CArgDescDefault    : virtual CArgDescOptional
//
//    CArgDescSynopsis
//


///////////////////////////////////////////////////////
//  CArgDesc::

CArgDesc::CArgDesc(const string& name, const string& comment)
    : m_Name(name), m_Comment(comment)
{
    if (!VerifyArgumentName(m_Name)) {
        HBN_ERR("Invalid argument name: ", m_Name.c_str());
    }
}


CArgDesc::~CArgDesc(void)
{
    return;
}


void CArgDesc::VerifyDefault(void) const
{
    return;
}


void CArgDesc::SetConstraint(const CArgAllow* constraint,
                             CArgDescriptions::EConstraintNegate)
{
    CConstRef<CArgAllow> safe_delete(const_cast<CArgAllow*>(constraint));
    string err_msg = string("No-value argument '")
                     +
                     GetName()
                     +
                     string("' may not be constained");
    if (constraint) err_msg += constraint->GetUsage();
    HBN_ERR("%s", err_msg.c_str());
}


const CArgAllow* CArgDesc::GetConstraint(void) const
{
    return 0;
}


string CArgDesc::GetUsageConstraint(void) const
{
    if (GetFlags() & CArgDescriptions::fConfidential) {
        return kEmptyStr;
    }
    const CArgAllow* constraint = GetConstraint();
    if (!constraint)
        return kEmptyStr;
    string usage;
    if (IsConstraintInverted()) {
        usage = " NOT ";
    }
    usage += constraint->GetUsage();
    return usage;
}


//  Overload the comparison operator -- to handle "AutoPtr<CArgDesc>" elements
//  in "CArgs::m_Args" stored as "set< AutoPtr<CArgDesc> >"
//
inline bool operator< (const AutoPtr<CArgDesc>& x, const AutoPtr<CArgDesc>& y)
{
    return x->GetName() < y->GetName();
}

///////////////////////////////////////////////////////
//  CArgDescMandatory::

CArgDescMandatory::CArgDescMandatory(const string&            name,
                                     const string&            comment,
                                     CArgDescriptions::EType  type,
                                     CArgDescriptions::TFlags flags)
    : CArgDesc(name, comment),
      m_Type(type), m_Flags(flags),
      m_NegateConstraint(CArgDescriptions::eConstraint)
{
    string err_msg;
    // verify if "flags" "type" are matching
    switch ( type ) {
    case CArgDescriptions::eBoolean:
    case CArgDescriptions::eOutputFile:
    case CArgDescriptions::eIOFile:
        return;
    case CArgDescriptions::eInputFile:
        if((flags &
            (CArgDescriptions::fAllowMultiple | CArgDescriptions::fAppend | CArgDescriptions::fTruncate)) == 0)
            return;
        break;
    case CArgDescriptions::k_EType_Size:
        _TROUBLE;
        err_msg = string("Invalid type 'k_EType_Size' for argument '")
                         +
                         GetName()
                         +
                         string("'");
        HBN_ERR("%s", err_msg.c_str());
        /*NOTREACHED*/
        break;
    case CArgDescriptions::eDirectory:
        if ( (flags & ~CArgDescriptions::fCreatePath) == 0 )
            return;
        break;
    default:
        if ( (flags & CArgDescriptions::fFileFlags) == 0 )
            return;
    }

    err_msg = string("Argument type/flags mismatch (type = '")
                     +
                     CArgDescriptions::GetTypeName(type)
                     +
                     string("', flags = '")
                     +
                     NStr::UIntToString(flags)
                     +
                     string("')");
    HBN_ERR("%s", err_msg.c_str());
}


CArgDescMandatory::~CArgDescMandatory(void)
{
    return;
}


string CArgDescMandatory::GetUsageCommentAttr(void) const
{
    CArgDescriptions::EType type = GetType();
    // Print type name
    string str = CArgDescriptions::GetTypeName(type);

    if (type == CArgDescriptions::eDateTime) {
        str += ", format: \"Y-M-DTh:m:gZ\" or \"Y/M/D h:m:gZ\"";
    }
    // Print constraint info, if any
    string constr = GetUsageConstraint();
    if ( !constr.empty() ) {
        str += ", ";
        str += constr;
    }
    return str;
}


CArgValue* CArgDescMandatory::ProcessArgument(const string& value) const
{
    // Process according to the argument type
    unique_ptr<CArgValue> arg_value;
    switch ( GetType() ) {
    case CArgDescriptions::eString:
        arg_value.reset(new CArg_String(GetName(), value));
        break;
    case CArgDescriptions::eBoolean:
        arg_value.reset(new CArg_Boolean(GetName(), value));
        break;
    case CArgDescriptions::eInt8:
        arg_value.reset(new CArg_Int8(GetName(), value));
        break;
    case CArgDescriptions::eInteger:
        arg_value.reset(new CArg_Integer(GetName(), value));
        break;
    case CArgDescriptions::eIntId:
        arg_value.reset(new CArg_IntId(GetName(), value));
        break;
    case CArgDescriptions::eDouble:
        arg_value.reset(new CArg_Double(GetName(), value));
        break;
    case CArgDescriptions::eInputFile: {
        //arg_value = new CArg_InputFile(GetName(), value, GetFlags());
        break;
    }
    case CArgDescriptions::eOutputFile: {
        //arg_value = new CArg_OutputFile(GetName(), value, GetFlags());
        break;
    }
    case CArgDescriptions::eIOFile: {
        //arg_value = new CArg_IOFile(GetName(), value, GetFlags());
        break;
    }
    case CArgDescriptions::eDirectory: {
        //arg_value = new CArg_Dir(GetName(), value, GetFlags());
        break;
    }
    case CArgDescriptions::eDataSize:
        arg_value.reset(new CArg_DataSize(GetName(), value));
        break;
    case CArgDescriptions::eDateTime:
        //arg_value = new CArg_DateTime(GetName(), value);
        break;
    case CArgDescriptions::k_EType_Size: {
        _TROUBLE;
        //NCBI_THROW(CArgException, eArgType, s_ArgExptMsg(GetName(),
        //    "Unknown argument type", NStr::IntToString((int)GetType())));
        string err_msg = string("Unknown type '")
                         +
                         NStr::IntToString(GetType())
                         +
                         string("' for argument '")
                         +
                         GetName()
                         +
                         string("'");
        HBN_ERR("%s", err_msg.c_str());
    }
    } /* switch GetType() */


    // Check against additional (user-defined) constraints, if any imposed
    if ( m_Constraint ) {
        bool err = false;
        try {
            bool check = m_Constraint->Verify(value);
            if (m_NegateConstraint == CArgDescriptions::eConstraintInvert) {
                err = check;
            } else {
                err = !check;
            }

        } catch (...) {
            err = true;
        }

        if (err) {
            if (GetFlags() & CArgDescriptions::fConfidential) {
                //NCBI_THROW(CArgException, eConstraint, s_ArgExptMsg(GetName(),
                //           "Disallowed value",value));
                string err_msg = string("Disallowed value '")
                                 +
                                 value
                                 +
                                 string("' for argument '")
                                 +
                                 GetName()
                                 +
                                 string("'");
                HBN_ERR("%s", err_msg.c_str());
            } else {
                string err_msg;
                if (m_NegateConstraint == CArgDescriptions::eConstraintInvert) {
                    err_msg = GetName() + ": Illegal value, unexpected ";
                } else {
                    err_msg = GetName() + ": Illegal value, expected ";
                }
                err_msg += m_Constraint->GetUsage();
                HBN_ERR("%s", err_msg.c_str());
                //NCBI_THROW(CArgException, eConstraint, s_ArgExptMsg(GetName(),
                //           err_msg + m_Constraint->GetUsage(),value));
            }
        }
    }

    const CArgDescDefault* dflt = dynamic_cast<const CArgDescDefault*> (this);
    if (dflt) {
        arg_value->x_SetDefault(dflt->GetDefaultValue(), false);
    }
    return arg_value.release();
}


CArgValue* CArgDescMandatory::ProcessDefault(void) const
{
    //NCBI_THROW(CArgException, eNoArg,
    //           s_ArgExptMsg(GetName(), "Mandatory value is missing",
    //                        GetUsageCommentAttr()));
    string err_msg = "Mandatory value for argument '"
                     +
                     GetName()
                     +
                     "' is missing";
    HBN_ERR("%s", err_msg.c_str());
}


void CArgDescMandatory::SetConstraint
(const CArgAllow*                    constraint,
 CArgDescriptions::EConstraintNegate negate)
{
    //m_Constraint       = constraint;
    m_Constraint.reset(const_cast<CArgAllow*>(constraint));
    m_NegateConstraint = negate;
}


const CArgAllow* CArgDescMandatory::GetConstraint(void) const
{
    return m_Constraint.get();
}


bool CArgDescMandatory::IsConstraintInverted() const
{
    return m_NegateConstraint == CArgDescriptions::eConstraintInvert;
}


///////////////////////////////////////////////////////
//  CArgDescOptional::


CArgDescOptional::CArgDescOptional(const string&            name,
                                   const string&            comment,
                                   CArgDescriptions::EType  type,
                                   CArgDescriptions::TFlags flags)
    : CArgDescMandatory(name, comment, type, flags),
      m_Group(0)
{
    return;
}


CArgDescOptional::~CArgDescOptional(void)
{
    return;
}


CArgValue* CArgDescOptional::ProcessDefault(void) const
{
    return new CArg_NoValue(GetName());
}




///////////////////////////////////////////////////////
//  CArgDescDefault::


CArgDescDefault::CArgDescDefault(const string&            name,
                                 const string&            comment,
                                 CArgDescriptions::EType  type,
                                 CArgDescriptions::TFlags flags,
                                 const string&            default_value,
                                 const string&            env_var,
                                 const char*              display_value)
    : CArgDescMandatory(name, comment, type, flags),
      CArgDescOptional(name, comment, type, flags),
      m_DefaultValue(default_value), m_EnvVar(env_var),
      m_use_display(display_value != nullptr)
{
    if (m_use_display) {
        m_DisplayValue = display_value;
    }
    return;
}


CArgDescDefault::~CArgDescDefault(void)
{
    return;
}

const string& CArgDescDefault::GetDefaultValue(void) const
{
#ifndef HBN_REMOVE_THIS
    if (!m_EnvVar.empty() && CNcbiApplication::Instance()) {
        const string& value =
            CNcbiApplication::Instance()->GetEnvironment().Get(m_EnvVar);
        if (!value.empty()) {
            return value;
        }
    }
#endif
    return m_DefaultValue;
}

const string& CArgDescDefault::GetDisplayValue(void) const
{
    return m_use_display ?  m_DisplayValue : GetDefaultValue();
}

CArgValue* CArgDescDefault::ProcessDefault(void) const
{
    CArgValue* v = ProcessArgument(GetDefaultValue());
    if (v) {
        v->x_SetDefault(GetDefaultValue(), true);
    }
    return v;
}


void CArgDescDefault::VerifyDefault(void) const
{
    if (GetType() == CArgDescriptions::eInputFile  ||
        GetType() == CArgDescriptions::eOutputFile ||
        GetType() == CArgDescriptions::eIOFile ||
        GetType() == CArgDescriptions::eDirectory) {
        return;
    }

    // Process, then immediately delete
    CRef<CArgValue> arg_value(ProcessArgument(GetDefaultValue()));
}


///////////////////////////////////////////////////////
//  CArgDescSynopsis::


CArgDescSynopsis::CArgDescSynopsis(const string& synopsis)
    : m_Synopsis(synopsis)
{
    for (string::const_iterator it = m_Synopsis.begin();
         it != m_Synopsis.end();  ++it) {
        if (*it != '_'  &&  !isalnum((unsigned char)(*it))) {
            //NCBI_THROW(CArgException,eSynopsis,
            string err_msg = (
                "Argument synopsis must be alphanumeric: "+ m_Synopsis);
            HBN_ERR("%s", err_msg.c_str());
        }
    }
}



/////////////////////////////////////////////////////////////////////////////
//  CArgDesc_***::   classes for argument descriptions
//    CArgDesc_Flag    : CArgDesc
//
//    CArgDesc_Key     : virtual CArgDescMandatory
//    CArgDesc_KeyOpt  : CArgDesc_Key, virtual CArgDescOptional
//    CArgDesc_KeyDef  : CArgDesc_Key, CArgDescDefault
//
//    CArgDesc_Pos     : virtual CArgDescMandatory
//    CArgDesc_PosOpt  : CArgDesc_Pos, virtual CArgDescOptional
//    CArgDesc_PosDef  : CArgDesc_Pos, CArgDescDefault
//


///////////////////////////////////////////////////////
//  CArgDesc_Flag::


CArgDesc_Flag::CArgDesc_Flag(const string& name,
                             const string& comment,
                             bool  set_value)

    : CArgDesc(name, comment),
      m_Group(0),
      m_SetValue(set_value)
{
    return;
}


CArgDesc_Flag::~CArgDesc_Flag(void)
{
    return;
}


string CArgDesc_Flag::GetUsageSynopsis(bool /*name_only*/) const
{
    return "-" + GetName();
}


string CArgDesc_Flag::GetUsageCommentAttr(void) const
{
    return kEmptyStr;
}


CArgValue* CArgDesc_Flag::ProcessArgument(const string& /*value*/) const
{
    CArgValue* v = new CArg_Flag(GetName(), m_SetValue);
    if (v) {
        v->x_SetDefault(NStr::BoolToString(!m_SetValue), false);
    }
    return v;
}


CArgValue* CArgDesc_Flag::ProcessDefault(void) const
{
    CArgValue* v = new CArg_Flag(GetName(), !m_SetValue);
    if (v) {
        v->x_SetDefault(NStr::BoolToString(!m_SetValue), true);
    }
    return v;
}



///////////////////////////////////////////////////////
//  CArgDesc_Pos::


CArgDesc_Pos::CArgDesc_Pos(const string&            name,
                           const string&            comment,
                           CArgDescriptions::EType  type,
                           CArgDescriptions::TFlags flags)
    : CArgDescMandatory(name, comment, type, flags)
{
    return;
}


CArgDesc_Pos::~CArgDesc_Pos(void)
{
    return;
}


string CArgDesc_Pos::GetUsageSynopsis(bool /*name_only*/) const
{
    return GetName().empty() ? s_ExtraName : GetName();
}


///////////////////////////////////////////////////////
//  CArgDesc_Opening::


CArgDesc_Opening::CArgDesc_Opening(const string&            name,
                           const string&            comment,
                           CArgDescriptions::EType  type,
                           CArgDescriptions::TFlags flags)
    : CArgDescMandatory(name, comment, type, flags)
{
    return;
}


CArgDesc_Opening::~CArgDesc_Opening(void)
{
    return;
}


string CArgDesc_Opening::GetUsageSynopsis(bool /*name_only*/) const
{
    return GetName().empty() ? s_ExtraName : GetName();
}



///////////////////////////////////////////////////////
//  CArgDesc_PosOpt::


CArgDesc_PosOpt::CArgDesc_PosOpt(const string&            name,
                                 const string&            comment,
                                 CArgDescriptions::EType  type,
                                 CArgDescriptions::TFlags flags)
    : CArgDescMandatory (name, comment, type, flags),
      CArgDescOptional  (name, comment, type, flags),
      CArgDesc_Pos      (name, comment, type, flags)
{
    return;
}


CArgDesc_PosOpt::~CArgDesc_PosOpt(void)
{
    return;
}



///////////////////////////////////////////////////////
//  CArgDesc_PosDef::


CArgDesc_PosDef::CArgDesc_PosDef(const string&            name,
                                 const string&            comment,
                                 CArgDescriptions::EType  type,
                                 CArgDescriptions::TFlags flags,
                                 const string&            default_value,
                                 const string&            env_var,
                                 const char*              display_value)
    : CArgDescMandatory (name, comment, type, flags),
      CArgDescOptional  (name, comment, type, flags),
      CArgDescDefault   (name, comment, type, flags, default_value, env_var, display_value),
      CArgDesc_PosOpt   (name, comment, type, flags)
{
    return;
}


CArgDesc_PosDef::~CArgDesc_PosDef(void)
{
    return;
}



///////////////////////////////////////////////////////
//  CArgDesc_Key::


CArgDesc_Key::CArgDesc_Key(const string&            name,
                           const string&            comment,
                           CArgDescriptions::EType  type,
                           CArgDescriptions::TFlags flags,
                           const string&            synopsis)
    : CArgDescMandatory(name, comment, type, flags),
      CArgDesc_Pos     (name, comment, type, flags),
      CArgDescSynopsis(synopsis)
{
    return;
}


CArgDesc_Key::~CArgDesc_Key(void)
{
    return;
}


inline string s_KeyUsageSynopsis(const string& name, const string& synopsis,
                                 bool name_only,
                                 CArgDescriptions::TFlags flags)
{
    if ( name_only ) {
        return '-' + name;
    } else {
        char separator =
            (flags & CArgDescriptions::fMandatorySeparator) ? '=' : ' ';
        return '-' + name + separator + synopsis;
    }
}


string CArgDesc_Key::GetUsageSynopsis(bool name_only) const
{
    return s_KeyUsageSynopsis(GetName(), GetSynopsis(), name_only, GetFlags());
}



///////////////////////////////////////////////////////
//  CArgDesc_KeyOpt::


CArgDesc_KeyOpt::CArgDesc_KeyOpt(const string&            name,
                                 const string&            comment,
                                 CArgDescriptions::EType  type,
                                 CArgDescriptions::TFlags flags,
                                 const string&            synopsis)
    : CArgDescMandatory(name, comment, type, flags),
      CArgDescOptional (name, comment, type, flags),
      CArgDesc_PosOpt  (name, comment, type, flags),
      CArgDescSynopsis(synopsis)
{
    return;
}


CArgDesc_KeyOpt::~CArgDesc_KeyOpt(void)
{
    return;
}


string CArgDesc_KeyOpt::GetUsageSynopsis(bool name_only) const
{
    return s_KeyUsageSynopsis(GetName(), GetSynopsis(), name_only, GetFlags());
}



///////////////////////////////////////////////////////
//  CArgDesc_KeyDef::


CArgDesc_KeyDef::CArgDesc_KeyDef(const string&            name,
                                 const string&            comment,
                                 CArgDescriptions::EType  type,
                                 CArgDescriptions::TFlags flags,
                                 const string&            synopsis,
                                 const string&            default_value,
                                 const string&            env_var,
                                 const char*              display_value)
    : CArgDescMandatory(name, comment, type, flags),
      CArgDescOptional (name, comment, type, flags),
      CArgDesc_PosDef  (name, comment, type, flags, default_value, env_var, display_value),
      CArgDescSynopsis(synopsis)
{
    return;
}


CArgDesc_KeyDef::~CArgDesc_KeyDef(void)
{
    return;
}


string CArgDesc_KeyDef::GetUsageSynopsis(bool name_only) const
{
    return s_KeyUsageSynopsis(GetName(), GetSynopsis(), name_only, GetFlags());
}


///////////////////////////////////////////////////////
//  CArgDesc_Alias::

CArgDesc_Alias::CArgDesc_Alias(const string& alias,
                               const string& arg_name,
                               const string& comment)
    : CArgDesc(alias, comment),
      m_ArgName(arg_name),
      m_NegativeFlag(false)
{
}


CArgDesc_Alias::~CArgDesc_Alias(void)
{
}


const string& CArgDesc_Alias::GetAliasedName(void) const
{
    return m_ArgName;
}


string CArgDesc_Alias::GetUsageSynopsis(bool /*name_only*/) const
{
    return kEmptyStr;
}


string CArgDesc_Alias::GetUsageCommentAttr(void) const
{
    return kEmptyStr;
}

    
CArgValue* CArgDesc_Alias::ProcessArgument(const string& /*value*/) const
{
    return new CArg_NoValue(GetName());
}


CArgValue* CArgDesc_Alias::ProcessDefault(void) const
{
    return new CArg_NoValue(GetName());
}


///////////////////////////////////////////////////////
//  CArgDescriptions::
//


CArgDescriptions::CArgDescriptions(bool              auto_help)
    : m_ArgsType(eRegularArgs),
      m_nExtra(0),
      m_nExtraOpt(0),
      m_CurrentGroup(0),
      m_PositionalMode(ePositionalMode_Strict),
      m_MiscFlags(fMisc_Default),
      m_AutoHelp(auto_help)
{
    SetUsageContext("NCBI_PROGRAM", kEmptyStr);
    m_ArgGroups.push_back(kEmptyStr); // create default group #0
    AddFlag(kArgHelp,
            "Print USAGE and DESCRIPTION;  ignore all other parameters");
    AddFlag(kArgFullHelp,
            "Print USAGE, DESCRIPTION and ARGUMENTS;"
            " ignore all other parameters");
    AddFlag(kArgVersion, "Print version number; ignore other arguments");
}


CArgDescriptions::~CArgDescriptions(void)
{
    return;
}


void CArgDescriptions::SetArgsType(EArgSetType args_type)
{
    m_ArgsType = args_type;

    // Run args check for a CGI application
    if (m_ArgsType == eCgiArgs) {
        // Must have no named positional arguments
        if ( !m_PosArgs.empty() ) {
            string err_msg = 
                       "CGI application cannot have positional arguments, "
                       "name of the offending argument: '"
                       + *m_PosArgs.begin() + "'.";
            HBN_ERR("%s", err_msg.c_str());
        }

        // Must have no flag arguments
        ITERATE (TArgs, it, m_Args) {
            const CArgDesc& arg = **it;
            if ( s_IsFlag(arg) ) {
                const string& name = arg.GetName();

                if (name == kArgHelp || name == kArgFullHelp)  // help
                    continue;

                if (name == kArgVersion)  // version
                    continue; 

                string err_msg = 
                           "CGI application cannot have flag arguments, "
                           "name of the offending flag: '" + name + "'.";
                HBN_ERR("%s", err_msg.c_str());
            }
        }

        // Must have no unnamed positional arguments
        if (m_nExtra  ||  m_nExtraOpt) {
            string err_msg = 
                       "CGI application cannot have unnamed positional "
                       "arguments.";
            HBN_ERR("%s", err_msg.c_str());
        }
    }
}


const char* CArgDescriptions::GetTypeName(EType type)
{
    static const char* s_TypeName[k_EType_Size] = {
        "String",
        "Boolean",
        "Int8",
        "Integer",
        "IntId",
        "Real",
        "File_In",
        "File_Out",
        "File_IO",
        "Directory",
        "DataSize",
        "DateTime"
    };

    if (type == k_EType_Size) {
        _TROUBLE;
        string err_msg = 
                   "Invalid argument type: k_EType_Size";
        HBN_ERR("%s", err_msg.c_str());
    }
    return s_TypeName[(int) type];
}


void CArgDescriptions::AddKey
(const string& name,
 const string& synopsis,
 const string& comment,
 EType         type,
 TFlags        flags)
{
    unique_ptr<CArgDesc_Key> arg(new CArgDesc_Key(name,
        comment, type, flags, synopsis));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddOptionalKey
(const string& name,
 const string& synopsis,
 const string& comment,
 EType         type,
 TFlags        flags)
{
    unique_ptr<CArgDesc_KeyOpt> arg(new CArgDesc_KeyOpt(name,
        comment, type, flags, synopsis));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddDefaultKey
(const string& name,
 const string& synopsis,
 const string& comment,
 EType         type,
 const string& default_value,
 TFlags        flags,
 const string& env_var,
 const char*   display_value)
{
    unique_ptr<CArgDesc_KeyDef> arg(new CArgDesc_KeyDef(name,
        comment, type, flags, synopsis, default_value, env_var, display_value));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddFlag(
    const string& name,
    const string& comment,
    CBoolEnum<EFlagValue> set_value)
{
    unique_ptr<CArgDesc_Flag> arg(new CArgDesc_Flag(name, comment, set_value == eFlagHasValueIfSet));
    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddPositional(
    const string& name,
    const string& comment,
    EType         type,
    TFlags        flags)
{
    unique_ptr<CArgDesc_Pos> arg(new CArgDesc_Pos(name, comment, type, flags));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddOpening(
    const string& name,
    const string& comment,
    EType         type,
    TFlags        flags)
{
    unique_ptr<CArgDesc_Opening> arg(new CArgDesc_Opening(name, comment, type, flags));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddOptionalPositional(
    const string& name,
    const string& comment,
    EType         type,
    TFlags        flags)
{
    unique_ptr<CArgDesc_PosOpt> arg
        (new CArgDesc_PosOpt(name, comment, type, flags));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddDefaultPositional(
     const string& name,
     const string& comment,
     EType         type,
     const string& default_value,
     TFlags        flags,
     const string& env_var,
     const char*   display_value)
{
    unique_ptr<CArgDesc_PosDef> arg(new CArgDesc_PosDef(name,
        comment, type, flags, default_value, env_var, display_value));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddExtra(
    unsigned      n_mandatory,
    unsigned      n_optional,
    const string& comment,
    EType         type,
    TFlags        flags)
{
    if (!n_mandatory  &&  !n_optional) {
        string err_msg = (
            "Number of extra arguments cannot be zero");
        HBN_ERR("%s", err_msg.c_str());
    }
    if (n_mandatory > 4096) {
        string err_msg = (
            "Number of mandatory extra arguments is too big");
        HBN_ERR("%s", err_msg.c_str());
    }

    m_nExtra    = n_mandatory;
    m_nExtraOpt = n_optional;

    unique_ptr<CArgDesc_Pos> arg
        (m_nExtra ?
         new CArgDesc_Pos   (kEmptyStr, comment, type, flags) :
         new CArgDesc_PosOpt(kEmptyStr, comment, type, flags));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddAlias(const string& alias,
                                const string& arg_name)
{
    unique_ptr<CArgDesc_Alias>arg
        (new CArgDesc_Alias(alias, arg_name, kEmptyStr));

    x_AddDesc(*arg);
    arg.release();
}


void CArgDescriptions::AddNegatedFlagAlias(const string& alias,
                                           const string& arg_name,
                                           const string& comment)
{
    // Make sure arg_name describes a flag
    TArgsCI orig = x_Find(arg_name);
    if (orig == m_Args.end()  ||  !s_IsFlag(**orig)) {
        string err_msg = (
            "Attempt to negate a non-flag argument: " + arg_name);
        HBN_ERR("%s", err_msg.c_str());
    }

    unique_ptr<CArgDesc_Alias> arg(new CArgDesc_Alias(alias, arg_name, comment));
    arg->SetNegativeFlag(true);

    x_AddDesc(*arg);
    arg.release();
}

void CArgDescriptions::AddDependencyGroup(CArgDependencyGroup* dep_group)
{
    m_DependencyGroups.insert( CConstRef<CArgDependencyGroup>(dep_group));
}

void CArgDescriptions::SetConstraint(const string&      name, 
                                     const CArgAllow*   constraint,
                                     EConstraintNegate  negate)
{
    TArgsI it = x_Find(name);
    if (it == m_Args.end()) {
        CConstRef<CArgAllow> safe_delete(const_cast<CArgAllow*>(constraint));  // delete, if last ref
        string err_msg = (
            "Attempt to set constraint for undescribed argument: " + name);
        HBN_ERR("%s", err_msg.c_str());
    }
    (*it)->SetConstraint(constraint, negate);
}


void CArgDescriptions::SetConstraint(const string&      name,
                                     const CArgAllow&   constraint,
                                     EConstraintNegate  negate)
{
    CArgAllow* onheap = constraint.Clone();
    if ( !onheap ) {
        string err_msg = (
                   "Clone method not implemented for: " + name);
        HBN_ERR("%s", err_msg.c_str());
    }
    SetConstraint(name, onheap, negate);
}


void CArgDescriptions::SetDependency(const string& arg1,
                                     EDependency   dep,
                                     const string& arg2)
{
    m_Dependencies.insert(TDependencies::value_type(arg1,
        SArgDependency(arg2, dep)));
    if (dep == eExcludes) {
        // Exclusions must work in both directions
        m_Dependencies.insert(TDependencies::value_type(arg2,
            SArgDependency(arg1, dep)));
    }
}


void CArgDescriptions::SetCurrentGroup(const string& group)
{
    m_CurrentGroup = x_GetGroupIndex(group);
    if (m_CurrentGroup >= m_ArgGroups.size()) {
        m_ArgGroups.push_back(group);
        m_CurrentGroup = m_ArgGroups.size() - 1;
    }
}


bool CArgDescriptions::Exist(const string& name) const
{
    return (x_Find(name) != m_Args.end());
}


void CArgDescriptions::Delete(const string& name)
{
    {{ // ...from the list of all args
        TArgsI it = x_Find(name);
        if (it == m_Args.end()) {
            string err_msg = (
                "Argument description is not found");
            HBN_ERR("%s", err_msg.c_str());
        }
        m_Args.erase(it);
        if (name == kArgHelp) {
            m_AutoHelp = false;
        }

        // take special care of the extra args
        if ( name.empty() ) {
            m_nExtra = 0;
            m_nExtraOpt = 0;
            return;
        }
    }}

    {{ // ...from the list of key/flag args
        TKeyFlagArgs::iterator it =
            find(m_KeyFlagArgs.begin(), m_KeyFlagArgs.end(), name);
        if (it != m_KeyFlagArgs.end()) {
            m_KeyFlagArgs.erase(it);
            _ASSERT(find(m_KeyFlagArgs.begin(), m_KeyFlagArgs.end(), name) ==
                         m_KeyFlagArgs.end());
            _ASSERT(find(m_PosArgs.begin(), m_PosArgs.end(), name) ==
                         m_PosArgs.end());
            return;
        }
    }}

    {{ // ...from the list of positional args' positions
        TPosArgs::iterator it =
            find(m_PosArgs.begin(), m_PosArgs.end(), name);
        _ASSERT (it != m_PosArgs.end());
        m_PosArgs.erase(it);
        _ASSERT(find(m_PosArgs.begin(), m_PosArgs.end(), name) ==
                     m_PosArgs.end());
    }}
}


// Fake class to hold only a name -- to find in "m_Args"
class CArgDesc_NameOnly : public CArgDesc
{
public:
    CArgDesc_NameOnly(const string& name) :
        CArgDesc(name, kEmptyStr) {}
private:
    virtual string GetUsageSynopsis(bool/*name_only*/) const{return kEmptyStr;}
    virtual string GetUsageCommentAttr(void) const {return kEmptyStr;}
    virtual CArgValue* ProcessArgument(const string&) const {return 0;}
    virtual CArgValue* ProcessDefault(void) const {return 0;}
};

CArgDescriptions::TArgsCI CArgDescriptions::x_Find(const string& name,
                                                   bool* negative) const
{
    CArgDescriptions::TArgsCI arg =
        m_Args.find(AutoPtr<CArgDesc> (new CArgDesc_NameOnly(name)));
    if ( arg != m_Args.end() ) {
        const CArgDesc_Alias* al =
            dynamic_cast<const CArgDesc_Alias*>(arg->get());
        if ( al ) {
            if ( negative ) {
                *negative = al->GetNegativeFlag();
            }
            return x_Find(al->GetAliasedName(), negative);
        }
    }
    return arg;
}

CArgDescriptions::TArgsI CArgDescriptions::x_Find(const string& name,
                                                   bool* negative)
{
    CArgDescriptions::TArgsI arg =
        m_Args.find(AutoPtr<CArgDesc> (new CArgDesc_NameOnly(name)));
    if ( arg != m_Args.end() ) {
        const CArgDesc_Alias* al =
            dynamic_cast<const CArgDesc_Alias*>(arg->get());
        if ( al ) {
            if ( negative ) {
                *negative = al->GetNegativeFlag();
            }
            return x_Find(al->GetAliasedName(), negative);
        }
    }
    return arg;
}


size_t CArgDescriptions::x_GetGroupIndex(const string& group) const
{
    if ( group.empty() ) {
        return 0;
    }
    for (size_t i = 1; i < m_ArgGroups.size(); ++i) {
        if ( NStr::EqualNocase(m_ArgGroups[i], group) ) {
            return i;
        }
    }
    return m_ArgGroups.size();
}


void CArgDescriptions::x_PreCheck(void) const
{
    // Check for the consistency of positional args
    if ( m_nExtra ) {
        for (TPosArgs::const_iterator name = m_PosArgs.begin();
             name != m_PosArgs.end();  ++name) {
            TArgsCI arg_it = x_Find(*name);
            _ASSERT(arg_it != m_Args.end());
            CArgDesc& arg = **arg_it;

            if (dynamic_cast<const CArgDesc_PosOpt*> (&arg)) {
                string err_msg = (
                    "Having both optional named and required unnamed "
                    "positional arguments is prohibited");
                HBN_ERR("%s", err_msg.c_str());
            }
        }
    }

    // Check for the validity of default values.
    // Also check for conflict between no-separator and regular names
    for (TArgsCI it = m_Args.begin();  it != m_Args.end();  ++it) {
        CArgDesc& arg = **it;

        const string& name = arg.GetName();
        if (name.length() > 1  &&  m_NoSeparator.find(name[0]) != NPOS) {
            // find the argument with optional separator and check its flags
            for (TArgsCI i = m_Args.begin();  i != m_Args.end();  ++i) {
                CArgDesc& a = **i;
                const string& n = a.GetName();
                if (n.length() == 1 && n[0] == name[0] &&
                    (a.GetFlags() & CArgDescriptions::fOptionalSeparator)) {
                    if ((a.GetFlags() & CArgDescriptions::fOptionalSeparatorAllowConflict) == 0) {
                        //NCBI_THROW(CArgException, eInvalidArg,
                        string err_msg = (
                            string("'") + name[0] +
                            "' argument allowed to contain no separator conflicts with '" +
                            name + "' argument. To allow such conflicts, add" +
                            " CArgDescriptions::fOptionalSeparatorAllowConflict flag into" +
                            " description of '" + name[0] + "'.");
                        HBN_ERR("%s", err_msg.c_str());
                    }
                    break;
                }
            }
        }

/*
        if (dynamic_cast<CArgDescDefault*> (&arg) == 0) {
            continue;
        }
*/

        arg.VerifyDefault();
    }
}


CArgs* CArgDescriptions::CreateArgs(int argc, char* argv[]) const
{
    //const_cast<CArgDescriptions&>(*this).SetCurrentGroup(kEmptyStr);
    //return CreateArgs(argc, argv);
    return NULL;
}


void CArgDescriptions::x_CheckAutoHelp(const string& arg) const
{
    if (arg.compare(string("-") + kArgHelp) == 0) {
        if (m_AutoHelp) {
            //NCBI_THROW(CArgHelpException,eHelp,kEmptyStr);
        }
    } else if (arg.compare(string("-") + kArgFullHelp) == 0) {
        //NCBI_THROW(CArgHelpException,eHelpFull,kEmptyStr);
    }
}


// (return TRUE if "arg2" was used)
bool CArgDescriptions::x_CreateArg(const string& arg1,
                                   bool have_arg2, const string& arg2,
                                   unsigned* n_plain, CArgs& args) const
{
    // Argument name
    string name;
    bool is_keyflag = false;

    // Check if to start processing the args as positional
    if (*n_plain == kMax_UInt || m_PositionalMode == ePositionalMode_Loose) {
        // Check for the "--" delimiter
        if (arg1.compare("--") == 0) {
            if (*n_plain == kMax_UInt) {
                *n_plain = 0;  // pos.args started
            }
            return false;
        }
        size_t  argssofar = args.GetAll().size();
        // Check if argument has key/flag syntax
        if ((arg1.length() > 1)  &&  arg1[0] == '-') {
            name = arg1.substr(1);
            TArgsCI it = m_Args.end();
            //try {
                it = x_Find(name);
            //} catch (CArgException&) {
            //}
            if (it == m_Args.end()) {
                if (m_OpeningArgs.size() > argssofar) {
                    return x_CreateArg(arg1, m_OpeningArgs[argssofar], have_arg2, arg2, *n_plain, args);
                }
            }
            // Check for '=' in the arg1
            size_t eq = name.find('=');
            if (eq != NPOS) {
                name = name.substr(0, eq);
            }
            if (m_PositionalMode == ePositionalMode_Loose) {
                is_keyflag = x_Find(name) != m_Args.end();
                // If not a valid key/flag, treat it as a positional value
                if (!VerifyArgumentName(name)  ||  !is_keyflag) {
                    if (*n_plain == kMax_UInt) {
                        *n_plain = 0;  // pos.args started
                    }
                }
            }
        } else {
            if (m_OpeningArgs.size() > argssofar) {
                return x_CreateArg(arg1, m_OpeningArgs[argssofar], have_arg2, arg2, *n_plain, args);
            }
            if (*n_plain == kMax_UInt) {
                *n_plain = 0;  // pos.args started
            }
        }
    }

    // Whether the value of "arg2" is used
    bool arg2_used = false;

    // Extract name of positional argument
    if (*n_plain != kMax_UInt && !is_keyflag) {
        if (*n_plain < m_PosArgs.size()) {
            name = m_PosArgs[*n_plain];  // named positional argument
        } else {
            name = kEmptyStr;  // unnamed (extra) positional argument
        }
        (*n_plain)++;

        // Check for too many positional arguments
        if (kMax_UInt - m_nExtraOpt > m_nExtra + m_PosArgs.size()  &&
            *n_plain > m_PosArgs.size() + m_nExtra + m_nExtraOpt) {
            //NCBI_THROW(CArgException,eSynopsis,
            string err_msg = (
                "Too many positional arguments (" +
                NStr::UIntToString(*n_plain) +
                "), the offending value: "+ arg1);
            HBN_ERR("%s", err_msg.c_str());
        }
    }

    arg2_used = x_CreateArg(arg1, name, have_arg2, arg2, *n_plain, args);

    // Success (also indicate whether one or two "raw" args have been used)
    return arg2_used;
}


bool CArgDescriptions::x_CreateArg(const string& arg1,
                                   const string& name_in, 
                                   bool          have_arg2,
                                   const string& arg2,
                                   unsigned      n_plain,
                                   CArgs&        args,
                                   bool          update,
                                   CArgValue**   new_value) const
{
#ifndef HBN_REMOVE_THIS
    if (new_value)
        *new_value = 0;

    string name(name_in);
    bool arg2_used = false;
    bool no_separator = false;
    bool eq_separator = false;
    bool negative = false;

    // Get arg. description
    TArgsCI it;
    try {
        it = x_Find(name, &negative);
    } catch (CArgException&) {
        // Suppress overzealous "invalid argument name" exceptions
        // in the no-separator case.
        if (m_NoSeparator.find(name[0]) != NPOS) {
            it = m_Args.end(); // avoid duplicating the logic below
        } else {
            throw;
        }
    }

    // Check for '/' in the arg1
    bool confidential = it != m_Args.end() &&
        ((*it)->GetFlags() & CArgDescriptions::fConfidential) != 0;
    char conf_method = confidential ? 't' : '\0';
    size_t dash = name.rfind('-');
    if (it == m_Args.end() && dash != NPOS && dash != 0) {
        string test(name.substr(0, dash));
        string suffix(name.substr(dash+1));
        if (NStr::strcasecmp(suffix.c_str(), "file") == 0 ||
            NStr::strcasecmp(suffix.c_str(), "verbatim") == 0)
        {
            try {
                it = x_Find(test);
            } catch (CArgException&) {
                it = m_Args.end();
            }
            if (it != m_Args.end()) {
// verify that it has Confidential flag
// and there is something after dash
                if (((*it)->GetFlags() & CArgDescriptions::fConfidential) &&
                    name.size() > (dash+1)) {
                    confidential = true;
                    conf_method = name[dash+1];
                    name = test;
                }
            }
        }
    }


    if (it == m_Args.end()  &&  m_NoSeparator.find(name[0]) != NPOS) {
        it = x_Find(name.substr(0, 1), &negative);
        _ASSERT(it != m_Args.end());
        no_separator = true;
    }
    if (it == m_Args.end()) {
        if ( name.empty() ) {
            NCBI_THROW(CArgException,eInvalidArg,
                    "Unexpected extra argument, at position # " +
                    NStr::UIntToString(n_plain));
        } else {
            NCBI_THROW(CArgException,eInvalidArg,
                    "Unknown argument: \"" + name + "\"");
        }
    }
    _ASSERT(*it);

    const CArgDesc& arg = **it;

    if ( s_IsFlag(arg) ) {
        x_CheckAutoHelp(arg1);
    }

    // Check value separated by '='
    string arg_val;
    if ( s_IsKey(arg) && !confidential) {
        eq_separator = arg1.length() > name.length()  &&
            (arg1[name.length() + 1] == '=');
        if ( !eq_separator ) {
            if ((arg.GetFlags() & fMandatorySeparator) != 0) {
                NCBI_THROW(CArgException,eInvalidArg,
                    "Invalid argument: " + arg1);
            }
            no_separator |= (arg.GetFlags() & fOptionalSeparator) != 0  &&
                name.length() == 1  &&  arg1.length() > 2;
        }
    }

    // Get argument value
    string value;
    if ( !eq_separator  &&  !no_separator ) {
        if ( !s_IsKey(arg)  || (confidential && conf_method == 't')) {
            value = arg1;
        }
        else {
            // <key> <value> arg  -- advance from the arg.name to the arg.value
            if ( !have_arg2 ) {

                // if update specified we try to add default value
                //  (mandatory throws an exception out of the ProcessDefault())
                if (update) {
                    CRef<CArgValue> arg_value(arg.ProcessDefault());
                    // Add the value to "args"
                    args.Add(arg_value, update);
                    return arg2_used;
                }

                NCBI_THROW(CArgException,eNoArg,s_ArgExptMsg(arg1,
                    "Value is missing", kEmptyStr));
            }
            value = arg2;
            arg2_used = true;
        }
    }
    else {
        _ASSERT(s_IsKey(arg));
        if ( no_separator ) {
            arg_val = arg1.substr(2);
        }
        else {
            arg_val = arg1.substr(name.length() + 2);
        }
        value = arg_val;
    }

    if (confidential) {
        switch (conf_method) {
        default:
            break;
        case 'f':
        case 'F':
            value = (value != "-") ? s_CArgs_ReadFromFile(name, value)
                                   : s_CArgs_ReadFromStdin(name, eNoEcho, "");
            break;
        case 't':
        case 'T':
            value = s_CArgs_ReadFromConsole(name, eNoEcho, nullptr);
            break;
        }
    }

    CArgValue* av = 0;
    try {
        // Process the "raw" argument value into "CArgValue"
        if ( negative  &&  s_IsFlag(arg) ) {
            // Negative flag - use default value rather than
            // normal one.
            av = arg.ProcessDefault();
        }
        else {
            av = arg.ProcessArgument(value);
        }
    }
    catch (CArgException) {
        const CArgErrorHandler* err_handler = arg.GetErrorHandler();
        if ( !err_handler ) {
            err_handler = m_ErrorHandler.GetPointerOrNull();
        }
        _ASSERT(err_handler);
        av = err_handler->HandleError(arg, value);
    }

    if ( !av ) {
        return arg2_used;
    }
    CRef<CArgValue> arg_value(av);

    if (new_value) {
        *new_value = av;
    }

    bool allow_multiple = false;
    const CArgDescMandatory* adm = 
        dynamic_cast<const CArgDescMandatory*>(&arg);

    if (adm) {
        allow_multiple = 
            (adm->GetFlags() & CArgDescriptions::fAllowMultiple) != 0;
    }

    // Add the argument value to "args"
    args.Add(arg_value, update, allow_multiple);

    return arg2_used;
#else
    return false;
#endif
}


bool CArgDescriptions::x_IsMultiArg(const string& name) const
{
    TArgsCI it = x_Find(name);
    if (it == m_Args.end()) {
        return false;
    }
    const CArgDesc& arg = **it;
    const CArgDescMandatory* adm = 
        dynamic_cast<const CArgDescMandatory*>(&arg);

    if (!adm) {
        return false;
    }
    return (adm->GetFlags() & CArgDescriptions::fAllowMultiple) != 0;
}

#ifndef HBN_REMOVE_THIS
void CArgDescriptions::x_PostCheck(CArgs&           args,
                                   unsigned int     n_plain,
                                   EPostCheckCaller caller)
    const
{
    // If explicitly specified, printout usage and exit in case there
    // was no args passed to the application
    if (IsSetMiscFlag(fUsageIfNoArgs)  &&  args.IsEmpty()) {
        NCBI_THROW(CArgHelpException, eHelp, kEmptyStr);
    }

    // Check dependencies, create set of exclusions
    unsigned int nExtra = m_nExtra;
    string nameReq, nameExc;
    unsigned int nExtraReq = 0;
    unsigned int nExtraExc = 0;
    unsigned int nEx = 0;
    set<string> exclude;
    map<string,string> require;
    ITERATE(TDependencies, dep, m_Dependencies) {
        // Skip missing and empty arguments
        if (!args.Exist(dep->first)  ||  !args[dep->first]) {
            continue;
        }
        switch ( dep->second.m_Dep ) {
        case eRequires:
            require.insert(make_pair(dep->second.m_Arg,dep->first));
            if (dep->second.m_Arg.at(0) == '#') {
                nEx = NStr::StringToUInt(
                    CTempString(dep->second.m_Arg.data()+1, dep->second.m_Arg.size()-1));
                if (nEx > nExtraReq) {
                    nExtraReq = nEx;
                    nameReq = dep->first;
                }
            }
            break;
        case eExcludes:
            // Excluded exists and is not empty?
            if (args.Exist(dep->second.m_Arg)  &&  args[dep->second.m_Arg]) {
                NCBI_THROW(CArgException, eConstraint,
                    s_ArgExptMsg(dep->second.m_Arg,
                    "Incompatible with argument", dep->first));
            }
            exclude.insert(dep->second.m_Arg);
            if (dep->second.m_Arg.at(0) == '#') {
                nEx = NStr::StringToUInt(
                    CTempString(dep->second.m_Arg.data()+1, dep->second.m_Arg.size()-1));
                if (nEx > nExtraExc) {
                    nExtraExc = nEx;
                    nameExc = dep->first;
                }
            }
            break;
        }
    }
    if (nExtraReq != 0 && nExtraExc != 0 && nExtraReq >= nExtraExc) {
        NCBI_THROW(CArgException,eSynopsis,
            "Conflicting dependencies on unnamed positional arguments: " + 
            nameReq + " requires #" + NStr::UIntToString(nExtraReq) + ", " +
            nameExc + " excludes #" + NStr::UIntToString(nExtraExc));
    }
    nExtra = max(nExtra, nExtraReq);
    if (nExtraExc > 0) {
        nExtra = max(nExtra, nExtraExc-1);
    }

    // Check that all opening args are provided
    ITERATE (TPosArgs, it, m_OpeningArgs) {
        if (!args.Exist(*it)) {
            NCBI_THROW(CArgException,eNoArg, "Opening argument not provided: " + *it);
        }
    }

    // Check if all mandatory unnamed positional arguments are provided
    // note that positional ones are filled first, no matter are they optional or not
    if (m_PosArgs.size() <= n_plain  &&
        n_plain < m_PosArgs.size() + nExtra){
        NCBI_THROW(CArgException,eNoArg,
            "Too few (" + NStr::NumericToString(n_plain - m_PosArgs.size()) +
            ") unnamed positional arguments. Must define at least " +
            NStr::NumericToString(nExtra));
    }

    // Compose an ordered list of args
    list<const CArgDesc*> def_args;
    ITERATE (TKeyFlagArgs, it, m_KeyFlagArgs) {
        const CArgDesc& arg = **x_Find(*it);
        def_args.push_back(&arg);
    }
    ITERATE (TPosArgs, it, m_PosArgs) {
        const CArgDesc& arg = **x_Find(*it);
        def_args.push_back(&arg);
    }

    for (set< CConstRef<CArgDependencyGroup> >::const_iterator i = m_DependencyGroups.begin();
        i != m_DependencyGroups.end(); ++i) {
        i->GetPointer()->Evaluate(args);
    }

    // Set default values (if available) for the arguments not defined
    // in the command line.
    ITERATE (list<const CArgDesc*>, it, def_args) {
        const CArgDesc& arg = **it;

        // Nothing to do if defined in the command-line
        if ( args.Exist(arg.GetName()) ) {
            continue;
        }

        if (require.find(arg.GetName()) != require.end() ) {
            string requester(require.find(arg.GetName())->second);
            // Required argument must be present
            NCBI_THROW(CArgException, eConstraint,
                s_ArgExptMsg(arg.GetName(),
                "Must be specified, as it is required by argument", requester));
        }

        if (exclude.find(arg.GetName()) != exclude.end()) {
            CRef<CArg_ExcludedValue> arg_exvalue(
                new CArg_ExcludedValue(arg.GetName()));
            // Add the excluded-value argument to "args"
            args.Add(arg_exvalue);
            continue;
        }
        // Use default argument value
        try {
            CRef<CArgValue> arg_value(arg.ProcessDefault());
            // Add the value to "args"
            args.Add(arg_value);
        } 
        catch (CArgException&) {
            // mandatory argument, for CGI can be taken not only from the
            // command line but also from the HTTP request
            if (GetArgsType() != eCgiArgs  ||  caller == eConvertKeys) {
                throw;
            }
        }
    }
}
#endif


void CArgDescriptions::SetUsageContext
(const string& usage_name,
 const string& usage_description,
 bool          usage_sort_args,
 SIZE_TYPE     usage_width)
{
    m_UsageName        = usage_name;
    m_UsageDescription = usage_description;
    usage_sort_args ? SetMiscFlags(fUsageSortArgs) : ResetMiscFlags(fUsageSortArgs);

    const SIZE_TYPE kMinUsageWidth = 30;
    if (usage_width < kMinUsageWidth) {
        usage_width = kMinUsageWidth;
        //ERR_POST_X(23, Warning <<
        //               "CArgDescriptions::SetUsageContext() -- usage_width=" <<
        //               usage_width << " adjusted to " << kMinUsageWidth);
    }
    m_UsageWidth = usage_width;
}

void CArgDescriptions::SetDetailedDescription( const string& usage_description)
{
    m_DetailedDescription = usage_description;
}


void CArgDescriptions::x_AddDesc(CArgDesc& arg)
{
    const string& name = arg.GetName();

    if ( Exist(name) ) {
        //NCBI_THROW(CArgException,eSynopsis,
        string err_msg = (
            "Argument with this name is already defined: " + name);
        HBN_ERR("%s", err_msg.c_str());
    }

    arg.SetGroup(m_CurrentGroup);

    if (s_IsKey(arg)  ||  s_IsFlag(arg)) {
        _ASSERT(find(m_KeyFlagArgs.begin(), m_KeyFlagArgs.end(), name)
                == m_KeyFlagArgs.end());
        m_KeyFlagArgs.push_back(name);
    } else if ( !s_IsAlias(arg)  &&  !name.empty() ) {
        TPosArgs& container = s_IsOpening(arg) ? m_OpeningArgs : m_PosArgs;
        _ASSERT(find(container.begin(), container.end(), name)
                == container.end());
        if ( s_IsOptional(arg) ) {
            container.push_back(name);
        } else {
            TPosArgs::iterator it;
            for (it = container.begin();  it != container.end();  ++it) {
                if ( s_IsOptional(**x_Find(*it)) )
                    break;
            }
            container.insert(it, name);
        }
    }
    
    if ((arg.GetFlags() & fOptionalSeparator) != 0  &&
        name.length() == 1  &&
        s_IsKey(arg)) {
        m_NoSeparator += arg.GetName();
    }

    //m_Args.insert(&arg);
    m_Args.insert(unique_ptr<CArgDesc>(&arg));
}


void CArgDescriptions::PrintUsageIfNoArgs(bool do_print)
{
    do_print ? SetMiscFlags(fUsageIfNoArgs) : ResetMiscFlags(fUsageIfNoArgs);
}



///////////////////////////////////////////////////////
//  CArgDescriptions::PrintUsage()


static void s_PrintCommentBody(list<string>& arr, const string& s,
                               SIZE_TYPE width)
{
    NStr::Wrap(s, width, arr, NStr::fWrap_Hyphenate, "   ");
}


void CArgDescriptions::x_PrintComment(list<string>&   arr,
                                      const CArgDesc& arg,
                                      SIZE_TYPE       width) const
{
    string intro = ' ' + arg.GetUsageSynopsis(true/*name_only*/);

    // Print type (and value constraint, if any)
    string attr = arg.GetUsageCommentAttr();
    if ( !attr.empty() ) {
        char separator =
            (arg.GetFlags() & CArgDescriptions::fMandatorySeparator) ? '=' : ' ';
        string t;
        t += separator;
        t += '<' + attr + '>';
        if (arg.GetFlags() &  CArgDescriptions::fConfidential) {
            arr.push_back( intro + "  - read value interactively from console");
            arr.push_back( intro + "-file <" +
                           CArgDescriptions::GetTypeName(CArgDescriptions::eInputFile) + "> - read value from file");
            t = "-verbatim";
            t += separator;
            t += '<' + attr + '>';
        }
        attr = t;
    }

    // Add aliases for non-positional arguments
    list<string> negatives;
    if ( !s_IsPositional(arg) ) {
        ITERATE(CArgDescriptions::TArgs, it, m_Args) {
            const CArgDesc_Alias* alias =
                dynamic_cast<const CArgDesc_Alias*>(it->get());
            if (!alias  ||  alias->GetAliasedName() != arg.GetName()) {
                continue;
            }
            if ( alias->GetNegativeFlag() ) {
                negatives.push_back(alias->GetName());
            }
            else {
                intro += ", -" + alias->GetName();
            }
        }
    }

    intro += attr;

    // Wrap intro if necessary...
    {{
        SIZE_TYPE indent = intro.find(", ");
        if (indent == NPOS  ||  indent > width / 2) {
            indent = intro.find(" <");
            if (indent == NPOS  ||  indent > width / 2) {
                indent = 0;
            }
        }
        NStr::Wrap(intro, width, arr, NStr::fWrap_Hyphenate,
                   string(indent + 2, ' '), kEmptyStr);
    }}

    // Print description
    s_PrintCommentBody(arr, arg.GetComment(), width);

    // Print default value, if any
    const CArgDescDefault* dflt = dynamic_cast<const CArgDescDefault*> (&arg);
    if ( dflt ) {
        s_PrintCommentBody
            (arr, "Default = `" + dflt->GetDisplayValue() + '\'', width);
    }

    // Print required/excluded args
    string require;
    string exclude;
    pair<TDependency_CI, TDependency_CI> dep_rg =
        m_Dependencies.equal_range(arg.GetName());
    for (TDependency_CI dep = dep_rg.first; dep != dep_rg.second; ++dep) {
        switch ( dep->second.m_Dep ) {
        case eRequires:
            if ( !require.empty() ) {
                require += ", ";
            }
            require += dep->second.m_Arg;
            break;
        case eExcludes:
            if ( !exclude.empty() ) {
                exclude += ", ";
            }
            exclude += dep->second.m_Arg;
            break;
        }
    }
    if ( !require.empty() ) {
        s_PrintCommentBody(arr, " * Requires:  " + require, width);
    }
    if ( !exclude.empty() ) {
        s_PrintCommentBody(arr, " * Incompatible with:  " + exclude, width);
    }
    if ( !negatives.empty() ) {
        string neg_info;
        ITERATE(list<string>, neg, negatives) {
            if ( !neg_info.empty() ) {
                neg_info += ", ";
            }
            neg_info += *neg;
        }
        SIZE_TYPE indent = neg_info.find(", ");
        if (indent == NPOS  ||  indent > width / 2) {
            indent = 0;
        }
        neg_info = " -" + neg_info;
        NStr::Wrap(neg_info, width, arr, NStr::fWrap_Hyphenate,
                string(indent + 2, ' '), kEmptyStr);

        // Print description
        string neg_comment = arg.GetComment();
        if ( neg_comment.empty() ) {
            neg_comment = "Negative for " + arg.GetName();
        }
        s_PrintCommentBody(arr, neg_comment, width);
    }
    if (s_IsFlag(arg)) {
        const CArgDesc_Flag* fl = dynamic_cast<const CArgDesc_Flag*>(&arg);
        if (fl && !fl->GetSetValue()) {
            s_PrintCommentBody(arr, "When the flag is present, its value is FALSE", width);
        }
    }
}

CArgDescriptions::CPrintUsage::CPrintUsage(const CArgDescriptions& desc)
    : m_desc(desc)
{
    typedef list<const CArgDesc*> TList;
    typedef TList::iterator       TListI;

    m_args.push_front(0);
    TListI it_pos = m_args.begin();

    // Opening
    for (TPosArgs::const_iterator name = desc.m_OpeningArgs.begin();
         name != desc.m_OpeningArgs.end();  ++name) {
        TArgsCI it = desc.x_Find(*name);
        _ASSERT(it != desc.m_Args.end());
        if (it->get()->GetFlags() & CArgDescriptions::fHidden)
        {
            continue;
        }
        m_args.insert(it_pos, it->get());
    }

    // Keys and Flags
    if ( desc.IsSetMiscFlag(fUsageSortArgs) ) {
        // Alphabetically ordered,
        // mandatory keys to go first, then flags, then optional keys
        TListI& it_opt_keys = it_pos;
        TListI it_keys  = m_args.insert(it_pos,nullptr);
        TListI it_flags = m_args.insert(it_pos,nullptr);

        for (TArgsCI it = desc.m_Args.begin();  it != desc.m_Args.end();  ++it) {
            const CArgDesc* arg = it->get();
            if (it->get()->GetFlags() & CArgDescriptions::fHidden)
            {
                continue;
            }

            if (dynamic_cast<const CArgDesc_KeyOpt*> (arg)  ||
                dynamic_cast<const CArgDesc_KeyDef*> (arg)) {
                m_args.insert(it_opt_keys, arg);
            } else if (dynamic_cast<const CArgDesc_Key*> (arg)) {
                m_args.insert(it_keys, arg);
            } else if (dynamic_cast<const CArgDesc_Flag*> (arg)) {
                if ((desc.m_AutoHelp &&
                    NStr::strcmp(kArgHelp.c_str(),     (arg->GetName()).c_str()) == 0) ||
                    NStr::strcmp(kArgFullHelp.c_str(), (arg->GetName()).c_str()) == 0)
                    m_args.push_front(arg);
                else
                    m_args.insert(it_flags, arg);
            }
        }
        m_args.erase(it_keys);
        m_args.erase(it_flags);
    } else {
        // Unsorted, just the order they were described by user
        for (TKeyFlagArgs::const_iterator name = desc.m_KeyFlagArgs.begin();
             name != desc.m_KeyFlagArgs.end();  ++name) {
            TArgsCI it = desc.x_Find(*name);
            _ASSERT(it != desc.m_Args.end());
            if (it->get()->GetFlags() & CArgDescriptions::fHidden)
            {
                continue;
            }

            m_args.insert(it_pos, it->get());
        }
    }

    // Positional
    for (TPosArgs::const_iterator name = desc.m_PosArgs.begin();
         name != desc.m_PosArgs.end();  ++name) {
        TArgsCI it = desc.x_Find(*name);
        _ASSERT(it != desc.m_Args.end());
        if (it->get()->GetFlags() & CArgDescriptions::fHidden)
        {
            continue;
        }
        const CArgDesc* arg = it->get();

        // Mandatory args to go first, then go optional ones
        if (dynamic_cast<const CArgDesc_PosOpt*> (arg)) {
            m_args.push_back(arg);
        } else if (dynamic_cast<const CArgDesc_Pos*> (arg)) {
            m_args.insert(it_pos, arg);
        }
    }
    m_args.erase(it_pos);

    // Extra
    {{
        TArgsCI it = desc.x_Find(kEmptyStr);
        if (it != desc.m_Args.end()) {
            if ((it->get()->GetFlags() & CArgDescriptions::fHidden) == 0)
            {
                m_args.push_back(it->get());
            }
        }
    }}
}

CArgDescriptions::CPrintUsage::~CPrintUsage()
{
}

void CArgDescriptions::CPrintUsage::AddSynopsis(list<string>& arr,
    const string& intro, const string& prefix) const
{
    list<const CArgDesc*>::const_iterator it;
    list<string> syn;
    if (m_desc.GetArgsType() == eCgiArgs) {
        for (it = m_args.begin();  it != m_args.end();  ++it) {
            const CArgDescSynopsis* as = 
                dynamic_cast<const CArgDescSynopsis*>(&**it);

            if (as) {
                const string& name  = (*it)->GetName();
                const string& synopsis  = as->GetSynopsis();
                syn.push_back(name+"="+synopsis);
            }
        } // for
        NStr::WrapList(
            syn, m_desc.m_UsageWidth, "&", arr, 0, "?", "  "+m_desc.m_UsageName+"?");

    } else { // regular application
        if (!intro.empty()) {
            syn.push_back(intro);
        }
        for (it = m_args.begin();  it != m_args.end();  ++it) {
            if ( s_IsOptional(**it) || s_IsFlag(**it) ) {
                syn.push_back('[' + (*it)->GetUsageSynopsis() + ']');
            } else if ( s_IsPositional(**it) || s_IsOpening(**it) ) {
                syn.push_back('<' + (*it)->GetUsageSynopsis() + '>');
            } else {
                syn.push_back((*it)->GetUsageSynopsis());
            }
        } // for
        NStr::WrapList(syn, m_desc.m_UsageWidth, " ", arr, 0, prefix, "  ");
    }
}

void CArgDescriptions::CPrintUsage::AddDescription(list<string>& arr, bool detailed) const
{
    if ( m_desc.m_UsageDescription.empty() ) {
        arr.push_back("DESCRIPTION    -- none");
    } else {
        arr.push_back("DESCRIPTION");
        s_PrintCommentBody(arr,
            (detailed  && !m_desc.m_DetailedDescription.empty()) ?
                m_desc.m_DetailedDescription : m_desc.m_UsageDescription,
            m_desc.m_UsageWidth);
    }
}

void CArgDescriptions::CPrintUsage::AddCommandDescription(list<string>& arr,
    const string& cmd, const map<string,string>* aliases,
    size_t max_cmd_len, bool detailed) const
{
    if (detailed) {
        arr.push_back(kEmptyStr);
    }
    string cmd_full(cmd);
    if (aliases) {
        map<string,string>::const_iterator a = aliases->find(cmd);
        if (a != aliases->end()) {
            cmd_full += " (" + a->second + ")";
        }
    }
    cmd_full += string( max_cmd_len - cmd_full.size(), ' ');
    cmd_full += "- ";
    cmd_full += m_desc.m_UsageDescription;
    arr.push_back(string("  ")+ cmd_full);
    if (detailed) {
        AddSynopsis(arr,string(max_cmd_len+3,' '),string(max_cmd_len+6,' '));
    }
}

void CArgDescriptions::CPrintUsage::AddDetails(list<string>& arr) const
{
    list<const CArgDesc*>::const_iterator it;
    list<string> req;
    list<string> opt;
    // Collect mandatory args
    for (it = m_args.begin();  it != m_args.end();  ++it) {
        if (s_IsOptional(**it)  ||  s_IsFlag(**it)) {
            continue;
        }
        m_desc.x_PrintComment(req, **it, m_desc.m_UsageWidth);
    }
    // Collect optional args
    for (size_t grp = 0;  grp < m_desc.m_ArgGroups.size();  ++grp) {
        list<string> grp_opt;
        bool group_not_empty = false;
        if ( !m_desc.m_ArgGroups[grp].empty() ) {
            NStr::Wrap(m_desc.m_ArgGroups[grp], m_desc.m_UsageWidth, grp_opt,
                NStr::fWrap_Hyphenate, " *** ");
        }
        for (it = m_args.begin();  it != m_args.end();  ++it) {
            if (!s_IsOptional(**it)  &&  !s_IsFlag(**it)) {
                continue;
            }
            if ((*it)->GetGroup() == grp) {
                m_desc.x_PrintComment(grp_opt, **it, m_desc.m_UsageWidth);
                group_not_empty = true;
            }
        }
        if ( group_not_empty ) {
            opt.insert(opt.end(), grp_opt.begin(), grp_opt.end());
            opt.push_back(kEmptyStr);
        }
    }
    if ( !req.empty() ) {
        arr.push_back(kEmptyStr);
        arr.push_back("REQUIRED ARGUMENTS");
        arr.splice(arr.end(), req);
    }
    if ( !m_desc.m_nExtra  &&  !opt.empty() ) {
        arr.push_back(kEmptyStr);
        arr.push_back("OPTIONAL ARGUMENTS");
        arr.splice(arr.end(), opt);
    }

    // # of extra arguments
    if (m_desc.m_nExtra  ||  (m_desc.m_nExtraOpt != 0  && m_desc.m_nExtraOpt != kMax_UInt)) {
        string str_extra = "NOTE:  Specify ";
        if ( m_desc.m_nExtra ) {
            str_extra += "at least ";
            str_extra += NStr::UIntToString(m_desc.m_nExtra);
            if (m_desc.m_nExtraOpt != kMax_UInt) {
                str_extra += ", and ";
            }
        }
        if (m_desc.m_nExtraOpt != kMax_UInt) {
            str_extra += "no more than ";
            str_extra += NStr::UIntToString(m_desc.m_nExtra + m_desc.m_nExtraOpt);
        }
        str_extra +=
            " argument" + string(&"s"[m_desc.m_nExtra + m_desc.m_nExtraOpt == 1]) +
            " in \"....\"";
        s_PrintCommentBody(arr, str_extra, m_desc.m_UsageWidth);
    }
    if ( m_desc.m_nExtra  &&  !opt.empty() ) {
        arr.push_back(kEmptyStr);
        arr.push_back("OPTIONAL ARGUMENTS");
        arr.splice(arr.end(), opt);
    }

    if (!m_desc.m_DependencyGroups.empty()) {
        arr.push_back(kEmptyStr);
        arr.push_back("DEPENDENCY GROUPS");
        for (set< CConstRef<CArgDependencyGroup> >::const_iterator i = m_desc.m_DependencyGroups.begin();
            i != m_desc.m_DependencyGroups.end(); ++i) {
            //i->GetPointer()->PrintUsage(arr, 0);
            i->get()->PrintUsage(arr, 0);
        }
    }
}

string& CArgDescriptions::PrintUsage(string& str, bool detailed) const
{
    CPrintUsage x(*this);
    list<string> arr;

    // SYNOPSIS
    arr.push_back("USAGE");
    x.AddSynopsis(arr, m_UsageName,"    ");

    // DESCRIPTION
    arr.push_back(kEmptyStr);
    x.AddDescription(arr, detailed);

    // details
    if (detailed) {
        x.AddDetails(arr);
    } else {
        arr.push_back(kEmptyStr);
        arr.push_back("Use '-help' to print detailed descriptions of command line arguments");
    }

    str += NStr::Join(arr, "\n");
    str += "\n";
    return str;
}

string& CArgDescriptions::HbnPrintUsage(const string& main_usage, string& str, bool detailed) const
{
    CPrintUsage x(*this);
    list<string> arr;

    arr.push_back("OPTIONAL ARGUMENTS");
    x.AddSynopsis(arr, kEmptyStr, "");

    // DESCRIPTION
    arr.push_back(kEmptyStr);
    x.AddDescription(arr, detailed);

    // details
    if (detailed) {
        x.AddDetails(arr);
    } else {
        arr.push_back(kEmptyStr);
        arr.push_back("Use '-help' to print detailed descriptions of command line arguments");
    }

    str = "USAGE:\n  " + main_usage + "\n\n";
    str += NStr::Join(arr, "\n");
    str += "\n";
    return str;
}

CArgDescriptions::CPrintUsageXml::CPrintUsageXml(const CArgDescriptions& desc, CNcbiOstream& out)
    : m_desc(desc), m_out(out)
{
#ifndef HBN_REMOVE_THIS
    m_out << "<?xml version=\"1.0\"?>" << endl;
    m_out << "<" << "ncbi_application xmlns=\"ncbi:application\"" << endl
        << " xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"" << endl
        << " xs:schemaLocation=\"ncbi:application ncbi_application.xsd\"" << endl
        << ">" << endl;
    m_out << "<" << "program" << " type=\"";
    if (desc.GetArgsType() == eRegularArgs) {
        m_out << "regular";
    } else if (desc.GetArgsType() == eCgiArgs) {
        m_out << "cgi";
    } else {
        m_out << "UNKNOWN";
    }
    m_out << "\"" << ">" << endl;    
    s_WriteXmlLine(m_out, "name", desc.m_UsageName);
    s_WriteXmlLine(m_out, "version", 
        CNcbiApplication::Instance()->GetVersion().Print());
    s_WriteXmlLine(m_out, "description", desc.m_UsageDescription);
    s_WriteXmlLine(m_out, "detailed_description", desc.m_DetailedDescription);
    m_out << "</" << "program" << ">" << endl;
#endif
}

CArgDescriptions::CPrintUsageXml::~CPrintUsageXml()
{
    m_out << "</" << "ncbi_application" << ">" << endl;
}

void CArgDescriptions::CPrintUsageXml::PrintArguments(const CArgDescriptions& desc) const
{
#ifndef HBN_REMOVE_THIS
    m_out << "<" << "arguments";
    if (desc.GetPositionalMode() == ePositionalMode_Loose) {
        m_out << " positional_mode=\"loose\"";
    }
    m_out << ">" << endl;

    string tag;

// opening
    ITERATE(TPosArgs, p, desc.m_OpeningArgs) {
        ITERATE (TArgs, a, desc.m_Args) {
            if ((**a).GetName() == *p) {
                tag = (*a)->PrintXml(m_out);
                m_out << "</" << tag << ">" << endl;
            }
        }
    }
// positional
    ITERATE(TPosArgs, p, desc.m_PosArgs) {
        ITERATE (TArgs, a, desc.m_Args) {
            if ((**a).GetName() == *p) {
                tag = (*a)->PrintXml(m_out);
                desc.x_PrintAliasesAsXml(m_out, (*a)->GetName());
                m_out << "</" << tag << ">" << endl;
            }
        }
    }
// keys
    ITERATE (TArgs, a, desc.m_Args) {
        if (s_IsKey(**a)) {
            tag = (*a)->PrintXml(m_out);
            desc.x_PrintAliasesAsXml(m_out, (*a)->GetName());
            m_out << "</" << tag << ">" << endl;
        }
    }
// flags
    ITERATE (TArgs, a, desc.m_Args) {
        if (s_IsFlag(**a)) {
            tag = (*a)->PrintXml(m_out);
            desc.x_PrintAliasesAsXml(m_out, (*a)->GetName());
            desc.x_PrintAliasesAsXml(m_out, (*a)->GetName(), true);
            m_out << "</" << tag << ">" << endl;
        }
    }
// extra positional
    ITERATE (TArgs, a, desc.m_Args) {
        if (s_IsPositional(**a) && (**a).GetName().empty()) {
            tag = (*a)->PrintXml(m_out);
            s_WriteXmlLine(m_out, "min_occurs", NStr::UIntToString(desc.m_nExtra));
            s_WriteXmlLine(m_out, "max_occurs", NStr::UIntToString(desc.m_nExtraOpt));
            m_out << "</" << tag << ">" << endl;
        }
    }
    if (!desc.m_Dependencies.empty()) {
        m_out << "<" << "dependencies" << ">" << endl;
        ITERATE(TDependencies, dep, desc.m_Dependencies) {
            if (dep->second.m_Dep == eRequires) {
                m_out << "<" << "first_requires_second" << ">" << endl;
                s_WriteXmlLine(m_out, "arg1", dep->first);
                s_WriteXmlLine(m_out, "arg2", dep->second.m_Arg);
                m_out << "</" << "first_requires_second" << ">" << endl;
            }
        }
        ITERATE(TDependencies, dep, desc.m_Dependencies) {
            if (dep->second.m_Dep == eExcludes) {
                m_out << "<" << "first_excludes_second" << ">" << endl;
                s_WriteXmlLine(m_out, "arg1", dep->first);
                s_WriteXmlLine(m_out, "arg2", dep->second.m_Arg);
                m_out << "</" << "first_excludes_second" << ">" << endl;
            }
        }
        m_out << "</" << "dependencies" << ">" << endl;
    }

    for (set< CConstRef<CArgDependencyGroup> >::const_iterator i = m_desc.m_DependencyGroups.begin();
        i != m_desc.m_DependencyGroups.end(); ++i) {
        i->GetPointer()->PrintUsageXml(m_out);
    }
    m_out << "</" << "arguments" << ">" << endl;
#endif
}

void CArgDescriptions::PrintUsageXml(CNcbiOstream& out) const
{
    CPrintUsageXml x(*this,out);
    x.PrintArguments(*this);
}

#ifndef HBN_REMOVE_THIS
void CArgDescriptions::x_PrintAliasesAsXml( CNcbiOstream& out,
    const string& name, bool negated /* =false*/) const
{
    ITERATE (TArgs, a, m_Args) {
        if (s_IsAlias(**a)) {
            const CArgDesc_Alias& alias =
                dynamic_cast<const CArgDesc_Alias&>(**a);
            if (negated == alias.GetNegativeFlag()) {
                string tag = negated ? "negated_alias" : "alias";
                if (alias.GetAliasedName() == name) {
                    s_WriteXmlLine(out, tag, alias.GetName());
                }
            }
        }
    }
}
#endif

/////////////////////////////////////////////////////////////////////////////

CRef<CArgDependencyGroup> CArgDependencyGroup::Create(
        const string& name, const string& description)
{
    CRef<CArgDependencyGroup> gr(new CArgDependencyGroup());
    gr->m_Name = name;
    gr->m_Description = description;
    return gr;
}

CArgDependencyGroup::CArgDependencyGroup()
    : m_MinMembers(0), m_MaxMembers(0)
{
}

CArgDependencyGroup::~CArgDependencyGroup(void)
{
}

CArgDependencyGroup& CArgDependencyGroup::SetMinMembers(size_t min_members)
{
    m_MinMembers = min_members;
    return *this;
}

CArgDependencyGroup& CArgDependencyGroup::SetMaxMembers(size_t max_members)
{
    m_MaxMembers = max_members;
    return *this;
}

CArgDependencyGroup& CArgDependencyGroup::Add(const string& arg_name, EInstantSet  instant_set)
{
    m_Arguments[arg_name] = instant_set;
    return *this;
}

CArgDependencyGroup& CArgDependencyGroup::Add(
    CArgDependencyGroup* dep_group, EInstantSet instant_set)
{
    m_Groups[ CConstRef<CArgDependencyGroup>(dep_group)] = instant_set;
    return *this;
}

void CArgDependencyGroup::Evaluate( const CArgs& args) const
{
    x_Evaluate(args, nullptr, nullptr);
}

bool CArgDependencyGroup::x_Evaluate( const CArgs& args, string* arg_set, string* arg_unset) const
{
    bool top_level = !arg_set || !arg_unset;
    bool has_instant_set = false;
    size_t count_set = 0;
    set<string> names_set, names_unset;
    string args_set, args_unset;

    for (map< CConstRef<CArgDependencyGroup>, EInstantSet>::const_iterator i = m_Groups.begin();
        i != m_Groups.end(); ++i) {
        string msg_set, msg_unset;
        if (i->first.get()->x_Evaluate(args, &msg_set, &msg_unset)) {
            ++count_set;
            has_instant_set = has_instant_set || (i->second == eInstantSet);
            names_set.insert(msg_set);
        } else {
            names_unset.insert(msg_unset);
        }
    }
    for (map<string, EInstantSet>::const_iterator i = m_Arguments.begin();
        i != m_Arguments.end(); ++i) {
        if (args.Exist(i->first)) {
            ++count_set;
            has_instant_set = has_instant_set || (i->second == eInstantSet);
            names_set.insert(i->first);
        } else {
            names_unset.insert(i->first);
        }
    }
    size_t count_total = m_Groups.size() + m_Arguments.size();
    size_t count_max = m_MaxMembers != 0 ? m_MaxMembers : count_total;

    if (names_set.size() > 1) {
        args_set = "(" + NStr::Join(names_set, ", ") + ")";
    } else if (names_set.size() == 1) {
        args_set = *names_set.begin();
    }

    if (names_unset.size() > 1) {
        args_unset = "(" + NStr::Join(names_unset, m_MinMembers <= 1 ? " | " : ", ") + ")";
    } else if (names_unset.size() == 1) {
        args_unset = *names_unset.begin();
    }

    bool result = count_set != 0 || top_level;
    if (result) {
        if (count_set > count_max) {
            string msg("Argument conflict: ");
            msg += args_set + " may not be specified simultaneously";
            //NCBI_THROW(CArgException, eConstraint, msg);
            HBN_ERR("%s", msg.c_str());
        }
        if (!has_instant_set && count_set < m_MinMembers) {
            string msg("Argument has no value: ");
            if (count_total != count_max) {
                msg += (m_MinMembers - count_set > 1) ? "some" : "one";                                  
                msg += " of ";
            }
            msg += args_unset + " must be specified";
            //NCBI_THROW(CArgException,eNoValue, msg);
            HBN_ERR("%s", msg.c_str());
        }
    }
    if (arg_set) {
        *arg_set = args_set;
    }
    if (arg_unset) {
        *arg_unset = args_unset;
    }
    return result;
}

void CArgDependencyGroup::PrintUsage(list<string>& arr, size_t offset) const
{
    arr.push_back(kEmptyStr);
    string off(2*offset,' ');
    string msg(off);
    msg += m_Name + ": {";

    bool first = true;
    list<string> instant;
    for (map< CConstRef<CArgDependencyGroup>, EInstantSet>::const_iterator i = m_Groups.begin();
        i != m_Groups.end(); ++i) {
        if (!first) {
            msg += ",";
        }
        first = false;
        msg += i->first.get()->m_Name;
        if (i->second == eInstantSet) {
            instant.push_back(i->first.get()->m_Name);
        }
    }
    for (map<string, EInstantSet>::const_iterator i = m_Arguments.begin();
        i != m_Arguments.end(); ++i) {
        if (!first) {
            msg += ",";
        }
        first = false;
        msg += i->first;
        if (i->second == eInstantSet) {
            instant.push_back(i->first);
        }
    }
    msg += "}";
    arr.push_back(msg);
    if (!m_Description.empty()) {
        msg = off;
        msg += m_Description;
        arr.push_back(msg);
    }
    size_t count_total = m_Groups.size() + m_Arguments.size();
    size_t count_max = m_MaxMembers != 0 ? m_MaxMembers : count_total;

    msg = off + "in which ";
    size_t count = m_MinMembers;
    if (m_MinMembers == count_max) {
        msg += "exactly ";
        //msg += NStr::NumericToString(m_MinMembers);
        msg += NStr::UInt8ToString(m_MinMembers);
    } else if (count_max == count_total && m_MinMembers != 0) {
        msg += "at least ";
        //msg += NStr::NumericToString(m_MinMembers);
        msg += NStr::UInt8ToString(m_MinMembers);
    } else if (count_max != count_total && m_MinMembers == 0) {
        msg += "no more than ";
        //msg += NStr::NumericToString(m_MaxMembers);
        msg += NStr::UInt8ToString(m_MaxMembers);
        count = m_MaxMembers;
    } else {
        //msg += NStr::NumericToString(m_MinMembers);
        msg += NStr::UInt8ToString(m_MinMembers);
        msg += " to ";
        //msg += NStr::NumericToString(m_MaxMembers);
        msg += NStr::UInt8ToString(m_MaxMembers);
        count = m_MaxMembers;
    }
    msg += " element";
    if (count != 1) {
        msg += "s";
    }
    msg += " must be set";
    arr.push_back(msg);

    if (!instant.empty()) {
        msg = off;
        msg += "Instant set: ";
        msg += NStr::Join(instant, ",");
        arr.push_back(msg);
    }
    for (map< CConstRef<CArgDependencyGroup>, EInstantSet>::const_iterator i = m_Groups.begin();
        i != m_Groups.end(); ++i) {
        i->first.get()->PrintUsage(arr, offset+1);
    }
}

void CArgDependencyGroup::PrintUsageXml(CNcbiOstream& out) const
{
    out << "<" << "dependencygroup" << ">" << endl;
    out << "<" << "name" << ">" << m_Name << "</" << "name" << ">" << endl;
    out << "<" << "description" << ">" << m_Description << "</" << "description" << ">" << endl;

    for (map< CConstRef<CArgDependencyGroup>, EInstantSet>::const_iterator i = m_Groups.begin();
        i != m_Groups.end(); ++i) {
        out << "<" << "group";
        if (i->second == eInstantSet) {
            out << " instantset=\"true\"";
        }
        out << ">" << i->first.get()->m_Name << "</" << "group" << ">" << endl;
    }
    for (map<string, EInstantSet>::const_iterator i = m_Arguments.begin();
        i != m_Arguments.end(); ++i) {
        out << "<" << "argument";
        if (i->second == eInstantSet) {
            out << " instantset=\"true\"";
        }
        out << ">" << i->first << "</" << "argument" << ">" << endl;
    }
    out << "<" << "minmembers" << ">" << m_MinMembers << "</" << "minmembers" << ">" << endl;
    out << "<" << "maxmembers" << ">" << m_MaxMembers << "</" << "maxmembers" << ">" << endl;
    for (map< CConstRef<CArgDependencyGroup>, EInstantSet>::const_iterator i = m_Groups.begin();
        i != m_Groups.end(); ++i) {
        i->first.get()->PrintUsageXml(out);
    }
    out << "</" << "dependencygroup" << ">" << endl;
}

END_NCBI_SCOPE