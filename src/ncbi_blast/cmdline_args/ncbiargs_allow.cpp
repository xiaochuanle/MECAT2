#include "ncbiargs_allow.hpp"

#include "../str_util/ncbistr.hpp"

BEGIN_NCBI_SCOPE

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
// CArgAllow::
//   CArgAllow_Symbols::
//   CArgAllow_String::
//   CArgAllow_Strings::
//   CArgAllow_Int8s::
//   CArgAllow_Integers::
//   CArgAllow_Doubles::
//


///////////////////////////////////////////////////////
//  CArgAllow::
//

CArgAllow::~CArgAllow(void)
{
}

void CArgAllow::PrintUsageXml(CNcbiOstream& ) const
{
}

CArgAllow* CArgAllow::Clone(void) const
{
    return NULL;
}

///////////////////////////////////////////////////////
//  s_IsSymbol() -- check if the symbol belongs to one of standard character
//                  classes from <ctype.h>, or to user-defined symbol set
//

inline bool s_IsAllowedSymbol(unsigned char                   ch,
                              CArgAllow_Symbols::ESymbolClass symbol_class,
                              const string&                   symbol_set)
{
    switch ( symbol_class ) {
    case CArgAllow_Symbols::eAlnum:   return isalnum(ch) != 0;
    case CArgAllow_Symbols::eAlpha:   return isalpha(ch) != 0;
    case CArgAllow_Symbols::eCntrl:   return iscntrl(ch) != 0;
    case CArgAllow_Symbols::eDigit:   return isdigit(ch) != 0;
    case CArgAllow_Symbols::eGraph:   return isgraph(ch) != 0;
    case CArgAllow_Symbols::eLower:   return islower(ch) != 0;
    case CArgAllow_Symbols::ePrint:   return isprint(ch) != 0;
    case CArgAllow_Symbols::ePunct:   return ispunct(ch) != 0;
    case CArgAllow_Symbols::eSpace:   return isspace(ch) != 0;
    case CArgAllow_Symbols::eUpper:   return isupper(ch) != 0;
    case CArgAllow_Symbols::eXdigit:  return isxdigit(ch) != 0;
    case CArgAllow_Symbols::eUser:
        return symbol_set.find_first_of(ch) != NPOS;
    }
    _TROUBLE;  return false;
}


static string s_GetUsageSymbol(CArgAllow_Symbols::ESymbolClass symbol_class,
                               const string&                   symbol_set)
{
    switch ( symbol_class ) {
    case CArgAllow_Symbols::eAlnum:   return "alphanumeric";
    case CArgAllow_Symbols::eAlpha:   return "alphabetic";
    case CArgAllow_Symbols::eCntrl:   return "control symbol";
    case CArgAllow_Symbols::eDigit:   return "decimal";
    case CArgAllow_Symbols::eGraph:   return "graphical symbol";
    case CArgAllow_Symbols::eLower:   return "lower case";
    case CArgAllow_Symbols::ePrint:   return "printable";
    case CArgAllow_Symbols::ePunct:   return "punctuation";
    case CArgAllow_Symbols::eSpace:   return "space";
    case CArgAllow_Symbols::eUpper:   return "upper case";
    case CArgAllow_Symbols::eXdigit:  return "hexadecimal";
    case CArgAllow_Symbols::eUser:
        return "'" + NStr::PrintableString(symbol_set) + "'";
    }
    _TROUBLE;  return kEmptyStr;
}

static string s_GetSymbolClass(CArgAllow_Symbols::ESymbolClass symbol_class)
{
    switch ( symbol_class ) {
    case CArgAllow_Symbols::eAlnum:   return "Alnum";
    case CArgAllow_Symbols::eAlpha:   return "Alpha";
    case CArgAllow_Symbols::eCntrl:   return "Cntrl";
    case CArgAllow_Symbols::eDigit:   return "Digit";
    case CArgAllow_Symbols::eGraph:   return "Graph";
    case CArgAllow_Symbols::eLower:   return "Lower";
    case CArgAllow_Symbols::ePrint:   return "Print";
    case CArgAllow_Symbols::ePunct:   return "Punct";
    case CArgAllow_Symbols::eSpace:   return "Space";
    case CArgAllow_Symbols::eUpper:   return "Upper";
    case CArgAllow_Symbols::eXdigit:  return "Xdigit";
    case CArgAllow_Symbols::eUser:    return "User";
    }
    _TROUBLE;  return kEmptyStr;
}



///////////////////////////////////////////////////////
//  CArgAllow_Symbols::
//

CArgAllow_Symbols::CArgAllow_Symbols(ESymbolClass symbol_class)
    : CArgAllow()
{
    Allow( symbol_class);
    return;
}


CArgAllow_Symbols::CArgAllow_Symbols(const string& symbol_set)
    : CArgAllow()
{
    Allow( symbol_set);
    return;
}


bool CArgAllow_Symbols::Verify(const string& value) const
{
    if (value.length() != 1)
        return false;

    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        if (s_IsAllowedSymbol(value[0], pi->first, pi->second)) {
            return true;
        }
    }
    return false;
}


string CArgAllow_Symbols::GetUsage(void) const
{
    string usage;
    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        if (!usage.empty()) {
            usage += ", or ";
        }
        usage += s_GetUsageSymbol(pi->first, pi->second);
    }

    return "one symbol: " + usage;
}

void CArgAllow_Symbols::PrintUsageXml(CNcbiOstream& out) const
{
#ifndef HBN_REMOVE_THIS
    out << "<" << "Symbols" << ">" << endl;
    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        if (pi->first != eUser) {
            s_WriteXmlLine( out, "type", s_GetSymbolClass(pi->first).c_str());
        } else {
            ITERATE( string, p, pi->second) {
                string c;
                s_WriteXmlLine( out, "value", c.append(1,*p).c_str());
            }
        }
    }
    out << "</" << "Symbols" << ">" << endl;
#endif
}

CArgAllow_Symbols&
CArgAllow_Symbols::Allow(CArgAllow_Symbols::ESymbolClass symbol_class)
{
    m_SymClass.insert( make_pair(symbol_class, kEmptyStr) );
    return *this;
}

CArgAllow_Symbols& CArgAllow_Symbols::Allow(const string& symbol_set)
{
    m_SymClass.insert( make_pair(eUser, symbol_set ));
    return *this;
}

CArgAllow* CArgAllow_Symbols::Clone(void) const
{
    CArgAllow_Symbols* clone = new CArgAllow_Symbols;
    clone->m_SymClass = m_SymClass;
    return clone;
}



///////////////////////////////////////////////////////
//  CArgAllow_String::
//

CArgAllow_String::CArgAllow_String(ESymbolClass symbol_class)
    : CArgAllow_Symbols(symbol_class)
{
    return;
}


CArgAllow_String::CArgAllow_String(const string& symbol_set)
    : CArgAllow_Symbols(symbol_set)
{
    return;
}


bool CArgAllow_String::Verify(const string& value) const
{
    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        string::const_iterator it;
        for (it = value.begin();  it != value.end(); ++it) {
            if ( !s_IsAllowedSymbol(*it, pi->first, pi->second) )
                break;;
        }
        if (it == value.end()) {
            return true;
        }
    }
    return false;
}


string CArgAllow_String::GetUsage(void) const
{
    string usage;
    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        if (!usage.empty()) {
            usage += ", or ";
        }
        usage += s_GetUsageSymbol(pi->first, pi->second);
    }

    return "to contain only symbols: " + usage;
}


void CArgAllow_String::PrintUsageXml(CNcbiOstream& out) const
{
#ifndef HBN_REMOVE_THIS
    out << "<" << "String" << ">" << endl;
    ITERATE( set< TSymClass >, pi,  m_SymClass) {
        if (pi->first != eUser) {
            s_WriteXmlLine( out, "type", s_GetSymbolClass(pi->first).c_str());
        } else {
            s_WriteXmlLine( out, "charset", pi->second.c_str());
        }
    }
    out << "</" << "String" << ">" << endl;
#endif
}

CArgAllow* CArgAllow_String::Clone(void) const
{
    CArgAllow_String* clone = new CArgAllow_String;
    clone->m_SymClass = m_SymClass;
    return clone;
}



///////////////////////////////////////////////////////
//  CArgAllow_Strings::
//

CArgAllow_Strings::CArgAllow_Strings(NStr::ECase use_case)
    : CArgAllow(),
      m_Strings(PNocase_Conditional(use_case))
{
    return;
}


CArgAllow_Strings* CArgAllow_Strings::Allow(const string& value)
{
    m_Strings.insert(value);
    return this;
}


bool CArgAllow_Strings::Verify(const string& value) const
{
    TStrings::const_iterator it = m_Strings.find(value);
    return it != m_Strings.end();
}


string 
CArgAllow_Strings::GetUsage(void) const
{
    if ( m_Strings.empty() ) {
        return "ERROR:  Constraint with no values allowed(?!)";
    }

    string str;
    TStrings::const_iterator it = m_Strings.begin();
    for (;;) {
        str += "`";
        str += *it;

        ++it;
        if (it == m_Strings.end()) {
            str += "'";
            if ( m_Strings.key_comp()("a", "A") ) {
                str += "  {case insensitive}";
            }
            break;
        }
        str += "', ";
    }
    return str;
}


void CArgAllow_Strings::PrintUsageXml(CNcbiOstream& out) const
{
#ifndef HBN_REMOVE_THIS
    out << "<" << "Strings";
    out << " case_sensitive=\"";
    if ( m_Strings.key_comp()("a", "A") ) {
        out << "false";
    } else {
        out << "true";
    }
    out << "\">" << endl;
    ITERATE( TStrings, p, m_Strings) {
        s_WriteXmlLine( out, "value", (*p).c_str());
    }
    out << "</" << "Strings" << ">" << endl;
#endif
}


CArgAllow_Strings& CArgAllow_Strings::AllowValue(const string& value)
{
    return *Allow(value);
}

CArgAllow* CArgAllow_Strings::Clone(void) const
{
    CArgAllow_Strings* clone = new CArgAllow_Strings(m_Strings.key_comp().GetCase());
    clone->m_Strings = m_Strings;
    return clone;
}


///////////////////////////////////////////////////////
//  CArgAllow_Int8s::
//

CArgAllow_Int8s::CArgAllow_Int8s(Int8 x_)
    : CArgAllow()
{
    Allow( x_);
}

CArgAllow_Int8s::CArgAllow_Int8s(Int8 x_min, Int8 x_max)
    : CArgAllow()
{
    AllowRange(x_min, x_max);
}

extern Int8 s_StringToInt8(const string& value);

bool CArgAllow_Int8s::Verify(const string& value) const
{
    Int8 val = s_StringToInt8(value);
    ITERATE( set< TInterval >, pi, m_MinMax) {
        if (pi->first <= val && val<= pi->second) {
            return true;
        }
    }
    return false;
}


string CArgAllow_Int8s::GetUsage(void) const
{
    if (m_MinMax.size() == 1) {
        Int8 x_min = m_MinMax.begin()->first;
        Int8 x_max = m_MinMax.begin()->second;
        if (x_min == x_max) {
            return NStr::Int8ToString(x_min);
        } else if (x_min == kMin_I8 && x_max != kMax_I8) {
            return string("less or equal to ") + NStr::Int8ToString(x_max);
        } else if (x_min != kMin_I8 && x_max == kMax_I8) {
            return string("greater or equal to ") + NStr::Int8ToString(x_min);
        } else if (x_min == kMin_I8 && x_max == kMax_I8) {
            return kEmptyStr;
        }
    }
    string usage;
    ITERATE( set< TInterval >, pi, m_MinMax) {
        if (!usage.empty()) {
            usage += ", ";
        }
        if (pi->first == pi->second) {
            usage += NStr::Int8ToString(pi->first);
        } else {
            usage += NStr::Int8ToString(pi->first) + ".." + NStr::Int8ToString(pi->second);
        }

    }
    return usage;
}


void CArgAllow_Int8s::PrintUsageXml(CNcbiOstream& out) const
{
#ifndef HBN_REMOVE_THIS
    string tag("Int8s");
    if (dynamic_cast<const CArgAllow_Integers*>(this) != 0) {
        tag = "Integers";
    }
    out << "<" << tag << ">" << endl;
    ITERATE( set< TInterval >, pi, m_MinMax) {
        s_WriteXmlLine( out, "min", NStr::Int8ToString(pi->first).c_str());
        s_WriteXmlLine( out, "max", NStr::Int8ToString(pi->second).c_str());
    }
    out << "</" << tag << ">" << endl;
#endif
}

CArgAllow_Int8s& CArgAllow_Int8s::AllowRange(Int8 from, Int8 to)
{
    m_MinMax.insert( make_pair(from,to) );
    return *this;
}

CArgAllow_Int8s& CArgAllow_Int8s::Allow(Int8 value)
{
    m_MinMax.insert( make_pair(value,value) );
    return *this;
}

CArgAllow* CArgAllow_Int8s::Clone(void) const
{
    CArgAllow_Int8s* clone = new CArgAllow_Int8s;
    clone->m_MinMax = m_MinMax;
    return clone;
}



///////////////////////////////////////////////////////
//  CArgAllow_Integers::
//

CArgAllow_Integers::CArgAllow_Integers(int x_)
    : CArgAllow_Int8s(x_)
{
}

CArgAllow_Integers::CArgAllow_Integers(int x_min, int x_max)
    : CArgAllow_Int8s(x_min, x_max)
{
}

string CArgAllow_Integers::GetUsage(void) const
{
    if (m_MinMax.size() == 1) {
        Int8 x_min = m_MinMax.begin()->first;
        Int8 x_max = m_MinMax.begin()->second;
        if (x_min == x_max) {
            return NStr::Int8ToString(x_min);
        } else if (x_min == kMin_Int && x_max != kMax_Int) {
            return string("less or equal to ") + NStr::Int8ToString(x_max);
        } else if (x_min != kMin_Int && x_max == kMax_Int) {
            return string("greater or equal to ") + NStr::Int8ToString(x_min);
        } else if (x_min == kMin_Int && x_max == kMax_Int) {
            return kEmptyStr;
        }
    }
    return CArgAllow_Int8s::GetUsage();;
}

CArgAllow* CArgAllow_Integers::Clone(void) const
{
    CArgAllow_Integers* clone = new CArgAllow_Integers;
    clone->m_MinMax = m_MinMax;
    return clone;
}


///////////////////////////////////////////////////////
//  CArgAllow_Doubles::
//

CArgAllow_Doubles::CArgAllow_Doubles(double x_value)
    : CArgAllow()
{
    Allow(x_value);
}

CArgAllow_Doubles::CArgAllow_Doubles(double x_min, double x_max)
    : CArgAllow()
{
    AllowRange( x_min, x_max );
}


bool CArgAllow_Doubles::Verify(const string& value) const
{
    double val = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
    ITERATE( set< TInterval >, pi, m_MinMax) {
        if (pi->first <= val && val<= pi->second) {
            return true;
        }
    }
    return false;
}


string CArgAllow_Doubles::GetUsage(void) const
{
    if (m_MinMax.size() == 1) {
        double x_min = m_MinMax.begin()->first;
        double x_max = m_MinMax.begin()->second;
        if (x_min == x_max) {
            return NStr::DoubleToString(x_min);
        } else if (x_min == kMin_Double && x_max != kMax_Double) {
            return string("less or equal to ") + NStr::DoubleToString(x_max);
        } else if (x_min != kMin_Double && x_max == kMax_Double) {
            return string("greater or equal to ") + NStr::DoubleToString(x_min);
        } else if (x_min == kMin_Double && x_max == kMax_Double) {
            return kEmptyStr;
        }
    }
    string usage;
    ITERATE( set< TInterval >, pi, m_MinMax) {
        if (!usage.empty()) {
            usage += ", ";
        }
        if (pi->first == pi->second) {
            usage += NStr::DoubleToString(pi->first);
        } else {
            usage += NStr::DoubleToString(pi->first) + ".." + NStr::DoubleToString(pi->second);
        }

    }
    return usage;
}


void CArgAllow_Doubles::PrintUsageXml(CNcbiOstream& out) const
{
#ifndef HBN_REMOVE_THIS
    out << "<" << "Doubles" << ">" << endl;
    ITERATE( set< TInterval >, pi, m_MinMax) {
        s_WriteXmlLine( out, "min", NStr::DoubleToString(pi->first).c_str());
        s_WriteXmlLine( out, "max", NStr::DoubleToString(pi->second).c_str());
    }
    out << "</" << "Doubles" << ">" << endl;
#endif
}

CArgAllow_Doubles& CArgAllow_Doubles::AllowRange(double from, double to)
{
    m_MinMax.insert( make_pair(from,to) );
    return *this;
}

CArgAllow_Doubles& CArgAllow_Doubles::Allow(double value)
{
    m_MinMax.insert( make_pair(value,value) );
    return *this;
}

CArgAllow* CArgAllow_Doubles::Clone(void) const
{
    CArgAllow_Doubles* clone = new CArgAllow_Doubles;
    clone->m_MinMax = m_MinMax;
    return clone;
}

END_NCBI_SCOPE