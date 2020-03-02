#include "ncbistr_util.hpp"

#include "../../corelib/hbn_aux.h"

BEGIN_NCBI_SCOPE

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//  CTempString (deprecated constructors, defined out of line to cut down
//  on spurious warnings when building with compilers that warn on
//  definition rather than merely, and arguably more sensibly, on usage).



void CTempStringList::Join(string* s) const
{
    s->reserve(GetSize());
    *s = m_FirstNode.str;
    for (const SNode* node = m_FirstNode.next.get();  node != NULL;
         node = node->next.get()) {
        s->append(node->str.data(), node->str.size());
    }
}


void CTempStringList::Join(CTempString* s) const
{
    CTempStringEx str;
    Join(&str);
    *s = str;
}


void CTempStringList::Join(CTempStringEx* s) const
{
    if (m_FirstNode.next.get() == NULL) {
        *s = m_FirstNode.str;
    } else {
        if ( !m_Storage ) {
            //NCBI_THROW2(CStringException, eBadArgs,
            HBN_ERR("CTempStringList::Join(): non-NULL storage required");
        }
        SIZE_TYPE n = GetSize();
        char* buf = m_Storage->Allocate(n + 1);
        char* p = buf;
        for (const SNode* node = &m_FirstNode;  node != NULL;
             node = node->next.get()) {
            memcpy(p, node->str.data(), node->str.size());
            p += node->str.size();
        }
        *p = '\0';
        s->assign(buf, n);
    }
}


SIZE_TYPE CTempStringList::GetSize(void) const
{
    SIZE_TYPE total = m_FirstNode.str.size();
    for (const SNode* node = m_FirstNode.next.get();  node != NULL;
         node = node->next.get()) {
        total += node->str.size();
    }
    return total;
}


CTempString_Storage::CTempString_Storage(void)
{
}


CTempString_Storage::~CTempString_Storage(void)
{
    NON_CONST_ITERATE(TData, it, m_Data) {
        delete[] (*it);
        *it = 0;
    }
}


char* CTempString_Storage::Allocate(CTempString::size_type len)
{
    m_Data.push_back(new char[len]);
    return m_Data.back();
}

bool CStrTokenizeBase::Advance(CTempStringList* part_collector, SIZE_TYPE* ptr_part_start, SIZE_TYPE* ptr_delim_pos)
{
    SIZE_TYPE pos, part_start, delim_pos = 0, quote_pos = 0;
    bool      found_text = false, done;
    char      active_quote = '\0';

    // Skip leading delimiters.
    // NOTE: We cannot process 
    if (!m_Pos  &&  (m_Flags & NStr::fSplit_Truncate_Begin) != 0) {
        ;
    }
    pos = part_start = m_Pos;
    done = (pos == NPOS);
    // save part start position
    if (ptr_part_start) {
        *ptr_part_start = part_start;
    }

    // ChecksSkipDelims
    if (pos >= m_Str.size()) {
        pos = NPOS;
        done = true;
    }
    if (ptr_delim_pos) {
        *ptr_delim_pos = NPOS;
    }

    // Each chunk covers the half-open interval [part_start, delim_pos).

    while ( !done  &&
            ((delim_pos = m_Str.find_first_of(m_InternalDelim, pos)) != NPOS)) {

        SIZE_TYPE next_start = pos = delim_pos + 1;
        bool      handled    = false;
        char      c          = m_Str[delim_pos];

        if ((m_Flags & NStr::fSplit_CanEscape) != 0  &&  c == '\\') {
            // treat the following character literally
            if (++pos > m_Str.size()) {
                //NCBI_THROW2(CStringException, eFormat, "Unescaped trailing \\", delim_pos);
                HBN_ERR("Unescaped trailing \\%zu", delim_pos);
            }
            handled = true;

        } else if ((m_Flags & NStr::fSplit_CanQuote) != 0) {
            if (active_quote != '\0') {
                if (c == active_quote) {
                    if (pos < m_Str.size()  &&  m_Str[pos] == active_quote) {
                        // count a doubled quote as one literal occurrence
                        ++pos;
                    } else {
                        active_quote = '\0';
                    }
                } else {
                    continue; // not actually a boundary
                }
                handled = true;
            } else if (((m_Flags & NStr::fSplit_CanSingleQuote) != 0  &&  c == '\'') || 
                       ((m_Flags & NStr::fSplit_CanDoubleQuote) != 0  &&  c == '"')) {
                active_quote = c;
                quote_pos    = delim_pos;
                handled = true;
            }
        }

        if ( !handled ) {
            if ((m_Flags & NStr::fSplit_ByPattern) != 0) {
                if (delim_pos + m_Delim.size() <= m_Str.size()
                    &&  (memcmp(m_Delim.data() + 1, m_Str.data() + pos,
                                m_Delim.size() - 1) == 0)) {
                    done = true;
                    next_start = pos = delim_pos + m_Delim.size();
                } else {
                    continue;
                }
            } else {
                done = true;
            }
            // save delimiter position
            if (ptr_delim_pos) {
                *ptr_delim_pos = delim_pos;
            }
        }

        if (delim_pos > part_start) {
            found_text = true;
            if (part_collector != NULL) {
                part_collector->Add(m_Str.substr(part_start, delim_pos - part_start));
            }
        }
        part_start = next_start;
    }

    if (active_quote != '\0') {
        //NCBI_THROW2(CStringException, eFormat, string("Unbalanced ") + active_quote, quote_pos);
        HBN_ERR("Unbalanced %c at position %zu", active_quote, quote_pos);
    }

    if (delim_pos == NPOS) {
        found_text = true;
        if (part_collector != NULL) {
            part_collector->Add(m_Str.substr(part_start));
        }
        m_Pos = NPOS;
    } else {
        m_Pos = pos;
        MergeDelims();
    }
    return found_text;
}


void CStrTokenizeBase::x_SkipDelims(bool force_skip)
{
    _ASSERT(NStr::fSplit_MergeDelimiters);

    if ( !force_skip  &&  (m_Flags & NStr::fSplit_MergeDelimiters) == 0 ) {
        return;
    }
    // skip all delimiters, starting from the current position
    if ((m_Flags & NStr::fSplit_ByPattern) == 0) {
        m_Pos = m_Str.find_first_not_of(m_Delim, m_Pos);
    } else {
        while (m_Pos != NPOS 
                &&  m_Pos + m_Delim.size() <= m_Str.size()
                && (memcmp(m_Delim.data(), m_Str.data() + m_Pos,
                            m_Delim.size()) == 0)) {
            m_Pos += m_Delim.size();
        }
    }
}


void CStrTokenizeBase::x_ExtendInternalDelim()
{
    if ( !(m_Flags & (NStr::fSplit_CanEscape | NStr::fSplit_CanQuote)) ) {
        return; // Nothing to do
    }
    SIZE_TYPE n = m_InternalDelim.size();
    char* buf = m_DelimStorage.Allocate(n + 3);
    char *s = buf;
    memcpy(s, m_InternalDelim.data(), n);
    if ((m_Flags & NStr::fSplit_CanEscape) != 0) {
        s[n++] = '\\';
    }
    if ((m_Flags & NStr::fSplit_CanSingleQuote) != 0) {
        s[n++] = '\'';
    }
    if ((m_Flags & NStr::fSplit_CanDoubleQuote) != 0) {
        s[n++] = '"';
    }
    m_InternalDelim.assign(buf, n);
}

END_NCBI_SCOPE