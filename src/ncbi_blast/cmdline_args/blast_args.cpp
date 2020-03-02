/* $Id: blast_args.cpp 581735 2019-03-05 16:42:54Z ivanov $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: blast_args.cpp

Author: Jason Papadopoulos

******************************************************************************/

/** @file blast_args.cpp
 * convert blast-related command line
 * arguments into blast options
*/

#include "blast_args.hpp"
#include "../../corelib/hbn_package_version.h"

#include <sstream>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

void
IBlastCmdLineArgs::ExtractAlgorithmOptions(const CArgs& /* cmd_line_args */,
                                           CBlastOptions& /* options */)
{}

CProgramDescriptionArgs::CProgramDescriptionArgs(const string& program_name,
                                                 const string& program_desc)
    : m_ProgName(program_name), m_ProgDesc(program_desc)
{}

void
CProgramDescriptionArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // program description
    arg_desc.SetUsageContext(m_ProgName, m_ProgDesc + " " + HBN_PACKAGE_VERSION);
}

CTaskCmdLineArgs::CTaskCmdLineArgs(const set<string>& supported_tasks,
                                   const string& default_task)
: m_SupportedTasks(supported_tasks), m_DefaultTask(default_task)
{
    _ASSERT( !m_SupportedTasks.empty() );
    if ( !m_DefaultTask.empty() ) {
        _ASSERT(m_SupportedTasks.find(m_DefaultTask) != m_SupportedTasks.end());
    }
}

void
CTaskCmdLineArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("General search options");
    if ( !m_DefaultTask.empty() ) {
        arg_desc.AddDefaultKey(kTask, "task_name", "Task to execute",
                               CArgDescriptions::eString, m_DefaultTask);
    } else {
        arg_desc.AddKey(kTask, "task_name", "Task to execute",
                        CArgDescriptions::eString);
    }
    arg_desc.SetConstraint(kTask, new CArgAllowStringSet(m_SupportedTasks));
    arg_desc.SetCurrentGroup("");

}

void
CTaskCmdLineArgs::ExtractAlgorithmOptions(const CArgs& /* cmd_line_args */,
                                          CBlastOptions& /* options */)
{
    // N.B.: handling of tasks occurs at the application level to ensure that
    // only relevant tasks are added (@sa CBlastnAppArgs)
}

void
CMTArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    // number of threads
    arg_desc.SetCurrentGroup("Miscellaneous options");
    const int kMinValue = static_cast<int>(kDfltNumThreads);

    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                           "Number of threads (CPUs) to use in the BLAST search",
                           CArgDescriptions::eInteger,
                           NStr::IntToString(kDfltNumThreads));
    arg_desc.SetConstraint(kArgNumThreads,
                           new CArgAllowValuesGreaterThanOrEqual(kMinValue));

    arg_desc.SetCurrentGroup("");
}

void
CMTArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& /* opts */)
{
    const int kMaxValue = static_cast<int>(hbn_get_cpu_count());

    if (args.Exist(kArgNumThreads) &&
        args[kArgNumThreads].HasValue()) {  // could be cancelled by the exclusion in CRemoteArgs

        // use the minimum of the two: user requested number of threads and
        // number of available CPUs for number of threads
        int num_threads = args[kArgNumThreads].AsInteger();
        if (num_threads > kMaxValue) {
            m_NumThreads = kMaxValue;

            string warn_msg = (string)"Number of threads was reduced to " +
                     NStr::IntToString((unsigned int)m_NumThreads) +
                     " to match the number of available CPUs";
            HBN_WARN("%s", warn_msg.c_str());
        }
        else {
            m_NumThreads = num_threads;
        }
    }
}

void CGridArgs::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    arg_desc.SetCurrentGroup("Miscellaneous options");
    arg_desc.AddDefaultKey(kArgGrid, "GRID_options",
                "Node configuration when multiple computing nodes are used,"
                " (Format: 'node_id num_nodes')",
                CArgDescriptions::eString,
                NStr::IntToString(kDfltNodeId) + " " + NStr::IntToString(kDfltNumNodes));
    arg_desc.SetCurrentGroup("");
}

void CGridArgs::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& /* opts */)
{
    if (args.Exist(kArgGrid) && args[kArgGrid].HasValue()) {
        string gridstr = args[kArgGrid].AsString();
        CTempString kDelim(" ");
        vector<string> components;
        NStr::Split(gridstr, " ", components);
        if (components.size() != 2) 
            HBN_ERR("Invalid value '%s' to argument '%s'", gridstr.c_str(), kArgGrid.c_str());
        m_NodeId = NStr::StringToInt(components[0]);
        m_NumNodes = NStr::StringToInt(components[1]);

        if (m_NodeId < 0) HBN_ERR("node id must be >=0: %s", gridstr.c_str());
        if (m_NumNodes <= 0) HBN_ERR("number of nodes must be >0: %s", gridstr.c_str());
        if (m_NodeId >= m_NumNodes) HBN_ERR("node index (%d) must be smaller than number of nodes (%d)", m_NodeId, m_NumNodes);
    }
}

string FindProgramDisplayName(const string& exec_path)
{
    hbn_assert(!exec_path.empty());
    string::size_type from = exec_path.size();
    while (from) {
        if (exec_path[from-1] == '/') break;
        --from;
    }
    string appname = exec_path.substr(from);
    return appname;
}

string PrintProgramVersion(const string& appname)
{
    ostringstream os;
    char build_info[512];
    hbn_build_info(build_info);

    os << appname << ": " << HBN_PACKAGE_VERSION << endl;
    os << "  Package: " << HBN_PACKAGE_NAME 
       << ' ' << HBN_PACKAGE_VERSION
       << ", build " << build_info << endl;;
    
    string v = os.str();
    return v;
}

END_SCOPE(blast)
END_NCBI_SCOPE