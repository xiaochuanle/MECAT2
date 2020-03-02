/*  $Id: blast_args.hpp 579216 2019-01-31 16:18:17Z ivanov $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
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
 * ===========================================================================
 *
 * Author:  Jason Papadopoulos
 *
 */

/** @file blast_args.hpp
 * Interface for converting blast-related command line
 * arguments into blast options
 */

#ifndef ALGO_BLAST_BLASTINPUT___BLAST_ARGS__HPP
#define ALGO_BLAST_BLASTINPUT___BLAST_ARGS__HPP

#if 0
#include <corelib/ncbistd.hpp>
#include <corelib/ncbiargs.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_options.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/igblast/igblast.hpp>
#include <algo/blast/api/setup_factory.hpp> // for CThreadable
#include <algo/blast/blastinput/cmdline_flags.hpp>
#include <algo/blast/blastinput/blast_input_aux.hpp>

#include <objmgr/scope.hpp>     // for CScope
#include <objects/seqloc/Na_strand.hpp>
#include <objects/scoremat/PssmWithParameters.hpp>

#include <util/compress/stream_util.hpp>
#endif

#include "ncbiargs_desc.hpp"
#include "cmdline_flags.hpp"

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

struct CBlastOptions {};

/**
 * BLAST Command line arguments design
 * The idea is to have several small objects (subclasses of IBlastCmdLineArgs)
 * which can do two things:
 * 1) On creation, add flags/options/etc to a CArgs object
 * 2) When passed in a CBlastOptions object, call the appropriate methods based
 * on the CArgs options set when the NCBI application framework parsed the
 * command line. If data collected by the small object (from the command line)
 * cannot be applied to the CBlastOptions object, then it's provided to the
 * application via some other interface methods.
 *
 * Each command line application will have its own argument class (e.g.:
 * CPsiBlastAppArgs), which will contain several of the aformentioned small
 * objects. It will create and hold a reference to a CArgs class as well as
 * a CBlastOptionsHandle object, which will pass to each of its small objects
 * aggregated as data members and then return it to the caller (application)
 *
 * Categories of data to extract from command line options
 * 1) BLAST algorithm options
 * 2) Input/Output files, and their modifiers (e.g.: believe query defline)
 * 3) BLAST database information (names, limitations, num db seqs)
 * 4) Formatting options (html, display formats, etc)
*/

/** Interface definition for a generic command line option for BLAST
 */
class NCBI_BLASTINPUT_EXPORT IBlastCmdLineArgs : public CObject
{
public:
    /** Our virtual destructor */
    virtual ~IBlastCmdLineArgs() {}

    /** Sets the command line descriptions in the CArgDescriptions object
     * relevant to the subclass
     * @param arg_desc the argument descriptions object [in|out]
     */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc) = 0;

    /** Extracts BLAST algorithmic options from the command line arguments into
     * the CBlastOptions object. Default implementation does nothing.
     * @param cmd_line_args Command line arguments parsed by the NCBI
     * application framework [in]
     * @param options object to which the appropriate options will be set
     * [in|out]
     */
    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);
};

/** Argument class to populate an application's name and description */
class NCBI_BLASTINPUT_EXPORT CProgramDescriptionArgs : public IBlastCmdLineArgs
{
public:
    /**
     * @brief Constructor
     *
     * @param program_name application's name [in]
     * @param program_description application's description [in]
     */
    CProgramDescriptionArgs(const string& program_name,
                            const string& program_description);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);

protected:
    string m_ProgName;  ///< Application's name
    string m_ProgDesc;  ///< Application's description
};

/// Argument class to specify the supported tasks a given program
class NCBI_BLASTINPUT_EXPORT CTaskCmdLineArgs : public IBlastCmdLineArgs
{
public:
    /** Constructor
     * @param supported_tasks list of supported tasks [in]
     * @param default_task One of the tasks above, to be displayed as
     * default in the command line arguments (cannot be empty or absent from
     * the set above) [in]
     */
    CTaskCmdLineArgs(const set<string>& supported_tasks,
                     const string& default_task);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);
private:
    /// Set of supported tasks by this command line argument
    const set<string> m_SupportedTasks;
    /// Default task for this command line argument
    string m_DefaultTask;
};

/** Argument class for collecting filtering options */
class NCBI_BLASTINPUT_EXPORT CFilteringArgs : public IBlastCmdLineArgs
{
public:
    /**
     * @brief Constructor
     *
     * @param query_is_protein is the query sequence(s) protein? [in]
     * @param filter_by_default should filtering be applied by default? [in]
     */
    CFilteringArgs(bool query_is_protein = true,
                   bool filter_by_default = true)
        : m_QueryIsProtein(query_is_protein),
          m_FilterByDefault(filter_by_default) {}

    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);
private:
    bool m_QueryIsProtein;  /**< true if the query is protein */
    bool m_FilterByDefault; /**< Should filtering be applied by default? */

    /**
     * @brief Auxiliary method to tokenize the filtering string.
     *
     * @param filtering_args string to tokenize [in]
     * @param output vector with tokens [in|out]
     */
    void x_TokenizeFilteringArgs(const string& filtering_args,
                                 vector<string>& output) const;
};

/// Defines values for match and mismatch in nucleotide comparisons as well as
/// non-greedy extension
class NCBI_BLASTINPUT_EXPORT CNuclArgs : public IBlastCmdLineArgs
{
public:
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);
};

/// Argument class to collect multi-threaded arguments
class NCBI_BLASTINPUT_EXPORT CMTArgs : public IBlastCmdLineArgs
{
public:
    /// Default Constructor
    CMTArgs(size_t default_num_threads = kDfltNumThreads) :
    	m_NumThreads(default_num_threads)
    {
    }
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);
    /** Interface method, \sa IBlastCmdLineArgs::SetArgumentDescriptions */
    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);

    /// Get the number of threads to spawn
    size_t GetNumThreads() const { return m_NumThreads; }
protected:
    size_t m_NumThreads;        ///< Number of threads to spawn
};

class CGridArgs : public IBlastCmdLineArgs
{
public:
    CGridArgs(int node_id = kDfltNodeId, int num_nodes = kDfltNumNodes)
        : m_NodeId(node_id), m_NumNodes(num_nodes) {}


    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);

    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args,
                                         CBlastOptions& options);

    int GetNodeId() const { return m_NodeId; }
    int GetNumNodes() const { return m_NumNodes; }
    
protected:
    int m_NodeId;
    int m_NumNodes;
};

/// Type definition of a container of IBlastCmdLineArgs
typedef vector< CRef<IBlastCmdLineArgs> > TBlastCmdLineArgs;



/**
 * @brief Create a CArgDescriptions object and invoke SetArgumentDescriptions
 * for each of the TBlastCmdLineArgs in its argument list
 *
 * @param args arguments to configure the return value [in]
 *
 * @return a CArgDescriptions object with the command line options set
 */
NCBI_BLASTINPUT_EXPORT
CArgDescriptions*
SetUpCommandLineArguments(TBlastCmdLineArgs& args);

string PrintProgramVersion(const string& appname);

string FindProgramDisplayName(const string& exec_path);

END_SCOPE(blast)
END_NCBI_SCOPE

#endif  /* ALGO_BLAST_BLASTINPUT___BLAST_ARGS__HPP */
