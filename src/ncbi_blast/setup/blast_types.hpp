/*  $Id: blast_types.hpp 573790 2018-11-01 15:46:35Z ivanov $
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
 * Author:  Ilya Dondoshansky
 *
 */

/** @file blast_types.hpp
 * Definitions of special type used in BLAST
 */

#ifndef ALGO_BLAST_API___BLAST_TYPE__HPP
#define ALGO_BLAST_API___BLAST_TYPE__HPP

#if 0
#include <corelib/ncbistd.hpp>
#include <objects/seqalign/Seq_align_set.hpp>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_message.h>a
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_filter.h>
#endif

#include "../ncbi_blast_aux.hpp"
#include "blast_message.h"
#include "blast_program.h"
#include "blast_types.h"

#include <set>
#include <vector>

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

using std::set;
using std::vector;

#if 0
/// This enumeration is to evolve into a task/program specific list that 
/// specifies sets of default parameters to easily conduct searches using
/// BLAST.
/// @todo EProgram needs to be renamed to denote a task (similar to those
/// exposed by the BLAST web page) rather than a program type
/// N.B.: When making changes to this enumeration, please update 
/// blast::ProgramNameToEnum (blast_aux.[ch]pp), blast::GetNumberOfFrames
/// (blast_setup_cxx.cpp) and BlastNumber2Program and BlastProgram2Number
/// (blast_util.c)
enum EProgram {
    eBlastNotSet = 0,   ///< Not yet set.
    eBlastn,            ///< Nucl-Nucl (traditional blastn)
    eBlastp,            ///< Protein-Protein
    eBlastx,            ///< Translated nucl-Protein
    eTblastn,           ///< Protein-Translated nucl
    eTblastx,           ///< Translated nucl-Translated nucl
    eRPSBlast,          ///< protein-pssm (reverse-position-specific BLAST)
    eRPSTblastn,        ///< nucleotide-pssm (RPS blast with translated query)
    eMegablast,         ///< Nucl-Nucl (traditional megablast)
    eDiscMegablast,     ///< Nucl-Nucl using discontiguous megablast
    ePSIBlast,          ///< PSI Blast
    ePSITblastn,        ///< PSI Tblastn
    ePHIBlastp,         ///< Protein PHI BLAST
    ePHIBlastn,         ///< Nucleotide PHI BLAST
    eDeltaBlast,        ///< Delta Blast
    eVecScreen,         ///< Vector screening
    eMapper,            ///< Jumper alignment for mapping
    eKBlastp,            ///< KMER screening and BLASTP
    eBlastProgramMax    ///< Undefined program
};
#endif

/** Convert a EProgram enumeration value to a task name (as those used in the
 * BLAST command line binaries)
 * @param p EProgram enumeration value to convert [in]
 */
NCBI_XBLAST_EXPORT 
string EProgramToTaskName(EProgram p);

/// Map a string into an element of the ncbi::blast::EProgram enumeration 
/// (except eBlastProgramMax).
/// @param program_name [in]
/// @return an element of the ncbi::blast::EProgram enumeration, except
/// eBlastProgramMax
/// @throws CBlastException if the string does not map into any of the EProgram
/// elements
NCBI_XBLAST_EXPORT
EProgram ProgramNameToEnum(const std::string& program_name);

/// Validates that the task provided is indeed a valid task, otherwise throws a
/// CBlastException
/// @param task task name to validate [in]
NCBI_XBLAST_EXPORT
void ThrowIfInvalidTask(const string& task);

/// Convert EProgram to EBlastProgramType.
/// @param p Program expressed as an api layer EProgram.
/// @return Same program using the core enumeration.
NCBI_XBLAST_EXPORT
EBlastProgramType
EProgramToEBlastProgramType(EProgram p);

/// Sets of tasks for the command line BLAST binaries
enum ETaskSets {
    eNuclNucl,      ///< Nucleotide-nucleotide tasks
    eProtProt,      ///< Protein-protein tasks
    eMapping,       ///< Mapping tasks
    eAll            ///< Retrieve all available tasks
};

/// Retrieve the set of supported tasks
set<string> GetTasks(ETaskSets choice = eAll);

END_SCOPE(blast)
END_NCBI_SCOPE

#endif  /* ALGO_BLAST_API___BLAST_TYPE__HPP */
