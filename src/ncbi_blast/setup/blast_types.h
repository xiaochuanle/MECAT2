#ifndef __BLAST_TYPES_H
#define __BLAST_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

/// This enumeration is to evolve into a task/program specific list that 
/// specifies sets of default parameters to easily conduct searches using
/// BLAST.
/// @todo EProgram needs to be renamed to denote a task (similar to those
/// exposed by the BLAST web page) rather than a program type
/// N.B.: When making changes to this enumeration, please update 
/// blast::ProgramNameToEnum (blast_aux.[ch]pp), blast::GetNumberOfFrames
/// (blast_setup_cxx.cpp) and BlastNumber2Program and BlastProgram2Number
/// (blast_util.c)
typedef enum EProgram {
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
} EProgram;

#ifdef __cplusplus
}
#endif

#endif // __BLAST_TYPES_H