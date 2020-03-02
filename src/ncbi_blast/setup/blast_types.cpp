#include "blast_types.hpp"

#include "../str_util/ncbistr.hpp"

BEGIN_NCBI_SCOPE
BEGIN_SCOPE(blast)

using namespace std;

string EProgramToTaskName(EProgram p)
{
    string retval;
    switch (p) {
    case eBlastn:           retval.assign("blastn"); break;
    case eMegablast:        retval.assign("megablast"); break;
    case eDiscMegablast:    retval.assign("dc-megablast"); break;
    case eBlastp:           retval.assign("blastp"); break;
    case eBlastx:           retval.assign("blastx"); break;
    case eTblastn:          retval.assign("tblastn"); break;
    case eTblastx:          retval.assign("tblastx"); break;
    case eRPSBlast:         retval.assign("rpsblast"); break;
    case eRPSTblastn:       retval.assign("rpstblastn"); break;
    case ePSIBlast:         retval.assign("psiblast"); break;
    case ePSITblastn:       retval.assign("psitblastn"); break;
    case ePHIBlastp:        retval.assign("phiblastp"); break;
    case ePHIBlastn:        retval.assign("phiblastn"); break;
    case eDeltaBlast:       retval.assign("deltablast"); break;
    case eVecScreen:        retval.assign("vecscreen"); break;
    case eMapper:           retval.assign("mapr2g"); break;
    default:
        cerr << "Invalid EProgram value: " << (int)p << endl;
        abort();
    }

#if _DEBUG
    ThrowIfInvalidTask(retval);
#endif

    return retval;
}

EBlastProgramType
EProgramToEBlastProgramType(EProgram p)
{
    switch (p) {
    case eBlastn:
    case eMegablast:
    case eDiscMegablast:
    case eVecScreen:
        return eBlastTypeBlastn;

    case eMapper:
        return eBlastTypeMapping;
        
    case eBlastp:
        return eBlastTypeBlastp;
        
    case eBlastx:
        return eBlastTypeBlastx;
        
    case eTblastn:
        return eBlastTypeTblastn;
        
    case eTblastx:
        return eBlastTypeTblastx;
        
    case eRPSBlast:
        return eBlastTypeRpsBlast;
        
    case eRPSTblastn:
        return eBlastTypeRpsTblastn;
        
    case ePSIBlast:
    case eDeltaBlast:
        return eBlastTypePsiBlast;
        
    case ePSITblastn:
        return eBlastTypePsiTblastn;

    case ePHIBlastp:
        return eBlastTypePhiBlastp;
        
    case ePHIBlastn:
        return eBlastTypePhiBlastn;
        
    default:
        return eBlastTypeUndefined;
    }
}

EProgram ProgramNameToEnum(const std::string& program_name)
{
    _ASSERT( !program_name.empty() );

    string lowercase_program_name(program_name);
    lowercase_program_name = NStr::ToLower(lowercase_program_name);

#if _DEBUG
    ThrowIfInvalidTask(lowercase_program_name);
#endif

    if (NStr::StartsWith(lowercase_program_name, "blastn")) {
        return eBlastn;
    // -RMH- support new toolkit program
    } else if (NStr::StartsWith(lowercase_program_name, "rmblastn")) {
        return eBlastn;
    } else if (NStr::StartsWith(lowercase_program_name, "blastp")) {
        return eBlastp;
    } else if (lowercase_program_name == "blastx") {
        return eBlastx;
    } else if (lowercase_program_name == "tblastn") {
        return eTblastn;
    } else if (lowercase_program_name == "tblastx") {
        return eTblastx;
    } else if (lowercase_program_name == "rpsblast") {
        return eRPSBlast;
    } else if (lowercase_program_name == "rpstblastn") {
        return eRPSTblastn;
    } else if (lowercase_program_name == "megablast") {
        return eMegablast; 
    } else if (lowercase_program_name == "psiblast") {
        return ePSIBlast;
    } else if (lowercase_program_name == "psitblastn") {
        return ePSITblastn;
    } else if (lowercase_program_name == "dc-megablast") {
        return eDiscMegablast;
    } else if (lowercase_program_name == "deltablast") {
        return eDeltaBlast;
    } else if (lowercase_program_name == "vecscreen") {
        return eVecScreen;
    // FIXME: mapper is used in core as a single program name for all tasks,
    // we may need a better approach to mapping tasks with fewer program names
    } else if (lowercase_program_name == "mapper") {
        return eMapper;
    } else if (lowercase_program_name == "mapr2g") {
        return eMapper;
    } else if (lowercase_program_name == "mapr2r") {
        return eMapper;
    } else if (lowercase_program_name == "mapg2g") {
        return eMapper;
    } else {
        HBN_ERR("Program type '%s' not supported", program_name.c_str());
    }
}

set<string>
GetTasks(ETaskSets choice /* = eAll */)
{
    set<string> retval;
    if (choice == eNuclNucl || choice == eAll) {
        retval.insert("blastn");
        retval.insert("blastn-short");
        retval.insert("megablast");
        retval.insert("dc-megablast");
        retval.insert("vecscreen");
        // -RMH-
        retval.insert("rmblastn");
    }

    if (choice == eProtProt || choice == eAll) {
        retval.insert("blastp");
        retval.insert("blastp-short");
        retval.insert("blastp-fast");
        // retval.insert("kblastp");
    }

    if (choice == eAll) {
        retval.insert("psiblast");
        //retval.insert("phiblastn"); // not supported yet
        retval.insert("phiblastp");
        retval.insert("rpsblast");
        retval.insert("rpstblastn");
        retval.insert("blastx");
        retval.insert("blastx-fast");
        retval.insert("deltablast");
        retval.insert("tblastn");
        retval.insert("tblastn-fast");
        retval.insert("psitblastn");
        retval.insert("tblastx");
        retval.insert("kblastp");
    }

    if (choice == eMapping || choice == eAll) {
        retval.insert("mapper");
        retval.insert("mapr2g");
        retval.insert("mapr2r");
        retval.insert("mapg2g");
    }

    return retval;
}

void ThrowIfInvalidTask(const string& task)
{
#if 0
    set<string> valid_tasks;
    if (valid_tasks.empty()) {
        valid_tasks = CBlastOptionsFactory::GetTasks();
    }

    if (valid_tasks.find(task) == valid_tasks.end()) {
        ostringstream os;
        os << "'" << task << "' is not a supported task";
        NCBI_THROW(CBlastException, eInvalidArgument, os.str());
    }
#endif
    set<string> valid_tasks = GetTasks();
    if (valid_tasks.find(task) == valid_tasks.end()) {
        HBN_ERR("'%s' is not a supported task", task.c_str());
    }
}

END_SCOPE(blast)
END_NCBI_SCOPE