#include "cmdline_args.h"

#include "../../../corelib/seqdb.h"
#include "../../../ncbi_blast/cmdline_args/blast_args.hpp"
#include "../../../ncbi_blast/setup/blast_types.hpp"

#include <limits>
#include <sstream>

BEGIN_NCBI_SCOPE
USING_SCOPE(blast);
USING_SCOPE(align_format);

using namespace std;

const string kArgMinSize("min_size");
const int kDfltMinSize = 1000;

const string kDfltOutput("-");

static const char* kGroupGeneralSearchOptions = "General search options";

class CommandLineArguments : public IBlastCmdLineArgs
{
public:
    CommandLineArguments(HbnProgramOptions* options);
    ~CommandLineArguments();

    virtual void SetArgumentDescriptions(CArgDescriptions& arg_desc);

    virtual void ExtractAlgorithmOptions(const CArgs& cmd_line_args, CBlastOptions& options);

private:
    HbnProgramOptions*              m_Options;
};

CommandLineArguments::CommandLineArguments(HbnProgramOptions* options)
{
    m_Options = options;
}

CommandLineArguments::~CommandLineArguments()
{

}

void CommandLineArguments::SetArgumentDescriptions(CArgDescriptions& arg_desc)
{
    /// create the groups so that the ordering is established
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);

    /// General search options
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);

    arg_desc.AddDefaultKey(kArgOutput, "File_Out",
                "Output file name",
                CArgDescriptions::eString,
                kDfltOutput);

    arg_desc.AddDefaultKey(kArgMinSize, "int_value",
                "Minimum length of trimmed sequences",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMinSize));
    arg_desc.SetConstraint(kArgMinSize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                "Number of threads (CPUs) to use",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltNumThreads));
    arg_desc.SetConstraint(kArgNumThreads, CArgAllowValuesGreaterThanOrEqual(1));
}

void CommandLineArguments::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& options)
{
    if (args.Exist(kArgMinSize) && args[kArgMinSize].HasValue()) {
        m_Options->min_size = args[kArgMinSize].AsInteger();
    }

    if (args.Exist(kArgOutput) && args[kArgOutput].HasValue()) {
        m_Options->output = strdup(args[kArgOutput].AsString().c_str());
    }

    /// misc options
    const int kMaxValue = static_cast<int>(hbn_get_cpu_count());
    if (args.Exist(kArgNumThreads) &&
        args[kArgNumThreads].HasValue()) {  // could be cancelled by the exclusion in CRemoteArgs

        // use the minimum of the two: user requested number of threads and
        // number of available CPUs for number of threads
        int num_threads = args[kArgNumThreads].AsInteger();
        if (num_threads > kMaxValue) {
            m_Options->num_threads = kMaxValue;

            string warn_msg = (string)"Number of threads was reduced to " +
                     NStr::IntToString((unsigned int)m_Options->num_threads) +
                     " to match the number of available CPUs";
            HBN_WARN("%s", warn_msg.c_str());
        }
        else {
            m_Options->num_threads = num_threads;
        }
    }
}

void Init_HbnProgramOptions(HbnProgramOptions* opts)
{
    opts->db_dir = NULL;
    opts->pm4_dir = NULL;
    opts->output = strdup(kDfltOutput.c_str());
    opts->min_size = kDfltMinSize;
    opts->num_threads = kDfltNumThreads;
}

static void
s_PreCheckCmdLineArgs(int argc, char* argv[], CArgDescriptions* arg_desc)
{
    string kProgram = FindProgramDisplayName(argv[0]);
    string kHbnUsage = kProgram + " [OPTIONS] db_dir pm4_dir lcr_path";

    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') continue;
        if (NStr::CompareCase(argv[i] + 1, kArgHelp) == 0) {
            string usage_info;
            //arg_desc->PrintUsage(usage_info);
            arg_desc->HbnPrintUsage(kHbnUsage, usage_info);
            cout << usage_info << endl;
            exit(0);
        } else if (NStr::CompareCase(argv[i] + 1, kArgFullHelp) == 0) {
            string usage_info;
            //arg_desc->PrintUsage(usage_info, true);
            arg_desc->HbnPrintUsage(kHbnUsage, usage_info, true);
            cout << usage_info << endl;
            exit(0);
        } else if (NStr::CompareCase(argv[i] + 1, kArgVersion) == 0) {
            string appname = FindProgramDisplayName(argv[0]);
            string version = PrintProgramVersion(appname);
            cout << version << endl;
            exit(0);
        }
    }
}

extern "C"
void ParseHbnProgramCmdLineArguments(int argc, char* argv[], HbnProgramOptions* opts)
{
    Init_HbnProgramOptions(opts);
    TBlastCmdLineArgs arg_list;

    /// setup description
    CRef<IBlastCmdLineArgs> arg;
    string kProgram = FindProgramDisplayName(argv[0]);
    string kProgramDescription("Find largest clear ranges of sequences using their overlapping results");
    arg.reset(new CProgramDescriptionArgs(kProgram, kProgramDescription));
    arg_list.push_back(arg);

    arg.reset(new CommandLineArguments(opts));
    arg_list.push_back(arg);

    unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg_iter, arg_list) {
        (*arg_iter)->SetArgumentDescriptions(*arg_desc);
    }

    /// examine trivial arguments (-help, -h, -version)
    s_PreCheckCmdLineArgs(argc, argv, arg_desc.get());

    /// process nontrivial arguments
    unique_ptr<CArgs> cmd_args(new CArgs);
    int argv_idx = 1;
    auto& supported_args = (*arg_desc).GetArgs();
    while (argv_idx < argc) {
        if (argv[argv_idx][0] != '-') break;
        string argname = argv[argv_idx] + 1;
        //cout << "process " << argname << endl;
        //if (!(*arg_desc).Exist(argname)) HBN_ERR("unrecognised argument '%s'", argname.c_str());
        bool negative = false;
        auto it = (*arg_desc).Find(argname, &negative);
        if (it == supported_args.end()) HBN_ERR("unrecognised argument '%s'", argname.c_str());
        CArgDesc& arg = **it;
        CArgValue* av = nullptr;
        if (ArgDescIsFlag(arg)) {
            //av = arg.ProcessDefault();
            av = arg.ProcessArgument(kEmptyStr);
            ++argv_idx;
        } else {
            if (argv_idx + 1 >= argc) HBN_ERR("Mandatory value to argument '%s' is missing", argname.c_str());
            av = arg.ProcessArgument(argv[argv_idx + 1]);
            argv_idx += 2;
        }
        cmd_args->Add(av, true, true);
    }

    CBlastOptions cblastopts;
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg_iter, arg_list) {
        (*arg_iter)->ExtractAlgorithmOptions(*cmd_args, cblastopts);
    }

    /// query, subject
    if (argc - argv_idx < 3) {
        HBN_ERR("The database directory, m4 dirctory and lcr path must be specified");
    } else if (argc - argv_idx > 3) {
        string err = "Too many database directory, m4 directory and lcr path values: '";
        for (int i = argv_idx; i < argc; ++i) {
            err += argv[i];
            if (i != argc - 1) err += ' ';
        }
        err += "'";
        HBN_ERR("%s", err.c_str());
    } else {
        opts->db_dir = argv[argv_idx];
        opts->pm4_dir = argv[argv_idx+1];
        opts->lcr_path = argv[argv_idx+2];
    }
}

#define os_one_option_value(name, value) os << '-' << name << ' ' << value << ' '
#define os_one_flag_option(name) os << '-' << name << ' ';

extern "C"
char* HbnProgramOptions2String(const HbnProgramOptions* opts)
{
    ostringstream os;
    os_one_option_value(kArgMinSize, opts->min_size);
    os_one_option_value(kArgNumThreads, opts->num_threads);

    os << endl
       << "database directory: " << opts->db_dir << endl
       << "m4 directory: " << opts->pm4_dir << endl;

    string os_str = os.str();
    return strdup(os_str.c_str());
}

END_NCBI_SCOPE