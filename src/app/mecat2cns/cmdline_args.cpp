#include "cmdline_args.h"

#include "../../corelib/seqdb.h"
#include "../../ncbi_blast/cmdline_args/blast_args.hpp"
#include "../../ncbi_blast/setup/blast_types.hpp"

#include <limits>
#include <sstream>

BEGIN_NCBI_SCOPE
USING_SCOPE(blast);
USING_SCOPE(align_format);

using namespace std;

const string kArgMemScKmerSize("memsc_kmer_size");
const int kDfltMemScKmerSize = 10;
const int kMemScKmerSizeMin = 8;
const int kMemScKmerSizeMax = 10;
const string kArgMemScKmerWindow("memsc_kmer_window");
const int kDfltMemScKmerWindow = 10;
const string kArgMemScMemSize("memsc_mem_size");
const int kDfltMemScMemSize = 15;
const string kArgMemScMemScore("memsc_mem_score");
const int kDfltMemScMemScore = 100;

const string kArgBatchSize("batch_size");
const int kDfltBatchSize = 5000;
const string kArgUseBatchMode("use_batch_mode");
const bool kDfltUseBatchMode = false;
const string kArgOvlpCovPerc("ovlp_cov_perc");
const double kDfltOvlpCovPerc = 60.0;
const string kArgOvlpCovRes("ovlp_cov_res");
const int kDfltOvlpCovRes = 1000;
const string kArgMinCov("min_cov");
const int kDfltMinCov = 4;
const string kArgMinSize("min_size");
const int kDfltMinSize = 2000;
const double kDfltPercIdentity = 70.0;

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

    arg_desc.AddDefaultKey(kArgBatchSize, "int_value",
                "Each time correct this number of raw reads",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltBatchSize));
    arg_desc.SetConstraint(kArgBatchSize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddFlag(kArgUseBatchMode, "Load raw reads into main memory in batches", true);

    arg_desc.AddDefaultKey(kArgOvlpCovPerc, "real_value",
                "Overlaps should cover at least this fraction of residues of the template or supporting read",
                CArgDescriptions::eDouble,
                NStr::DoubleToString(kDfltOvlpCovPerc));
    arg_desc.SetConstraint(kArgOvlpCovPerc, CArgAllowValuesBetween(0.0, 100.0));

    arg_desc.AddDefaultKey(kArgOvlpCovRes, "int_value",
                "Overlaps should cover at least this number of residues",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltOvlpCovRes));
    arg_desc.SetConstraint(kArgOvlpCovRes, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddDefaultKey(kArgPercentIdentity, "real_value",
                "Minimal identity percentage of overlaps",
                CArgDescriptions::eDouble,
                NStr::DoubleToString(kDfltPercIdentity));
    arg_desc.SetConstraint(kArgPercentIdentity, CArgAllowValuesBetween(0.0, 100.0));

    arg_desc.AddDefaultKey(kArgMinCov, "int_value",
                "Only correct subsequences that are covered by at least this value of supporting reads",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMinCov));
    arg_desc.SetConstraint(kArgMinCov, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMinSize, "int_value",
                "Minimal length of corrected sequences",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMinSize));
    arg_desc.SetConstraint(kArgMinSize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMemScKmerSize, "int_value",
                "Length of perfect matched kmers that are to be extended to MEMs",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScKmerSize));
    arg_desc.SetConstraint(kArgMemScKmerSize, CArgAllowValuesBetween(kMemScKmerSizeMin, kMemScKmerSizeMax));

    arg_desc.AddDefaultKey(kArgMemScKmerWindow, "int_value",
                "Kmer sampling window",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScKmerWindow));
    arg_desc.SetConstraint(kArgMemScKmerWindow, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMemScMemSize, "int_value",
                "Length of maximal exact match",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScMemSize));
    arg_desc.SetConstraint(kArgMemScMemSize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMemScMemScore, "int_value",
                "Minimum chaining score of a candidate",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMemScMemScore));
    arg_desc.SetConstraint(kArgMemScMemScore, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddDefaultKey(kArgNumThreads, "int_value",
                "Number of threads (CPUs) to use in the search",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltNumThreads));
    arg_desc.SetConstraint(kArgNumThreads, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgGrid, "GRID_options",
                "Index of this node when multiple computation nodes are used in the search\n"
                "(Format: 'node_id num_nodes')\n"
                "Default = '0 1'",
                CArgDescriptions::eString);
}

void CommandLineArguments::ExtractAlgorithmOptions(const CArgs& args, CBlastOptions& options)
{
    if (args.Exist(kArgBatchSize) && args[kArgBatchSize].HasValue()) {
        m_Options->batch_size = args[kArgBatchSize].AsInteger();
    }

    if (args.Exist(kArgUseBatchMode)) {
        m_Options->use_batch_mode = static_cast<bool>(args[kArgUseBatchMode]);
    }

    if (args.Exist(kArgOvlpCovPerc) && args[kArgOvlpCovPerc].HasValue()) {
        m_Options->ovlp_cov_perc = args[kArgOvlpCovPerc].AsDouble();
    }

    if (args.Exist(kArgOvlpCovRes) && args[kArgOvlpCovRes].HasValue()) {
        m_Options->ovlp_cov_res = args[kArgOvlpCovRes].AsInteger();
    }

    if (args.Exist(kArgPercentIdentity) && args[kArgPercentIdentity].HasValue()) {
        m_Options->perc_identity = args[kArgPercentIdentity].AsDouble();
    }

    if (args.Exist(kArgMinCov) && args[kArgMinCov].HasValue()) {
        m_Options->min_cov = args[kArgMinCov].AsInteger();
    }

    if (args.Exist(kArgMinSize) && args[kArgMinSize].HasValue()) {
        m_Options->min_size = args[kArgMinSize].AsInteger();
    }
 
    /// mem chaining scoring options
    if (args.Exist(kArgMemScKmerSize) && args[kArgMemScKmerSize].HasValue()) {
        m_Options->memsc_kmer_size = args[kArgMemScKmerSize].AsInteger();
    }

    if (args.Exist(kArgMemScKmerWindow) && args[kArgMemScKmerWindow].HasValue()) {
        m_Options->memsc_kmer_window = args[kArgMemScKmerWindow].AsInteger();
    }

    if (args.Exist(kArgMemScMemSize) && args[kArgMemScMemSize].HasValue()) {
        m_Options->memsc_mem_size = args[kArgMemScMemSize].AsInteger();
    }

    if (args.Exist(kArgMemScMemScore) && args[kArgMemScMemScore].HasValue()) {
        m_Options->memsc_score = args[kArgMemScMemScore].AsInteger();
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

    if (args.Exist(kArgGrid) && args[kArgGrid].HasValue()) {
        string gridstr = args[kArgGrid].AsString();
        CTempString kDelim(" ");
        vector<string> components;
        NStr::Split(gridstr, " ", components);
        if (components.size() != 2) 
            HBN_ERR("Invalid value '%s' to argument '%s'", gridstr.c_str(), kArgGrid.c_str());
        m_Options->node_id = NStr::StringToInt(components[0]);
        m_Options->num_nodes = NStr::StringToInt(components[1]);

        if (m_Options->node_id < 0) HBN_ERR("node id must be >=0: %s", gridstr.c_str());
        if (m_Options->num_nodes <= 0) HBN_ERR("number of nodes must be >0: %s", gridstr.c_str());
        if (m_Options->node_id >= m_Options->num_nodes) 
            HBN_ERR("node index (%d) must be smaller than number of nodes (%d)", 
                m_Options->node_id, m_Options->num_nodes);
    }
}

void Init_HbnProgramOptions(HbnProgramOptions* opts)
{
    opts->db_dir = NULL;
    opts->db_title = INIT_QUERY_DB_TITLE;
    opts->can_dir = NULL;
    opts->output = strdup(kDfltOutput.c_str());

    opts->batch_size = kDfltBatchSize;
    opts->use_batch_mode = kDfltUseBatchMode;
    opts->ovlp_cov_perc = kDfltOvlpCovPerc;
    opts->ovlp_cov_res = kDfltOvlpCovRes;
    opts->perc_identity = kDfltPercIdentity;
    opts->min_cov = kDfltMinCov;
    opts->min_size = kDfltMinSize;

    /// mem chaining scoring options
    opts->memsc_kmer_size = kDfltMemScKmerSize;
    opts->memsc_kmer_window = kDfltMemScKmerWindow;
    opts->memsc_mem_size = kDfltMemScMemSize;
    opts->memsc_score = kDfltMemScMemScore;

    /// misc options
    opts->num_threads = kDfltNumThreads;
    opts->node_id = kDfltNodeId;
    opts->num_nodes = kDfltNumNodes;
}

static void
s_PreCheckCmdLineArgs(int argc, char* argv[], CArgDescriptions* arg_desc)
{
    string kProgram = FindProgramDisplayName(argv[0]);
    string kHbnUsage = kProgram + " [OPTIONS] db_dir can_dir";

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
    string kProgramDescription("SMRT raw reads correction toolkit");
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
    if (argc - argv_idx  < 2) {
        HBN_ERR("The database directory and candidate dirctory must be specified");
    } else if (argc - argv_idx > 2) {
        string err = "Too many database directory and candidate directory values: '";
        for (int i = argv_idx; i < argc; ++i) {
            err += argv[i];
            if (i != argc - 1) err += ' ';
        }
        err += "'";
        HBN_ERR("%s", err.c_str());
    } else {
        opts->db_dir = argv[argv_idx];
        opts->can_dir = argv[argv_idx + 1];
    }
}

#define os_one_option_value(name, value) os << '-' << name << ' ' << value << ' '
#define os_one_flag_option(name) os << '-' << name << ' ';

extern "C"
char* HbnProgramOptions2String(const HbnProgramOptions* opts)
{
    ostringstream os;

    os_one_option_value(kArgBatchSize, opts->batch_size);
    if (opts->use_batch_mode) os_one_flag_option(kArgUseBatchMode);
    os_one_option_value(kArgOvlpCovPerc, opts->ovlp_cov_perc);
    os_one_option_value(kArgOvlpCovRes, opts->ovlp_cov_res);
    os_one_option_value(kArgPercentIdentity, opts->perc_identity);
    os_one_option_value(kArgMinCov, opts->min_cov);
    os_one_option_value(kArgMinSize, opts->min_size);

    /// mem chaining scoring options
    os_one_option_value(kArgMemScKmerSize, opts->memsc_kmer_size);
    os_one_option_value(kArgMemScKmerWindow, opts->memsc_kmer_window);
    os_one_option_value(kArgMemScMemSize, opts->memsc_mem_size);
    os_one_option_value(kArgMemScMemScore, opts->memsc_score);

    /// misc options
    os_one_option_value(kArgNumThreads, opts->num_threads);
    os << '-' << kArgGrid << ' ' << opts->node_id << ' ' << opts->num_nodes << ' ';

    os << endl
       << "database directory: " << opts->db_dir << endl
       << "candidate directory: " << opts->can_dir << endl;

    string os_str = os.str();
    return strdup(os_str.c_str());
}

END_NCBI_SCOPE