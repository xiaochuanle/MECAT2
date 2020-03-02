#include "cmdline_args.h"

#include "../../ncbi_blast/cmdline_args/blast_args.hpp"
#include "../../ncbi_blast/setup/blast_types.hpp"

#include <limits>
#include <sstream>

BEGIN_NCBI_SCOPE
USING_SCOPE(blast);
USING_SCOPE(align_format);

using namespace std;

const string kArgMinQuerySize("min_query_size");
const int kDfltMinQuerySize = 0;
const string kArgMaxQueryVolSeqs("max_query_vol_seqs");
const int kDfltMaxQueryVolSeqs = numeric_limits<int>::max();
const string kArgMaxQueryVolRes("max_query_vol_res");
const size_t kDfltMaxQueryVolRes = static_cast<size_t>(4000000000);
const string kArgMinSubjectSize("min_subject_size");
const int kDfltMinSubjectSize = 0;
const string kArgMaxSubjectVolSeqs("max_subject_vol_seqs");
const int kDfltMaxSubjectVolSeqs = numeric_limits<int>::max();
const string kArgMaxSubjectVolRes("max_subject_vol_res");
const size_t kDfltMaxSubjectVolRes = static_cast<size_t>(4000000000);
const string kArgDbDir("db_dir");
const string kDfltDbDir("hbndb");
const string kArgKeepDb("keep_db");
const bool kDfltKeepDb = false;

const string kArgKmerSize("kmer_size");
const int kDfltKmerSize = 15;
const string kArgKmerWindowSize("kmer_window");
const int kDfltKmerWindowSize = 10;
const string kArgMaxKmerOcc("max_kmer_occ");
const int kDfltMaxKmerOcc = 200;
const string kArgBlockSize("block_size");
const int kDfltBlockSize = 2000;
const string kArgMinDDFS("min_ddfs");
const int kDfltMinDDFS = 3;

const string kArgMemScKmerSize("memsc_kmer_size");
const int kDfltMemScKmerSize = 10;
const int kMemScKmerSizeMin = 8;
const int kMemScKmerSizeMax = 12;
const string kArgMemScMemSize("memsc_mem_size");
const int kDfltMemScMemSize = 15;
const string kArgMemScKmerWindow("memsc_kmer_window");
const int kDfltMemScKmerWindow = 10;
const string kArgMemScMemScore("memsc_mem_score");
const int kDfltMemScMemScore = 100;

const string kArgQueryCovHspRes("qcov_hsp_res");
const int kDfltQueryCovHspRes = 100;
const string kDfltHbnOutput("-");
const string kArgSkipOverhang("skip_overhang");
const bool kDfltSkipOverhang = false;

const string kArgOutputCigar("cigar");
const bool kDfltOutputCigar = false;
const string kArgOutputMd("md");
const bool kDfltOutputMd = false;

static const char* kGroupGeneralSearchOptions = "General search options";
static const char* kGroupInputQuery = "Input query options";
static const char* kGroupDbOptions = "Database options";
static const char* kGroupDDFSc = "DDF scoring options";
static const char* kGroupMemSc = "MEM chaining scoring options";
static const char* kGroupFormat = "Formatting options";
static const char* kGroupQueryFiltering = "Query Filtering Options";
static const char* kGroupRestrictSearch = "Restrict search or results";
static const char* kGroupStatistics = "Statistical options";
static const char* kGroupExtension = "Extension options";
static const char* kGroupMiscellaneous = "Miscellaneous options";

#define HBN_REWARD  (1)
#define HBN_PENELTY (-2)
#define HBN_GAP_OPEN    (0)
#define HBN_GAP_EXTEND  (2)

const int kDfltMaxHspsPerSubject = 5;
const EOutputFormat kDfltHbnOutfmt = eM4;
const EHbnTask kDfltTask = eHbnTask_rm;

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
    arg_desc.SetCurrentGroup(kGroupInputQuery);
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);
    arg_desc.SetCurrentGroup(kGroupDbOptions);
    arg_desc.SetCurrentGroup(kGroupDDFSc);
    arg_desc.SetCurrentGroup(kGroupMemSc);
    arg_desc.SetCurrentGroup(kGroupFormat);
    arg_desc.SetCurrentGroup(kGroupQueryFiltering);
    arg_desc.SetCurrentGroup(kGroupRestrictSearch);
    arg_desc.SetCurrentGroup(kGroupStatistics);
    arg_desc.SetCurrentGroup(kGroupMiscellaneous);

    /// General search options
    arg_desc.SetCurrentGroup(kGroupGeneralSearchOptions);

    string kDefaultTask = hbn_task_names[kDfltTask];
    arg_desc.AddDefaultKey(kTask, "task_name",
                "Task to execute",
                CArgDescriptions::eString,
                kDefaultTask);
    set<string> supported_tasks;
    supported_tasks.insert(string(hbn_task_names[eHbnTask_pm]));
    supported_tasks.insert(string(hbn_task_names[eHbnTask_rm]));
    arg_desc.SetConstraint(kTask, new CArgAllowStringSet(supported_tasks));

    arg_desc.AddDefaultKey(kArgOutput, "File_Out",
                "Output file name",
                CArgDescriptions::eString,
                kDfltHbnOutput);
    
    arg_desc.AddDefaultKey(kArgEvalue, "evalue",
                "Expectation value (E) threshold for saving hits",
                CArgDescriptions::eDouble,
                NStr::DoubleToString(BLAST_EXPECT_VALUE));

    arg_desc.AddOptionalKey(kArgMismatch, "penalty",
                "Penalty for a nucleotide mismatch",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMismatch, new CArgAllowValuesLessThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMatch, "reward",
                "Reward for a nucleotide match",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMatch, new CArgAllowValuesGreaterThanOrEqual(0));
    
    arg_desc.AddOptionalKey(kArgGapOpen, "open_penalty",
                "Cost to open a gap",
                CArgDescriptions::eInteger);
    
    arg_desc.AddOptionalKey(kArgGapExtend, "extend_penalty",
                "Cost to extend a gap",
                CArgDescriptions::eInteger);
                
    /// input query options
    arg_desc.SetCurrentGroup(kGroupInputQuery);

    arg_desc.AddDefaultKey(kArgStrand, "strand",
                "Query strand(s) to search against database/subject",
                CArgDescriptions::eString,
                kDfltArgStrand);
    arg_desc.SetConstraint(kArgStrand, &(* new CArgAllow_Strings, kDfltArgStrand, "plus", "minus"));

    /// database options
    arg_desc.SetCurrentGroup(kGroupDbOptions);

    arg_desc.AddOptionalKey(kArgDbDir, "directory",
                "Directory to store the query and subject database",
                CArgDescriptions::eString);

    arg_desc.AddFlag(kArgKeepDb, "Do not delete the database after search?", true);

    arg_desc.AddOptionalKey(kArgMinQuerySize, "int_value",
                "Skip query sequences shorter than this value",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinQuerySize, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxQueryVolSeqs, "int_value",
                "Maximum number of sequences in one volume of query database",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinQuerySize, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxQueryVolRes, "volume_size",
                "Maximum number of residues in one volume of query database",
                CArgDescriptions::eDataSize,
                NStr::UInt8ToString_DataSize(kDfltMaxQueryVolRes));
    arg_desc.SetConstraint(kArgMaxQueryVolRes, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgMinSubjectSize, "int_value",
                "Skip subject sequences shorter than this value",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinSubjectSize, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxSubjectVolSeqs, "int_value",
                "Maximum number of sequences in one volume of subject database",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxSubjectVolSeqs, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxSubjectVolRes, "volume_size",
                "Maximum number of residues in one volume of subject database",
                CArgDescriptions::eDataSize,
                NStr::UInt8ToString_DataSize(kDfltMaxSubjectVolRes));
    arg_desc.SetConstraint(kArgMaxSubjectVolRes, CArgAllowValuesGreaterThanOrEqual(1));

    /// DDF scoring options
    arg_desc.SetCurrentGroup(kGroupDDFSc);

    arg_desc.AddOptionalKey(kArgKmerSize, "int_value",
                "Kmer size for DDF scoring algorithm (length of best perfect match)",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgKmerSize, CArgAllowValuesBetween(10, 32));

    arg_desc.AddOptionalKey(kArgKmerWindowSize, "int_value",
                "Kmer sampling window size in subject sequences",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgKmerWindowSize, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgBlockSize, "int_value",
                "Split subject database into consecutive blocks, each having this number of residues",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgBlockSize, CArgAllowValuesGreaterThanOrEqual(500));

    arg_desc.AddOptionalKey(kArgMinDDFS, "int_value",
                "Minimum DDF score of a candidate",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMinDDFS, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddDefaultKey(kArgMaxKmerOcc, "int_value",
                "Filter out kmers occur larger than this value",
                CArgDescriptions::eInteger,
                NStr::IntToString(kDfltMaxKmerOcc));
    arg_desc.SetConstraint(kArgMaxKmerOcc, CArgAllowValuesGreaterThanOrEqual(1));

    /// output format
    arg_desc.SetCurrentGroup(kGroupFormat);

    string OutputFormatDescription = string(
        "alignment view options:\n"
        "  seqid  = candidate aligned sequence index,\n"
        "  seqidx = candidate aligned sequence index in binary format,\n"
        "  subseq = candidate aligned subsequence,\n"
        "  m4     = m4 alignment format,\n"  
        "  m4x    = binary m4 alignment format,\n"
        "  paf    = PAF format,\n"
        "  sam    = Sequence Alignment/Map (SAM)"
    );
    arg_desc.AddDefaultKey(kArgOutputFormat, "format",
                OutputFormatDescription,
                CArgDescriptions::eString,
                NStr::IntToString(kDfltArgOutputFormat));

    arg_desc.AddFlag(kArgOutputCigar, 
                "output CIGAR in m4 or paf format", 
                true);
    arg_desc.AddFlag(kArgOutputMd,
                "Output the MD tag in m4, paf or sam format",
                true);

    /// query filtering options
    arg_desc.SetCurrentGroup(kGroupQueryFiltering);

    arg_desc.AddDefaultKey(kArgDustFiltering, "DUST_options",
                "Filter query sequence with DUST " 
                "(Format: '" + kDfltArgApplyFiltering + "', " +
                "'level window linker', or '" + kDfltArgNoFiltering +
                "' to disable)",
                CArgDescriptions::eString,
                kDfltArgDustFiltering);

    arg_desc.AddFlag(kArgUseLCaseMasking, 
                "Use lower case filtering in query and subject sequence(s)?", true);

    /// restrict search or results
    arg_desc.SetCurrentGroup(kGroupRestrictSearch);

    arg_desc.AddOptionalKey(kArgPercentIdentity, "float_value",
                "Percent identity",
                CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgPercentIdentity, CArgAllow_Doubles(0.0, 100.0));

    arg_desc.AddOptionalKey(kArgQueryCovHspPerc, "float_value",
                "Percent query coverage per hsp",
                CArgDescriptions::eDouble);
    arg_desc.SetConstraint(kArgQueryCovHspPerc, CArgAllow_Doubles(0.0, 100.0));

    arg_desc.AddOptionalKey(kArgQueryCovHspRes, "int_value",
                "Residues query coverage per hsp",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgQueryCovHspRes, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgMaxHSPsPerSubject, "int_value",
                "Set maximum number of HSPs per subject sequence to save for each query",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxHSPsPerSubject, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddOptionalKey(kArgMaxTargetSequences, "num_sequences",
                "Maximum number of aligned sequences to keep \n"
                "(value of 5 or more is recommanded)\n"
                "Default = '" + NStr::IntToString(BLAST_HITLIST_SIZE) + "'",
                CArgDescriptions::eInteger);
    arg_desc.SetConstraint(kArgMaxTargetSequences, CArgAllowValuesGreaterThanOrEqual(1));

    arg_desc.AddFlag(kArgSubjectBestHit, "Turn on best hit per subject sequence", true);

    arg_desc.AddFlag(kArgSkipOverhang, "Skip overhangs in alignments", true);

    /// statistical options
    arg_desc.SetCurrentGroup(kGroupStatistics);

    arg_desc.AddOptionalKey(kArgEffSearchSpace, "int_value",
                "Effective length of the search space",
                CArgDescriptions::eInt8);
    arg_desc.SetConstraint(kArgEffSearchSpace, CArgAllowValuesGreaterThanOrEqual(0));

    arg_desc.AddOptionalKey(kArgDbSize, "num_letters",
                "Effective length of the database",
                CArgDescriptions::eInt8);

    /// miscellaneous options
    arg_desc.SetCurrentGroup(kGroupMiscellaneous);

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
    /// general search options
    if (args.Exist(kTask) && args[kTask].HasValue()) {
        string task_str = args[kTask].AsString();
        m_Options->align_task = string_to_task(task_str.c_str());
        if (m_Options->align_task == eHbnInvalidTask) {
            HBN_LOG("Invalid task: %s", task_str.c_str());
        }
    }

    if (args.Exist(kArgOutput) && args[kArgOutput].HasValue()) {
        m_Options->output = strdup(args[kArgOutput].AsString().c_str());
    }

    if (args.Exist(kArgEvalue) && args[kArgEvalue].HasValue()) {
        m_Options->expect_value = args[kArgEvalue].AsDouble();
    }

    if (args.Exist(kArgMismatch) && args[kArgMismatch].HasValue()) {
        m_Options->penalty = args[kArgMismatch].AsInteger();
    }

    if (args.Exist(kArgMatch) && args[kArgMatch].HasValue()) {
        m_Options->reward = args[kArgMatch].AsInteger();
    }

    if (args.Exist(kArgGapOpen) && args[kArgGapOpen].HasValue()) {
        m_Options->gap_open = args[kArgGapOpen].AsInteger();
    }

    if (args.Exist(kArgGapExtend) && args[kArgGapExtend].HasValue()) {
        m_Options->gap_extend = args[kArgGapExtend].AsInteger();
    }

    /// input query options
    if (args.Exist(kArgStrand) && args[kArgStrand].HasValue()) {
        string strand = args[kArgStrand].AsString();
        NStr::ToUpper(strand);
        if (strand == "BOTH") {
            m_Options->strand = F_R;
        } else if (strand == "PLUS") {
            m_Options->strand = FWD;
        } else if (strand == "MINUS") {
            m_Options->strand = REV;
        } else {
            HBN_ERR("Invalid strand value '%s'", strand.c_str());
        }
    }

    /// database options
    if (args.Exist(kArgDbDir) && args[kArgDbDir].HasValue()) {
        if (m_Options->db_dir) free((void*)m_Options->db_dir);
        string db_dir = args[kArgDbDir].AsString();
        m_Options->db_dir = strdup(db_dir.c_str());
    }

    if (args.Exist(kArgKeepDb))
        m_Options->keep_db = static_cast<bool>(args[kArgKeepDb]);

    if (args.Exist(kArgMinQuerySize) && args[kArgMinQuerySize].HasValue()) {
        m_Options->min_query_size = args[kArgMinQuerySize].AsInteger();
    }

    if (args.Exist(kArgMaxQueryVolSeqs) && args[kArgMaxQueryVolSeqs].HasValue()) {
        m_Options->max_query_vol_seqs = args[kArgMaxQueryVolSeqs].AsInteger();
    }

    if (args.Exist(kArgMaxQueryVolRes) && args[kArgMaxQueryVolRes].HasValue()) {
        m_Options->max_query_vol_res = args[kArgMaxQueryVolRes].AsInt8();
    }

    if (args.Exist(kArgMinSubjectSize) && args[kArgMinSubjectSize].HasValue()) {
        m_Options->min_subject_size = args[kArgMinSubjectSize].AsInteger();
    }

    if (args.Exist(kArgMaxSubjectVolSeqs) && args[kArgMaxSubjectVolSeqs].HasValue()) {
        m_Options->max_subject_vol_seqs = args[kArgMaxSubjectVolSeqs].AsInteger();
    }

    if (args.Exist(kArgMaxSubjectVolRes) && args[kArgMaxSubjectVolRes].HasValue()) {
        m_Options->max_subject_vol_res = args[kArgMaxSubjectVolRes].AsInt8();
    }

    /// ddf scoring options
    if (args.Exist(kArgKmerSize) && args[kArgKmerSize].HasValue()) {
        m_Options->kmer_size = args[kArgKmerSize].AsInteger();
    }    

    if (args.Exist(kArgKmerWindowSize) && args[kArgKmerWindowSize].HasValue()) {
        m_Options->kmer_window_size = args[kArgKmerWindowSize].AsInteger();
    }

    if (args.Exist(kArgBlockSize) && args[kArgBlockSize].HasValue()) {
        m_Options->block_size = args[kArgBlockSize].AsInteger();
    }

    if (args.Exist(kArgMinDDFS) && args[kArgMinDDFS].HasValue()) {
        m_Options->min_ddfs = args[kArgMinDDFS].AsInteger();
    }

    if (args.Exist(kArgMaxKmerOcc) && args[kArgMaxKmerOcc].HasValue()) {
        m_Options->max_kmer_occ = args[kArgMaxKmerOcc].AsInteger();
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

    /// output format
    if (args.Exist(kArgOutputFormat) && args[kArgOutputFormat].HasValue()) {
        string fmtstr = args[kArgOutputFormat].AsString();
        EOutputFormat outfmt = string_to_outfmt(fmtstr.c_str());
        if (outfmt == eInvalidFmt) {
            HBN_ERR("Invalid output format: %s", fmtstr.c_str());
        }
        m_Options->outfmt = outfmt;
    }

    if (args.Exist(kArgOutputCigar))
        m_Options->dump_cigar = static_cast<bool>(args[kArgOutputCigar]);

    if (args.Exist(kArgOutputMd))
        m_Options->dump_md = static_cast<bool>(args[kArgOutputMd]);

    /// query filtering options
    if (args.Exist(kArgDustFiltering) && args[kArgDustFiltering].HasValue()) {
        string duststr = args[kArgDustFiltering].AsString();
        CTempString cduststr(duststr);
        CTempString delim(" ");
        vector<string> dust_components;
        NStr::Split(cduststr, delim, dust_components);
        if (dust_components.size() == 1) {
            string s = dust_components[0];
            NStr::ToUpper(s);
            if (s == "NO") {
                m_Options->use_dust_masker = 0;
            } else if (s == "YES") {

            } else {
                HBN_ERR("Invalid value '%s' to argument '%s'", duststr.c_str(), kArgDustFiltering.c_str());
            }
        } else if (dust_components.size() == 3) {
            m_Options->dust_masker_level = NStr::StringToInt(dust_components[0]);
            m_Options->dust_masker_window = NStr::StringToInt(dust_components[1]);
            m_Options->dust_masker_linker = NStr::StringToInt(dust_components[2]);
        } else {
            HBN_ERR("Invalid value '%s' to argument '%s'", duststr.c_str(), kArgDustFiltering.c_str());
        }
    }

    if (args.Exist(kArgUseLCaseMasking))
        m_Options->use_lower_case_masker = static_cast<bool>(args[kArgUseLCaseMasking]);

    /// restrict search or results
    if (args.Exist(kArgPercentIdentity) && args[kArgPercentIdentity].HasValue()) {
        m_Options->perc_identity = args[kArgPercentIdentity].AsDouble();
    }

    if (args.Exist(kArgQueryCovHspPerc) && args[kArgQueryCovHspPerc].HasValue()) {
        m_Options->query_cov_hsp_perc = args[kArgQueryCovHspPerc].AsDouble();
    }

    if (args.Exist(kArgQueryCovHspRes) && args[kArgQueryCovHspRes].HasValue()) {
        m_Options->query_cov_hsp_res = args[kArgQueryCovHspRes].AsInteger();
    }

    if (args.Exist(kArgMaxHSPsPerSubject) && args[kArgMaxHSPsPerSubject].HasValue()) {
        m_Options->max_hsps_per_subject = args[kArgMaxHSPsPerSubject].AsInteger();
    }

    if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences].HasValue()) {
        m_Options->hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    }

    if (args.Exist(kArgSubjectBestHit))
        m_Options->keep_best_hsp_per_subject = static_cast<bool>(args[kArgSubjectBestHit]);

    if (args.Exist(kArgSkipOverhang))
        m_Options->skip_overhang = static_cast<bool>(args[kArgSkipOverhang]);

    /// statistical options
    if (args.Exist(kArgEffSearchSpace) && args[kArgEffSearchSpace].HasValue()) {
        m_Options->searchsp_eff = args[kArgEffSearchSpace].AsInt8();
    }

    if (args.Exist(kArgDbSize) && args[kArgDbSize].HasValue()) {
        m_Options->db_length = args[kArgDbSize].AsInt8();
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

void Init_HbnProgramOptions(HbnProgramOptions* opts, const EHbnTask task)
{
    /// general search options
    opts->align_task = task;
    opts->expect_value = BLAST_EXPECT_VALUE;
    opts->penalty = HBN_PENELTY;
    opts->reward = HBN_REWARD;
    opts->gap_open = HBN_GAP_OPEN;
    opts->gap_extend = HBN_GAP_EXTEND;

    /// input query options
    opts->strand = F_R;

    /// database options
    opts->db_dir = strdup(kDfltDbDir.c_str());
    opts->keep_db = kDfltKeepDb;
    opts->min_query_size = kDfltMinQuerySize;
    opts->max_query_vol_seqs = kDfltMaxQueryVolSeqs;
    opts->max_query_vol_res = kDfltMaxQueryVolRes;
    opts->min_subject_size = kDfltMinSubjectSize;
    opts->max_subject_vol_seqs = kDfltMaxSubjectVolSeqs;
    opts->max_subject_vol_res = kDfltMaxSubjectVolRes;

    /// ddf scoring options
    opts->kmer_size = kDfltKmerSize;
    opts->kmer_window_size = kDfltKmerWindowSize;
    opts->max_kmer_occ = kDfltMaxKmerOcc;
    opts->block_size = kDfltBlockSize;
    opts->min_ddfs = kDfltMinDDFS;

    /// mem chaining scoring options
    opts->memsc_kmer_size = kDfltMemScKmerSize;
    opts->memsc_kmer_window = kDfltMemScKmerWindow;
    opts->memsc_mem_size = kDfltMemScMemSize;
    opts->memsc_score = kDfltMemScMemScore;

    /// formatting options
    opts->outfmt = kDfltHbnOutfmt;
    opts->dump_cigar = kDfltOutputCigar;
    opts->dump_md = kDfltOutputMd;

    /// query filtering options
    opts->use_dust_masker = 1;
    opts->dust_masker_level = kDustLevel;
    opts->dust_masker_linker = kDustLinker;
    opts->dust_masker_window = kDustWindow;
    opts->use_lower_case_masker = 0;

    /// restrict search or results
    opts->perc_identity = 65.0;
    opts->query_cov_hsp_perc = 0.0;
    opts->query_cov_hsp_res = kDfltQueryCovHspRes;
    opts->max_hsps_per_subject = kDfltMaxHspsPerSubject;
    opts->hitlist_size = BLAST_HITLIST_SIZE;
    opts->keep_best_hsp_per_subject = 0;
    opts->skip_overhang = kDfltSkipOverhang;

    if (task == eHbnTask_pm) {
        opts->max_hsps_per_subject = 2;
        opts->hitlist_size = 100;
    } else {
        opts->max_hsps_per_subject = 5;
        opts->hitlist_size = 3;
    }

    /// statistical options
    opts->searchsp_eff = 0;
    opts->db_length = 0;
    opts->dbseq_num = 0;

    /// misc options
    opts->num_threads = kDfltNumThreads;
    opts->node_id = kDfltNodeId;
    opts->num_nodes = kDfltNumNodes;

    opts->query = NULL;
    opts->subject = NULL;
    opts->output = strdup(kDfltHbnOutput.c_str());
}

static EHbnTask
s_PreCheckCmdLineArgs(int argc, char* argv[], CArgDescriptions* arg_desc)
{
    string kProgram = FindProgramDisplayName(argv[0]);
    string kHbnUsage = kProgram + " [OPTIONS] query subject";
    EHbnTask task = kDfltTask;

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
        } else if (NStr::CompareCase(argv[i] + 1, kTask) == 0) {
            task = string_to_task(argv[i+1]);
            if (task == eHbnInvalidTask) {
                HBN_ERR("Invalid task: %s", argv[i+1]);
            }
        }
    }
    return task;
}

extern "C"
void ParseHbnProgramCmdLineArguments(int argc, char* argv[], HbnProgramOptions* opts)
{
    TBlastCmdLineArgs arg_list;

    /// setup description
    CRef<IBlastCmdLineArgs> arg;
    string kProgram = FindProgramDisplayName(argv[0]);
    string kProgramDescription("Nucleotide-Nucleotide sequence alignment toolkit");
    arg.reset(new CProgramDescriptionArgs(kProgram, kProgramDescription));
    arg_list.push_back(arg);

    arg.reset(new CommandLineArguments(opts));
    arg_list.push_back(arg);

    unique_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);
    NON_CONST_ITERATE(TBlastCmdLineArgs, arg_iter, arg_list) {
        (*arg_iter)->SetArgumentDescriptions(*arg_desc);
    }

    /// examine trivial arguments (-help, -h, -version)
    EHbnTask task = s_PreCheckCmdLineArgs(argc, argv, arg_desc.get());
    Init_HbnProgramOptions(opts, task);

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
        HBN_ERR("The query and subject must be specified");
    } else if (argc - argv_idx > 2) {
        string err = "Too many query and subject values: '";
        for (int i = argv_idx; i < argc; ++i) {
            err += argv[i];
            if (i != argc - 1) err += ' ';
        }
        err += "'";
        HBN_ERR("%s", err.c_str());
    } else {
        opts->query = argv[argv_idx];
        opts->subject = argv[argv_idx + 1];
    }
}

#define os_one_option_value(name, value) os << '-' << name << ' ' << value << ' '
#define os_one_flag_option(name) os << '-' << name << ' ';

extern "C"
char* HbnProgramOptions2String(const HbnProgramOptions* opts)
{
    ostringstream os;

    /// general search options
    os_one_option_value(kTask, hbn_task_names[opts->align_task]);
    os_one_option_value(kArgEvalue, opts->expect_value);
    os_one_option_value(kArgMismatch, opts->penalty);
    os_one_option_value(kArgMatch, opts->reward);
    os_one_option_value(kArgGapOpen, opts->gap_open);
    os_one_option_value(kArgGapExtend, opts->gap_extend);

    /// input query options
    const char* strand_name[] = { "plus", "minus", "both" };
    os_one_option_value(kArgStrand, strand_name[opts->strand]);

    /// database options
    os_one_option_value(kArgDbDir, opts->db_dir);
    if (opts->keep_db) os_one_flag_option(kArgKeepDb);
    if (opts->min_query_size) os_one_option_value(kArgMinQuerySize, opts->min_query_size);
    if (opts->max_query_vol_seqs != kDfltMaxQueryVolSeqs) os_one_option_value(kArgMaxQueryVolSeqs, opts->max_query_vol_seqs);
    string size_str = NStr::UInt8ToString_DataSize(opts->max_query_vol_res);
    os_one_option_value(kArgMaxQueryVolRes, size_str);
    if (opts->min_subject_size) os_one_option_value(kArgMinSubjectSize, opts->min_subject_size);
    if (opts->max_subject_vol_seqs != kDfltMaxSubjectVolSeqs) os_one_option_value(kArgMaxSubjectVolSeqs, opts->max_subject_vol_seqs);
    size_str = NStr::UInt8ToString_DataSize(opts->max_subject_vol_res);
    os_one_option_value(kArgMaxSubjectVolRes, size_str);

    /// ddf scoring
    os_one_option_value(kArgKmerSize, opts->kmer_size);
    os_one_option_value(kArgKmerWindowSize, opts->kmer_window_size);
    os_one_option_value(kArgMaxKmerOcc, opts->max_kmer_occ);
    os_one_option_value(kArgBlockSize, opts->block_size);
    os_one_option_value(kArgMinDDFS, opts->min_ddfs);

    /// mem chaining scoring options
    os_one_option_value(kArgMemScKmerSize, opts->memsc_kmer_size);
    os_one_option_value(kArgMemScKmerWindow, opts->memsc_kmer_window);
    os_one_option_value(kArgMemScMemSize, opts->memsc_mem_size);
    os_one_option_value(kArgMemScMemScore, opts->memsc_score);

    /// output format
    os_one_option_value(kArgOutputFormat, opts->outfmt);

    /// query filtering options
    if (opts->use_dust_masker) {
        os << '-' << kArgDustFiltering << ' '
           << "'" << opts->dust_masker_level 
           << " " << opts->dust_masker_window
           << " " << opts->dust_masker_linker
           << "' ";
    }
    if (opts->use_lower_case_masker) os_one_flag_option(kArgUseLCaseMasking);

    /// restrict search of results
    os_one_option_value(kArgPercentIdentity, opts->perc_identity);
    os_one_option_value(kArgQueryCovHspPerc, opts->query_cov_hsp_perc);
    os_one_option_value(kArgQueryCovHspRes, opts->query_cov_hsp_res);
    os_one_option_value(kArgMaxHSPsPerSubject, opts->max_hsps_per_subject);
    os_one_option_value(kArgMaxTargetSequences, opts->hitlist_size);
    if (opts->keep_best_hsp_per_subject) os_one_flag_option(kArgSubjectBestHit);

    /// misc options
    os_one_option_value(kArgNumThreads, opts->num_threads);
    os << '-' << kArgGrid << ' ' << opts->node_id << ' ' << opts->num_nodes << ' ';

    size_str = os.str();
    return strdup(size_str.c_str());
}

END_NCBI_SCOPE
