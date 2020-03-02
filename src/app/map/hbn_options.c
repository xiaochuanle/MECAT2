#include "hbn_options.h"

const char* outfmt_names[eInvalidFmt] = {
    "seqid",
    "seqidx",
    "subseq",
    "subseqx",
    "m4",
    "m4x",
    "paf",
    "sam"
};

EOutputFormat
string_to_outfmt(const char* str)
{
    for (int i = 0; i < eInvalidFmt; ++i) {
        if (strcmp(outfmt_names[i], str) == 0) return i;
    }
    return eInvalidFmt;
}

const char* hbn_task_names[eHbnInvalidTask] = {
    "pm",
    "rm"
};

EHbnTask
string_to_task(const char* str)
{
    for (int i = 0; i < eHbnInvalidTask; ++i) {
        if (strcmp(hbn_task_names[i], str) == 0) return i;
    }
    return eHbnInvalidTask;
}

HbnProgramOptions*
HbnProgramOptionsFree(HbnProgramOptions* opts)
{
    if (opts->db_dir) free((void*)opts->db_dir);
    opts->db_dir = NULL;

    if (opts->output) free((void*)opts->output);
    opts->output = NULL;

    free(opts);

    return NULL;
}