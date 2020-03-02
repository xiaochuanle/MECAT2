#include "hbn_build_seqdb.h"

#include "../../corelib/build_db.h"
#include "../../corelib/seqdb_summary.h"

#include <errno.h>
#include <sys/stat.h>

void 
hbn_build_seqdb(const HbnProgramOptions* opts, 
    const char* query_db_title, 
    const char* subject_db_title)
{
    if ((access(opts->db_dir, F_OK) != 0)
        &&
        (mkdir(opts->db_dir, S_IRWXU) != 0)) {
        HBN_ERR("Failed to create directory %s: %s", opts->db_dir, strerror(errno));
    }

    if (!hbndb_is_built(opts->query,
            opts->db_dir,
            query_db_title,
            opts->min_query_size,
            opts->max_query_vol_seqs,
            opts->max_query_vol_res,
            FALSE)) {
        build_db(opts->query,
            opts->db_dir,
            query_db_title,
            opts->min_query_size,
            opts->max_query_vol_seqs,
            opts->max_query_vol_res,
            FALSE);
        hbndb_make_built(opts->query,
            opts->db_dir,
            query_db_title,
            opts->min_query_size,
            opts->max_query_vol_seqs,
            opts->max_query_vol_res,
            FALSE);
        hbn_assert(hbndb_is_built(opts->query,
            opts->db_dir,
            query_db_title,
            opts->min_query_size,
            opts->max_query_vol_seqs,
            opts->max_query_vol_res,
            FALSE));
    }

    BOOL query_and_subject_are_the_same = (strcmp(opts->query, opts->subject) == 0);
    if ((!query_and_subject_are_the_same)
        &&
        (!hbndb_is_built(opts->subject,
            opts->db_dir,
            subject_db_title,
            opts->min_subject_size,
            opts->max_subject_vol_seqs,
            opts->max_subject_vol_res,
            FALSE))) {
        build_db(opts->subject,
            opts->db_dir,
            subject_db_title,
            opts->min_subject_size,
            opts->max_subject_vol_seqs,
            opts->max_subject_vol_res,
            FALSE);
        hbndb_make_built(opts->subject,
            opts->db_dir,
            subject_db_title,
            opts->min_subject_size,
            opts->max_subject_vol_seqs,
            opts->max_subject_vol_res,
            FALSE);
        hbn_assert(hbndb_is_built(opts->subject,
            opts->db_dir,
            subject_db_title,
            opts->min_subject_size,
            opts->max_subject_vol_seqs,
            opts->max_subject_vol_res,
            FALSE));
    }
}