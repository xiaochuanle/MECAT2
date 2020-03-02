#include "partition_aux.h"

#include <pthread.h>

void
make_partition_name(const char* data_dir, const char* prefix, const int pid, char path[])
{
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (prefix) strcat(path, prefix);
    char buf[64];
    u64_to_fixed_width_string_r(pid, buf, HBN_DIGIT_WIDTH);
    strcat(path, buf);
}

void
dump_partition_count(const char* data_dir, const char* prefix, const int np)
{
    char path[HBN_MAX_PATH_LEN];
    path[0] = '\0';
    if (data_dir) sprintf(path, "%s/", data_dir);
    if (prefix) {
        strcat(path, prefix);
        strcat(path, ".");
    }
    strcat(path, "np");
    hbn_dfopen(out, path, "w");
    fprintf(out, "%d\n", np);
    hbn_fclose(out);
}

int
load_partition_count(const char* data_dir, const char* prefix)
{
    char path[HBN_MAX_PATH_LEN];
    if (prefix) {
        sprintf(path, "%s/%s.np", data_dir, prefix);
    } else {
        sprintf(path, "%s/np", data_dir);
    }
    int np;
    hbn_dfopen(in, path, "r");
    HBN_SCANF(fscanf, in, 1, "%d", &np);
    hbn_fclose(in);
    return np;    
}

void* load_part_records(const char* path, const size_t record_size, size_t* n_record)
{
    size_t s = hbn_file_size(path);
    if (!s) return NULL;
    hbn_assert(s % record_size == 0);
    void* a = malloc(s);
    s /= record_size;
    hbn_dfopen(in, path, "rb");
    hbn_fread(a, record_size, s, in);
    hbn_fclose(in);
    *n_record = s;
    return a;    
}

typedef struct {
    FILE** out_list;
    int n;
} RecordWriter;

RecordWriter*
can_writer_new(const char* wrk_dir, const int pid_from, const int pid_to)
{
    RecordWriter* w = (RecordWriter*)malloc( sizeof(RecordWriter) );
    int n = pid_to - pid_from;
    w->out_list = (FILE**)malloc( sizeof(FILE*) * n );
    w->n = n;
    char path[HBN_MAX_PATH_LEN];
    for (int i = pid_from; i < pid_to; ++i) {
        make_partition_name(wrk_dir, DEFAULT_PART_PREFIX, i, path);
        int fid = i - pid_from;
        hbn_fopen(w->out_list[fid], path, "wb");
    }
    return w;
}

RecordWriter*
can_writer_free(RecordWriter* w)
{
    for (int i = 0; i < w->n; ++i) hbn_fclose(w->out_list[i]);
    free(w->out_list);
    free(w);
    return NULL;
}

static size_t
load_records(void* array, const size_t size, const size_t count, FILE* in, pthread_mutex_t* in_lock)
{
    size_t n = 0;
    pthread_mutex_lock(in_lock);
    n = fread(array, size, count, in);
    pthread_mutex_unlock(in_lock);
    return n;
}

typedef struct {
    FILE* in;
    pthread_mutex_t*            in_lock;
    RecordWriter*               w;
    pthread_mutex_t*            w_lock;
    qid_extract_func           get_qid;
    sid_extract_func           get_sid;
    change_record_roles_func    change_roles;
    normolize_sdir_func        normalise_sdir;
    record_sort_func           sort_records;
    int                         min_read_id;
    int                         max_read_id;
    int                         batch_size;
    size_t                      record_size;
} PartRecordData;

#define id_in_range(id, L, R) ((id) >= (L) && (id) < (R))

static void*
pcan_worker(void* param)
{
    PartRecordData* data = (PartRecordData*)(param);
    const size_t S = U64_ONE * 256 * 1024 * 1024;
    const size_t N = S / data->record_size;
    const size_t N2 = 2 * N;
    void* a = malloc( N2 * data->record_size );
    void* record = malloc( data->record_size );
    kv_dinit(vec_size_t, idx_range);

    size_t n, m;
    while ((n = load_records(a, data->record_size, N, data->in, data->in_lock))) {
        m = 0;
        for (size_t i = 0; i < n; ++i) {
            void* e = a + data->record_size * i;
            int qid = (*data->get_qid)(e);
            int sid = (*data->get_sid)(e);
            int r = id_in_range(qid, data->min_read_id, data->max_read_id) || 
                    id_in_range(sid, data->min_read_id, data->max_read_id);
            if (r) {
                (*data->normalise_sdir)(e);
                if (i > m) memcpy(a + m * data->record_size, a + i * data->record_size, data->record_size);
                ++m;
            }
        }
        if (m == 0) continue;
        n = m;
        for (size_t i = 0; i < m; ++i) {
            void* e = a + i * data->record_size;
            int sid_is_in = 0;
            int sid = (*data->get_sid)(e);
            if (id_in_range(sid, data->min_read_id, data->max_read_id)) sid_is_in = 1;
            int qid = (*data->get_qid)(e);
            if (id_in_range(qid, data->min_read_id, data->max_read_id)) {
                (*data->change_roles)(e, record);
                (*data->normalise_sdir)(record);
                if (sid_is_in) {
                    memcpy(a + n * data->record_size, record, data->record_size);
                    ++n;
                } else {
                    memcpy(a + i * data->record_size, record, data->record_size);
                }
            }
        }

        (*data->sort_records)(n, a);
        kv_clear(idx_range);
        kv_push(size_t, idx_range, 0);
        size_t i = 0;
        while (i < n) {
            void* e = a + i * data->record_size;
            const int sid = (*data->get_sid)(e);
            const int pid = sid / data->batch_size;
            const int sid_from = pid * data->batch_size;
            const int sid_to = sid_from + data->batch_size;
            size_t j = i + 1;
            while (j < n) {
                e = a + j * data->record_size;
                int id = (*data->get_sid)(e);
                if (id >= sid_to) break;
                ++j;
            }
            kv_push(size_t, idx_range, j);
            i = j;
        }
        hbn_assert(kv_back(idx_range) == n);

        pthread_mutex_lock(data->w_lock);
        for (i = 0; i < kv_size(idx_range) - 1; ++i) {
            size_t from = kv_A(idx_range, i);
            size_t to = kv_A(idx_range, i + 1);
            m = to - from;
            void* e = a + from * data->record_size;
            int sid = (*data->get_sid)(e);
            int fid = (sid - data->min_read_id) / data->batch_size;
            hbn_assert(fid < data->w->n);
            hbn_fwrite(a + from * data->record_size, data->record_size, m, data->w->out_list[fid]);
        }
        pthread_mutex_unlock(data->w_lock);
    }

    free(record);
    free(a);
    kv_destroy(idx_range);
    return NULL;
}

void
part_record_main(const char* part_wrk_dir,
    const char* record_path,
    const int num_batches,
    const int batch_size,
    const int num_threads,
    const int num_dumpped_files,
    const size_t record_size,
    qid_extract_func           get_qid,
    sid_extract_func           get_sid,
    change_record_roles_func    change_roles,
    normolize_sdir_func        normalise_sdir,
    record_sort_func           sort_records)
{
    pthread_mutex_t in_lock;
    pthread_mutex_init(&in_lock, NULL);
    pthread_mutex_t out_lock;
    pthread_mutex_init(&out_lock, NULL);
    pthread_t job_ids[num_threads];

    for (int fid = 0; fid < num_batches; fid += num_dumpped_files) {
        int sfid = fid;
        int efid = hbn_min(sfid + num_dumpped_files, num_batches);
        int min_read_id = sfid * batch_size;
        int max_read_id = efid * batch_size;
        RecordWriter* w = can_writer_new(part_wrk_dir, sfid, efid);
        hbn_dfopen(record_in, record_path, "rb");
        PartRecordData can_data = {
            record_in,
            &in_lock,
            w,
            &out_lock,
            get_qid,
            get_sid,
            change_roles,
            normalise_sdir,
            sort_records,
            min_read_id,
            max_read_id,
            batch_size,
            record_size
        };
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(job_ids + i, NULL, pcan_worker, &can_data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(job_ids[i], NULL);
        }
        hbn_fclose(record_in);
        can_writer_free(w);
    }
}