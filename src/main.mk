ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libhbn.a

SOURCES      := \
	./corelib/build_db.c \
	./corelib/cmd_arg.c \
	./corelib/cstr_util.c \
	./corelib/db_format.c \
	./corelib/fasta.c \
	./corelib/gapped_candidate.c \
	./corelib/hbn_aux.c \
	./corelib/hbn_format.c \
	./corelib/hbn_hit.c \
	./corelib/hbn_package_version.c \
	./corelib/kstring.c \
	./corelib/line_reader.c \
	./corelib/m4_record.c \
	./corelib/name2id_map.c \
	./corelib/partition_aux.c \
	./corelib/raw_reads.c \
	./corelib/seqdb_summary.c \
	./corelib/seqdb.c \
	./corelib/seq_tag.c \
	./corelib/seq_tag_report.cpp \
	./corelib/small_object_alloc.c \
	./corelib/string2hsp.c \
	./algo/chain_dp.c \
	./algo/diff_gapalign.cpp \
	./algo/hash_list_bucket_sort.c \
	./algo/kalloc.c \
	./algo/ksw2_extd2_sse.c \
	./algo/ksw2_extz2_sse.c \
	./algo/ksw2_wrapper.c \
	./algo/hbn_lookup_table.c \
	./algo/hbn_traceback_aux.c \
	./algo/word_finder.c \
	./ncbi_blast/c_ncbi_blast_aux.c \
	./ncbi_blast/ncbi_blast_aux.cpp \
	./ncbi_blast/cmdline_args/blast_args.cpp \
	./ncbi_blast/cmdline_args/cmdline_flags.cpp \
	./ncbi_blast/cmdline_args/format_flags.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_allow.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_desc.cpp \
	./ncbi_blast/cmdline_args/ncbiargs_types.cpp \
	./ncbi_blast/str_util/ncbistr_util.cpp \
	./ncbi_blast/str_util/ncbistr.cpp \
	./ncbi_blast/str_util/str_cmp.cpp \
	./ncbi_blast/str_util/numeric_str_interconv.cpp \
	./ncbi_blast/str_util/str_util.cpp \
	./ncbi_blast/setup/blast_encoding.c \
	./ncbi_blast/setup/blast_hits.c \
	./ncbi_blast/setup/blast_message.c \
	./ncbi_blast/setup/blast_options.c \
	./ncbi_blast/setup/blast_stat.c \
	./ncbi_blast/setup/blast_parameters.c \
	./ncbi_blast/setup/blast_program.c \
	./ncbi_blast/setup/blast_types.cpp \
	./ncbi_blast/setup/boost_erf.c \
	./ncbi_blast/setup/hsp2string.cpp \
	./ncbi_blast/setup/ncbi_math.c \
	./ncbi_blast/setup/blast_query_info.c \
	./ncbi_blast/setup/blast_sequence_blk.c \
	./ncbi_blast/setup/gapinfo.c

SRC_INCDIRS  :=

SUBMAKEFILES := \
	./app/mecat2pcan/pcan.mk \
	./app/map/main.mk \
	./app/mecat2cns/main.mk \
	./app/hbndb/viewhbndb.mk \
	./app/mecat2trim/1_largest_cover_range/main.mk \
	./app/mecat2trim/2_split_reads/main.mk \
	./app/mecat2trim/3_trim_bases/main.mk \
	./app/mecat2cns/main.mk \
	./app/mecat2pm4/main.mk \
	./app/mecat2extseqs/main.mk \
	./app/fsa/fsa.mk	\
	./app/fsa/filter.mk \
	./app/fsa/assemble.mk \
	./app/fsa/bridge.mk	\
	./app/fsa/rd_stat.mk \
	./pipeline/main.mk \