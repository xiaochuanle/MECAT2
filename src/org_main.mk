ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET       := libmecat.a

SOURCES      := common/alignment.cpp \
		common/buffer_line_iterator.cpp \
		common/defs.cpp \
		common/diff_gapalign.cpp \
		common/fasta_reader.cpp \
		common/gapalign.cpp \
		common/lookup_table.cpp \
		common/packed_db.cpp \
		common/sequence.cpp \
		common/split_database.cpp \
		common/xdrop_gapalign.cpp

SRC_INCDIRS  := common \

SUBMAKEFILES := mecat2pw/pw.mk \
		mecat2ref/mecat2ref.mk \
		mecat2cns/mecat2cns.mk \
		filter_reads/filter_reads.mk \
		./mecat2asm/v2pm/v2_make_volumes.mk \
		./mecat2asm/v2pm/v2_asmpm.mk \
		./mecat2asm/v2trim/pm4.mk \
		./mecat2asm/v2trim/largest_cover_range_main.mk \
		./mecat2asm/v2trim/split_reads_main.mk \
		./mecat2asm/v2trim/trim_bases.mk \
		./mecat2asm/v2elr/extract_sequences.mk \
		./fsa/fsa.mk	\
		./fsa/filter.mk \
		./fsa/assemble.mk \
		./fsa/bridge.mk	\
		./fsa/rd_stat.mk \
		./pipeline/main.mk
