ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2map
SOURCES  := \
	cmdline_args.cpp \
	hbn_align_one_volume.c \
	hbn_build_seqdb.c \
	hbn_extend_subseq_hit.c \
	hbn_find_subseq_hit.c \
	hbn_job_control.c \
	hbn_options.c \
	hbn_subseq_hit.c \
	hbn_task_struct.c \
	main.c \
	mecat_results.c \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=