ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2cns
SOURCES  := \
	cmdline_args.cpp \
	cns_aux.c \
	cns_one_part.c \
	cns_one_read.cpp \
	cns_options.c \
	fccns_align_tag.c \
	fccns_aux.c \
	fccns.c \
	hbn_task_struct.c \
	main.c \
	raw_reads_reader.c \

SRC_INCDIRS  := .

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lhbn
TGT_PREREQS := libhbn.a

SUBMAKEFILES :=