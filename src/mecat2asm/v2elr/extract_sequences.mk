ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := mecat2elr
SOURCES  := extract_sequences.c ../common/ontcns_aux.c

SRC_INCDIRS  := .

TGT_LDFLAGS :=
TGT_LDLIBS  :=
TGT_PREREQS :=

SUBMAKEFILES :=
