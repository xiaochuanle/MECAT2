ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := v2pm4
SOURCES  := pm4_aux.c pm4_main.c ../common/ontcns_aux.c ../common/ontcns_defs.c ../common/oc_assert.c ../klib/kstring.c m4_record.c 

SRC_INCDIRS  := .

TGT_LDFLAGS :=
TGT_LDLIBS  :=
TGT_PREREQS :=

SUBMAKEFILES :=
