ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := v2lcr
SOURCES  := pm4_aux.c largest_cover_range.c ../klib/kstring.c ../common/ontcns_aux.c range_list.c m4_record.c ../common/oc_assert.c largest_cover_range_main.c

SRC_INCDIRS  := .

TGT_LDFLAGS :=
TGT_LDLIBS  :=
TGT_PREREQS :=

SUBMAKEFILES :=
