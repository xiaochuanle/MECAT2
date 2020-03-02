ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := libfsa.a
SOURCES  := argument_parser.cpp getopt.c logger.cpp overlap.cpp read_store.cpp sequence.cpp\
             utility.cpp fasta_reader.cpp fastq_reader.cpp overlap_store.cpp \
   			 ./simple_align.cpp overlap_filter.cpp overlap_stat.cpp

TGT_CXXFLAGS := -U_GLIBCXX_PARALLEL -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64 
SRC_INCDIRS  := . 

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lontcns
TGT_PREREQS := libontcns.a

SUBMAKEFILES :=
