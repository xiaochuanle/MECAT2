libfsa_src = argument_parser.cpp getopt.c logger.cpp overlap.cpp read_store.cpp sequence.cpp\
             utility.cpp fasta_reader.cpp fastq_reader.cpp overlap_store.cpp \
			 overlap_filter.cpp overlap_stat.cpp simple_align.cpp

libfsa_obj:=$(patsubst %.cpp,%.o,$(libfsa_src))
libfsa_obj:=$(patsubst %.c,%.o,$(libfsa_obj))

fsa_ol_filter_src = fsa_ol_filter.cpp

fsa_ol_filter_obj:=$(patsubst %.cpp,%.o,$(fsa_ol_filter_src))
fsa_ol_filter_obj:=$(patsubst %.c,%.o,$(fsa_ol_filter_obj))

fsa_ol_improve_src = overlap_improve.cpp 

fsa_ol_improve_obj:=$(patsubst %.cpp,%.o,$(fsa_ol_improve_src))
fsa_ol_improve_obj:=$(patsubst %.c,%.o,$(fsa_ol_improve_obj))

fsa_ol_compare_src = overlap_compare.cpp fsa_ol_compare.cpp

fsa_ol_compare_obj:=$(patsubst %.cpp, %.o, $(fsa_ol_compare_src))
fsa_ol_compare_obj:=$(patsubst %.c, %.o, $(fsa_ol_compare_obj))

fsa_ol_stat_src = fsa_ol_stat.cpp
fsa_ol_stat_obj:=$(patsubst %.cpp, %.o, $(fsa_ol_stat_src))

fsa_ol_show_src = fsa_ol_show.cpp overlap_show.cpp
fsa_ol_show_obj:=$(patsubst %.cpp,%.o,$(fsa_ol_show_src))

fsa_assemble_src = fsa_assemble.cpp assembly.cpp graph.cpp string_graph.cpp path_graph.cpp
   

fsa_assemble_obj:=$(patsubst %.cpp,%.o,$(fsa_assemble_src))
fsa_assemble_obj:=$(patsubst %.c,%.o,$(fsa_assemble_obj))

fsa_ctg_bridge_src = fsa_ctg_bridge.cpp contig_bridge.cpp contig_graph.cpp contig_link.cpp contig_link_store.cpp

fsa_ctg_bridge_obj:=$(patsubst %.cpp,%.o,$(fsa_ctg_bridge_src))
fsa_ctg_bridge_obj:=$(patsubst %.c,%.o,$(fsa_ctg_bridge_obj))

fsa_rd_stat_src = fsa_rd_stat.cpp read_stat.cpp
fsa_rd_stat_obj:=$(patsubst %.cpp, %.o, $(fsa_rd_stat_src))


ifeq ($(BUILDDEBUG), 1)
	CFLAGS =  -Wall -O0  -D_FILE_OFFSET_BITS=64 -g
	CXXFLAGS = -std=c++11 -Wall -O0  -D_FILE_OFFSET_BITS=64 -g
else
	CFLAGS =  -Wall -O3  -D_FILE_OFFSET_BITS=64
	CXXFLAGS = -std=c++11 -Wall -O3 -D_FILE_OFFSET_BITS=64
endif

LDFLAGS = -pthread 


obj:=$(libfsa_obj) $(fsa_ol_filter_obj) $(fsa_ol_improve_obj) $(fsa_ol_compare_obj) \
	$(fsa_ol_show_obj) $(fsa_ol_stat_obj) $(fsa_assemble_obj) $(fsa_ctg_bridge_obj) \
	$(fsa_rd_stat_obj)

dep:=$(patsubst %.o,%.d, $(obj))


bin = libfsa.a fsa_ol_filter  fsa_ol_improve fsa_ol_compare fsa_ol_stat fsa_ol_show fsa_assemble fsa_ctg_bridge \
      fsa_rd_stat

all: $(bin)

ifneq ($(MAKECMDGOALS),clean)
-include $(dep)
endif

libfsa.a: $(libfsa_obj)
	ar -r $@ $^

fsa_ol_filter:$(fsa_ol_filter_obj) libfsa.a
	g++ $(LDFLAGS) $(fsa_ol_filter_obj) libfsa.a -o $@

fsa_ol_improve:$(fsa_ol_improve_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

fsa_ol_compare:$(fsa_ol_compare_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@
	
fsa_ol_stat:$(fsa_ol_stat_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

fsa_ol_show:$(fsa_ol_show_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

fsa_assemble:$(fsa_assemble_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

fsa_ctg_bridge:$(fsa_ctg_bridge_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

fsa_rd_stat:$(fsa_rd_stat_obj) libfsa.a
	g++ $(LDFLAGS) $^ -o $@

%.o:%.c
	gcc $(CFLAGS) -c -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.o)" -o "$@" "$<"

%.o:%.cpp
	g++ $(CXXFLAGS) -c -MMD -MP -MF "$(@:%.o=%.d)" -MT "$(@:%.o=%.o)" -o "$@" "$<"

clean:
	rm $(obj) $(dep) $(bin) -f

