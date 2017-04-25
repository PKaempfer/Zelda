AUTOMAKE_OPTIONS= foreign

CFLAGS = -Wall -Wextra -Wno-sign-compare -Wno-unused-result -Wno-strict-aliasing -Wno-unused-but-set-parameter -Wno-unused-but-set-variable
CFLAGS += -I.. -I../include -I/usr/include 
CFLAGS += -g -O3 -D_POSIX_C_SOURCE=199309L -std=c99
LIBS = -pthread -lrt -lm

bin_PROGRAMS = Zelda $(LINKFLAGS)

Zelda_SOURCES = Zelda.c
Zelda_SOURCES +=  DBGraph.c DBGraph_oa.c
Zelda_SOURCES += FileReader.c  kmer.c kmerHash.c kmerHash_oa.c readDB.c
Zelda_SOURCES += read_filter.c DBGraph_error.c
Zelda_SOURCES += DBGraph_reduced.c DBGraph_reduced_strong.c
Zelda_SOURCES += DBGraph_stringer_2.c DBGraph_stringer_refine.c
Zelda_SOURCES += DBGraph_scaffold.c
Zelda_SOURCES += ConsensusCaller.c ConsensusCaller2.c CC.c CC2.c
Zelda_SOURCES += misc_hash.c

#Zelda_OBJECTS = $(Zelda_SOURCES:.c=.o)

#%.o: src/%.c
