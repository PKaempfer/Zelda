#SHELL=/bin/bash -i
CC = gcc  
#CC = clang-3.5
#CFLAGS = -Wall -Wextra -I../include -I./Jellyfish-2.2.4/include -I. -I../uthash/uthash-master/src -I/usr/local/include -mcpu=power8 -mcx16
CFLAGS = -D_POSIX_C_SOURCE=199309L -std=c99 -I./include -I/usr/include  #-mcpu=power8 -mcx16 -D_POSIX_C_SOURCE=199309L 
CFLAGS += -Wall -Wextra  -Wno-unused-result -Wno-strict-aliasing -Wno-sign-compare #-Wno-unused-but-set-variable
#CFLAGS += -Wno-unused-parameter 
#CFLAGS += -Wno-unused-but-set-parameter 

LINKFLAGS = -pthread -lm #-lrt 

#CFLAGS += -O1 -fno-strict-aliasing -Wno-sign-compare -Wno-comment #-Wno-unused-but-set-variable
CFLAGS += -g -O3 
# Adress sanitizer seeks for memory leaeks in address room
#CFLAGS += -g -O1
#CFLAGS += -fsanitize=address
#CFLAGS += -fsanitize=leak

SOURCES =  src/DBGraph.c src/DBGraph_oa.c
SOURCES += src/FileReader.c  src/kmer.c src/kmerHash.c src/kmerHash_oa.c src/readDB.c
SOURCES += src/read_filter.c src/DBGraph_error.c
SOURCES += src/DBGraph_reduced.c src/DBGraph_reduced_strong.c
SOURCES += src/DBGraph_stringer_2.c src/DBGraph_stringer_refine.c
SOURCES += src/DBGraph_scaffold.c src/DBGraph_scaffold_touring.c
SOURCES += src/ConsensusCaller.c src/ConsensusCaller2.c src/CC.c src/CC2.c
SOURCES += src/misc_hash.c


OBJECTS = $(SOURCES:.c=.o)
BINDIR = bin

TODAY=`date +%Y-%m-%d`

all: Zelda

%.o: src/%.c
	$(CC) $(CFLAGS) -c $<
Zelda: Zelda.o $(BINDIR)
	$(CC) $(CFLAGS) src/Zelda.o $(OBJECTS) -o $(BINDIR)/Zelda $(LINKFLAGS) -v
Zelda.o: $(OBJECTS)
	$(CC) $(CFLAGS) -c src/Zelda.c -o src/Zelda.o
Blast:
	$(CC) $(CFLAGS) -o $(BINDIR)/BlastValidator src/BlastValidator.c $(LINKFLAGS) -v
$(BINDIR):
	if [ ! -d "./$(BINDIR)" ]; then\
		mkdir $(BINDIR);\
	fi
clean:
	rm -f src/*.o Zelda
