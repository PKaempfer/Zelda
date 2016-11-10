#SHELL=/bin/bash -i
#CC = gcc  
CC = clang-3.5
#CFLAGS = -Wall -Wextra -I../include -I./Jellyfish-2.2.4/include -I. -I../uthash/uthash-master/src -I/usr/local/include -mcpu=power8 -mcx16
CFLAGS = -Wall -Wextra -I./include -I/usr/local/include #-mcpu=power8 -mcx16 

LINKFLAGS = -pthread -lrt -lm 

#CFLAGS += -O3 -fno-strict-aliasing -Wno-sign-compare -Wno-unused-but-set-variable -Wno-comment
CFLAGS += -g -O1 -Wno-sign-compare
# Adress sanitizer seeks for memory leaeks in address room
CFLAGS += -fsanitize=address
CFLAGS += -fsanitize=leak

SOURCES =  src/DBGraph.c src/DBGraph_oa.c
SOURCES += src/FileReader.c  src/kmer.c src/kmerHash.c src/kmerHash_oa.c src/readDB.c 
SOURCES += src/DBGraph_error.c
SOURCES += src/DBGraph_reduced.c src/DBGraph_reduced_strong.c
SOURCES += src/DBGraph_stringer_2.c src/DBGraph_stringer_refine.c
SOURCES += src/DBGraph_scaffold.c
SOURCES += src/ConsensusCaller.c src/ConsensusCaller2.c
SOURCES += src/misc_hash.c

OBJECTS = $(SOURCES:.c=.o) 

TODAY=`date +%Y-%m-%d`

all: Zelda

%.o: src/%.c
	$(CC) $(CFLAGS) -c $<

Zelda: Zelda.o
	$(CC) $(CFLAGS) src/Zelda.o $(OBJECTS) -o bin/Zelda $(LINKFLAGS) -v
Zelda.o: $(OBJECTS)
	$(CC) $(CFLAGS) -c src/Zelda.c -o src/Zelda.o
clean:
	rm -f src/*.o Zelda