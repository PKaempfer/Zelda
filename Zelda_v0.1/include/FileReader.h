/*
 * FileReader.h
 *
 *  Created on: Oct 13, 2014
 *      Author: kaempfpp
 */

#ifndef FILEREADER_H_
#define FILEREADER_H_

#include <pthread.h>

#define CHAIN_0 0
#define CHAIN_5 1
#define CHAIN_10 2
#define CHAIN_20 3
#define CHAIN_100 4
#define CHAIN_MAX 5
//#define UT_HASH

/** Relative orientation of the reads in the graph to each other */
// Insert size:    |-----------|
#define FR 0 	//  --->   <--- (forward - reverse) -> Default
#define RF 1	//  <---   ---> (reverse - forward)
#define FF 2	//  --->   ---> (forward - forward) or
#define	 RR 2	// 	<---   <--- (reverse - reverse)
// FF =^ RR

static const char* peOri[3]={"FR","RF","FF"};


extern struct hashTable *kmere;
extern struct readLink *links;
extern int graphSize;
extern char **readList;
extern int numreads;
extern volatile int maxReadId;

struct readFiles{
	char *leftReads;		// path to the original left read sequence file
	char *rightReads;		// path to the original right read sequence file, NULL if SE Lib
	int insertSize;			// if PE Lib -> fragment length (avg -> start read 1 to start read 2)
	int maxInsert;			// Maximum Insert size
	int minInsert;			// Minimum Insert size
	int oriPE;				// orientation of the PE reads (FR (Default), RF, FF) -> see macro definition
	int libNum;				// ID of Lib
	int startId;			// First read ID of this Lib
	int endId;				// Last read ID of this Lib
};

struct hash_block{
	int pthr_id;
	int pthr_num;
	long start;
	long end;
	long file_size;
	char* fasta;
};

void printUsage();

struct readFiles* readCMDline(int argc, char *argv[], struct readFiles *files);

int readFile(struct readFiles* files);

int readFastA(char* inFile, int readNum, int jump);

int readFastQ(char* inFile, int readNum, int jump);

void createKmers(char* read, int readNum);// ,struct hashTable *kmere);

void cleanGraph();

void iterLink();

void iterHash();

void inHash(char* buffer);

// MultiThread Reader

void* mt_fileReader(void* block);

void* mt_fileReaderDB(void* block);

void fileScheduler(char* inFile, int pthr_num, pthread_t* threads);

struct readFiles* fileScheduler_DB(char* dbFile, int pthr_num, pthread_t* threads);

//void hash_chain_len_histogram(UT_hash_table *tbl);


#endif /* FILEREADER_H_ */
