/*
 ============================================================================
 Name        : FileReader.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Zelda File- and Parameter-Handling
 ============================================================================
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


extern struct hashTable *kmere;
extern struct readLink *links;
extern int graphSize;
extern char **readList;
extern int numreads;
extern volatile int maxReadId;

struct readFiles{
	char *leftReads;		// path to the original left read sequence file
	char *rightReads;		// path to the original right read sequence file, NULL if SE Lib
//	int insertSize;			// if PE Lib -> fragment length (avg -> start read 1 to start read 2)
	int maxInsert;			// Maximum Insert size
	int minInsert;			// Minimum Insert size
	int avgInsert;			// Average Insert size
	int oriPE;				// orientation of the PE reads (FR (Default), RF, FF) -> see macro definition
	int libNum;				// ID of Lib
	int startId;			// First read ID of this Lib
	int endId;				// Last read ID of this Lib
};

struct para{
	char* assemblyName;
	char* readDB;
	char* asemblyFolder;
	struct readFiles* files;
	unsigned int blocks;
	unsigned int threads;
	unsigned int kSize;
	unsigned int minOvlLen;
	char run;
	char prefilter;
};

struct hash_block{
	int pthr_id;
	int pthr_num;
	long start;
	long end;
	long file_size;
	char* fasta;
};

void finished(struct para* para);

struct para* readCMDline(int argc, char *argv[]);

struct readFiles* readCMDmakeDB(int argc, char *argv[], int libNum);

int readFile(struct readFiles* files);

int readFastA(char* inFile, int readNum, int jump);

int readFastQ(char* inFile, int readNum, int jump);

void cleanGraph();

void iterLink();

void iterHash();

void inHash(char* buffer);

// MultiThread Reader

void* mt_fileReader(void* block);

void* mt_fileReaderDB(void* block);

void fileScheduler(char* inFile, int pthr_num, pthread_t* threads);

struct readFiles* fileScheduler_DB(char* dbFile, int pthr_num, pthread_t* threads);

void freeFiles(struct para* para);

#endif /* FILEREADER_H_ */
