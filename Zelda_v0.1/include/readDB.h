/*
 ============================================================================
 Name        : readDB.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Zelda Database Handling
 ============================================================================
 */

#ifndef READDB_H_
#define READDB_H_

#include "FileReader.h"

static const char version[] = "Zelda_v0.1";
extern int maxReadLen;

/** FLAG meanings, important for Scaffolding */
#define SINGLE 		0	// -> All SE Reads get this flag if still part of the ovl-graph
#define WIDOW		1	// PE-Read non existent anymore, read became a widow
#define CONCORDANT 2	// PE-Read exists and is concordant with the given range around the insert size (If concordant, then both reads)
#define DISCORDANT 3	// PE-Read exists but is either in the wrong orientation or discordant with the given range around the insert size (if discordant, then both reads)
#define REJECTED	4	// Read is out of observation for some reason

struct reads{
	int ID;				// ID - not necessary but should be written in binary outFile
	int len;			// length of the read in bits (IMPORTANT not a byte number)
	char flag;			// Read annotation -> see macro definitions above
	char* seq;			// compressed read sequence in bit-representation (00->A,01->C,10->G;11->T) // Don't store reads containing gaps
	void* annotation;	// Void because can be either be of type j_anno* if read describes a JUNCTION in the graph or pc_anno* if it is PROPER or CONTAINED -> See Definition of struct j_anno, struct pc_anno
};

void runTest();

int makeDB(char* outDB, int blocks, struct readFiles* files);

void writeDB(char* outDB, int blocks, struct readFiles* files);

void write_filteredFasta(struct readFiles* files, struct reads* reads);

void write_filteredDB(char* outDB, int blocks, struct readFiles* files, struct reads* reads);

struct reads* readDB(char* outDB);

void freeDB(struct reads* reads);

int readFastA_DB(char* inFile, int readNum, int jump);

int readFastQ_DB(char* inFile, int readNum, int jump);

char* decompressRead(char* compRead, int len);

void decompressReadSt(char* compRead, char* read, int len);

char* compressRead(char* compRead);


#endif /* READDB_H_ */
