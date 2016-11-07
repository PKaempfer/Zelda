/*
 * kmerHash_oa.h
 *
 *  Created on: Jan 25, 2016
 *      Author: lkaempfpp
 */

#ifndef KMERHASH_OA_H_
#define KMERHASH_OA_H_

#include "kmer.h"

extern volatile unsigned char resize_mutex;
extern volatile unsigned char fin_mutex;
extern unsigned char pthr_runN;

struct readEnd{
	int32_t len;
	readID read;
	volatile struct readEnd* next;
};

struct hashkmer_oa{
	volatile KmerBitBuffer kmer;
	unsigned char trans;
	volatile unsigned char count;
	volatile struct readEnd* ends;
	int32_t index;  // If minus the direction has changed, take further the children instead of parents and vice versa
};

void createHashTable_oa();

void freeHashTable_oa();

char addReadEnd_oa(KmerBitBuffer, readID, int, int32_t);

#ifdef TYPE128
char addKmer128_oa(KmerBitBuffer); // 128 bit (dwCAS)
#else
char addKmer_oa(KmerBitBuffer); // 64 bit (CAS)
#endif

void hashStats_oa();

// MultiThread DBG

void mt_createKmers(char* read, int readNum);

void mt_createKmers_DB(char* read, int len, int readNum);

#endif /* KMERHASH_OA_H_ */
