/*
 * kmerHash.h
 *
 *  Created on: Mar 2, 2015
 *      Author: kaempfpp
 */

#ifndef KMERHASH_H_
#define KMERHASH_H_

volatile struct hashkmer_oa* dbHash_oa;

extern void **dbHashg;
extern uint32_t itemNum;
extern int bitnum;
extern uint32_t bitmask;
extern uint32_t expansionThreshold;

struct hashkmer {
	KmerBitBuffer kmer;
	char trans;
	unsigned char count;
	int index;  // If minus the direction has changed, take further the children instead of parents and vice versa
	struct hashkmer *next;
};

//struct hashTable{
//	KmerBitBuffer kmer;
//};

readID setID(int, int);

void createHashTable();

char insertKmer(KmerBitBuffer *kmer);

struct hashkmer* getKmer(KmerBitBuffer kmer);

void setIter();

struct hashkmer* iterKmer();

char reSize();

void hashStats();

uint32_t my_hash(KmerBitBuffer key);


#endif /* KMERHASH_H_ */
