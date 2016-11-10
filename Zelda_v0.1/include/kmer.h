/*
 * kmer.h
 *
 *  Created on: Oct 13, 2014
 *      Author: kaempfpp
 */

#ifndef KMER_H_
#define KMER_H_

#include <stdint.h>
//#include "uthash.h"
#include <math.h>

#ifndef NULL
#define NULL   ((void *) 0)
#endif

#define _max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define _min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#define ABS(a)       ( ((a) < (0)) ? (-(a)) : (a) )
#define INITHASHSIZE(a) pow(2,a) + max_reprobes // initial 1 MBuckets / 20 bit

#define OB -1
#define AB 0
#define CB 1
#define GB 2
#define TB 3

#define TYPE128

//typedef atomic_uint_fast64_t KmerBitBuffer;
#ifdef TYPE128
	typedef __uint128_t KmerBitBuffer;
#else
	typedef uint64_t KmerBitBuffer;
#endif

typedef uint32_t readID;

extern int maxRlen;
extern int nK;
extern int dtSize;
extern int *readLenList;
extern int *readStartList;

extern volatile int max_reprobes; 	// max number of occupied table entries before resizing the hash table
KmerBitBuffer empty;					// Empty hash bucket

static const uint32_t EMPTYBUCKET = 0xFFFFFFFF;
static const readID  SET_READ_END =  0x80000000;  // 1000 0000 0000 0000
static const readID  DEL_READ_END =  0x7FFFFFFF;  // 0111 1111 1111 1111
static const char TRANS_MASK  = 0x03;
static const KmerBitBuffer  FLAGS_ENDPOS_DISABLE = 0x7FFFFFFF;  // 0111 1111 1111 1111 1111 1111 1111 1111

static const int codes[256] = {
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, AB, OB, CB, OB, OB, OB, GB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, TB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, AB, OB, CB, OB, OB, OB, GB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, TB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB,
  OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB, OB
};

static const int compcodes[52] = {
		'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A',
		'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
		'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A'
};

static const char rev_codes[6] = { 'A', 'C', 'G', 'T', 'N', '-'};

struct hashTable {
	KmerBitBuffer kmer;
	char trans;
	int index;  // If minus the direction has changed, take further the children instead of parents and vice versa
//	UT_hash_handle hhb;
};

struct readLink{
	KmerBitBuffer kmer;
	readID *read;
//	UT_hash_handle hhb;
};

struct kmerList{
	int num;
	struct hashTable *s;
};

KmerBitBuffer toBuffer(char*);

void toBufferPtr(char* read, KmerBitBuffer *kmer);

char* toSeq(KmerBitBuffer);

KmerBitBuffer revKmer(KmerBitBuffer);

char getTransBase(KmerBitBuffer *kmer, int index);

char getTransBaseDown(KmerBitBuffer *kmer, int dir);

struct hashTable* hasChild(KmerBitBuffer *kmer, char base, char* dir);

struct hashTable* hasParent(KmerBitBuffer *kmer, char base, char* dir);

struct hashkmer* hasChild_2(KmerBitBuffer *kmer, char base, char *dir);

struct hashkmer* hasParent_2(KmerBitBuffer *kmer, char base, char *dir);

uint32_t hasChild_2_oa(KmerBitBuffer*, char, char*);

uint32_t hasParent_2_oa(KmerBitBuffer*, char, char*);

#endif /* KMER_H_ */