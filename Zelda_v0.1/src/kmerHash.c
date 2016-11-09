/*
 * kmerHash.c
 *
 *  Created on: Mar 2, 2015
 *      Author: kaempfpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kmer.h"
#include "kmerHash.h"
#include "FileReader.h"
#include <stdint.h>
#include "math.h"

void **dbHashg = NULL;
int bitnum = 25; // Initial hash size
uint32_t itemNum = 0;
uint32_t expansionThreshold = 0;
uint32_t bitmask = 0x0;

void createHashTable(){
	uint32_t i,j;
	j = (uint32_t)INITHASHSIZE(bitnum);
	dbHashg = (void**)malloc(sizeof(void*)*j);
	for(i=0;i<j;i++){
		dbHashg[i] = NULL;
	}
	for(i=0;i<bitnum;i++){
		bitmask = bitmask << 1;
		bitmask |= 1;
	}
	expansionThreshold = INITHASHSIZE(bitnum);
	printf("Create HashTable of size: %i (expTreshold: %i)\n",j,expansionThreshold);
}

char insertKmer(KmerBitBuffer *kmer){
	static uint32_t bucket;
	static struct hashkmer* oldkmer;
	KmerBitBuffer realkmer = (*kmer);
	bucket = my_hash(realkmer);
	if(dbHashg[bucket]){
		oldkmer = dbHashg[bucket];
		while(oldkmer){
			if(oldkmer->kmer == realkmer){
				if(oldkmer->count<255) oldkmer->count++;
				return 0;
			}
			oldkmer = oldkmer->next;
		}
		struct hashkmer* newkmer = (struct hashkmer*)malloc(sizeof(struct hashkmer));
		newkmer->kmer = realkmer;
		newkmer->index = 0;
		newkmer->trans = 0;
		newkmer->count = 1;
		newkmer->next = dbHashg[bucket];
		dbHashg[bucket] = newkmer;
		itemNum++;
		if(itemNum > expansionThreshold){
			reSize();
			expansionThreshold *= 2;
		}
		return 1;
	}
	else{
		struct hashkmer* newkmer = (struct hashkmer*)malloc(sizeof(struct hashkmer));
		newkmer->kmer = realkmer;
		newkmer->index = 0;
		newkmer->trans = 0;
		newkmer->next = NULL;
		newkmer->count = 1;
		dbHashg[bucket] = newkmer;
		itemNum++;
		return 1;
	}
	return 0;
}

struct hashkmer* getKmer(KmerBitBuffer kmer){
	static uint32_t bucket;
	static struct hashkmer *kmer_ptr;
	bucket = my_hash(kmer);
	if((kmer_ptr = dbHashg[bucket])){
		do {
			if(kmer_ptr->kmer == kmer) return kmer_ptr;
			kmer_ptr = kmer_ptr->next;
		} while(kmer_ptr);
	}
	return NULL;
}


void deleteKmer(){

}

void cleanTable(){
	printf("Clean HashTable\n");

	struct hashkmer* kmer;
	uint32_t i;
	uint32_t j = INITHASHSIZE(bitnum);

	for(i=0; i < j; i++) {
		while(dbHashg[i]){
			kmer = dbHashg[i];
			kmer = kmer->next;
			free(dbHashg[i]);
			dbHashg[i] = kmer;
		}
	}
	free(dbHashg);
}

/**
 * Table becomes to inefficient due to too many collisions and has to be resized
 * --> increase bitnum / add bit in bitmask
 * --> Hash function provides only 32 bit hash values, table can't be expanded if bitnum reaches 32
 * --> Double size of the Table
 * --> NULL new elements
 * --> Rehash all old elements in the table by same hash function but using 1 more bit as return value (half of entries do not have to be touched)
  */
char reSize(){

	bitnum++;
	bitmask = bitmask << 1;
	bitmask |= 1;

	if(bitnum > 32){
		printf("Can't reallocate memory: Max hash size reached (uint32_t)!\n");
		return 0;
	}

	void **temp;
	uint32_t j, i;
	j = INITHASHSIZE(bitnum);
	temp = realloc(dbHashg,sizeof(void*)*j);
	if(temp){
//		fflush(stdout);
//		printf("\n");
//		printf("Expand HashTable !!! \n");
		dbHashg = temp;
	}
	else{
		printf("Can't reallocate memory: Exit program!\n");
		free(dbHashg);
		exit(1);
	}
	for(i = INITHASHSIZE((bitnum-1));i<j;i++){
		dbHashg[i] = NULL;
	}

	j = INITHASHSIZE((bitnum-1));
	struct hashkmer *newkmer, *oldkmer, *tempkmer;
	uint32_t bucket;
	for(i=0;i<j;i++){
		while(dbHashg[i]){
			oldkmer = dbHashg[i];
			bucket = my_hash(oldkmer->kmer);
			if(bucket!=i){
				newkmer = oldkmer->next;
				oldkmer->next = dbHashg[bucket];
				dbHashg[bucket] = oldkmer;
				dbHashg[i] = newkmer;
			}
			else{
				newkmer = oldkmer->next;
				while(newkmer){
					bucket = my_hash(newkmer->kmer);
					if(bucket!=i){
						tempkmer = newkmer->next;
						newkmer->next = dbHashg[bucket];
						dbHashg[bucket] = newkmer;
						oldkmer->next = tempkmer;
						newkmer = tempkmer;
					}
					else{
						oldkmer = newkmer;
						newkmer = newkmer->next;
					}
				}
				break;
			}
		}
	}

	return 1;
}

readID setID(int readNum, int end){
	readID rtemp;
	rtemp=(readID)readNum;
	if(end) rtemp |= SET_READ_END;
	return rtemp;
}

void hashStats(){
	printf("Histogram\n");
	unsigned i, bkt_hist[CHAIN_MAX+1];
	uint32_t j = INITHASHSIZE(bitnum);
	int count, allcount = 0;
	double pct = 100.0/j;
	memset(bkt_hist,0,sizeof(bkt_hist));
	struct hashkmer* kmer;
	for(i=0; i < j; i++) {
		count=0;
		kmer = dbHashg[i];
		while(kmer){
			count++;
			kmer = kmer->next;
		}
		allcount += count;
		if (count == 0) bkt_hist[CHAIN_0]++;
		else if (count < 5) bkt_hist[CHAIN_5]++;
		else if (count < 10) bkt_hist[CHAIN_10]++;
		else if (count < 20) bkt_hist[CHAIN_20]++;
		else if (count < 100) bkt_hist[CHAIN_100]++;
		else bkt_hist[CHAIN_MAX]++;
	}
	fprintf(stderr, "Buckets with 0 items: %.1f%%\n", bkt_hist[CHAIN_0 ]*pct);
	fprintf(stderr, "Buckets with < 5 items: %.1f%%\n", bkt_hist[CHAIN_5 ]*pct);
	fprintf(stderr, "Buckets with < 10 items: %.1f%%\n", bkt_hist[CHAIN_10]*pct);
	fprintf(stderr, "Buckets with < 20 items: %.1f%%\n", bkt_hist[CHAIN_20]*pct);
	fprintf(stderr, "Buckets with < 100 items: %.1f%%\n", bkt_hist[CHAIN_100]*pct);
	fprintf(stderr, "Buckets with > 100 items: %.1f%%\n", bkt_hist[CHAIN_MAX]*pct);
	printf("Number of all Item: %i in %i buckets\n",allcount,j);


}


/**
 * my_hash use Jenkins hash algorithm
 * @param kmer is kmer sequence as bit representation
 * @return hash value
 */

uint32_t my_hash(KmerBitBuffer kmer)
{
	uint32_t hash, i;

    for(hash = i = 0; i < (uint32_t)(nK/4)+1; ++i){
        hash += (char) kmer;
        kmer = kmer >> 8;
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return (hash & bitmask);
}




uint32_t murmur3_32(const char *key) {
	static const uint32_t c1 = 0xcc9e2d51;
	static const uint32_t c2 = 0x1b873593;
	static const uint32_t r1 = 15;
	static const uint32_t r2 = 13;
	static const uint32_t m = 5;
	static const uint32_t n = 0xe6546b64;

	uint32_t hash = 0;
	uint32_t len = (nK/4)+1;

	const int nblocks = len / 4;
	const uint32_t *blocks = (const uint32_t *) key;
	int i;
	for (i = 0; i < nblocks; i++) {
		uint32_t k = blocks[i];
		k *= c1;
		k = (k << r1) | (k >> (32 - r1));
		k *= c2;

		hash ^= k;
		hash = ((hash << r2) | (hash >> (32 - r2))) * m + n;
	}

	const uint8_t *tail = (const uint8_t *) (key + nblocks * 4);
	uint32_t k1 = 0;

	switch (len & 3) {
	case 3:
		k1 ^= tail[2] << 16;
	case 2:
		k1 ^= tail[1] << 8;
	case 1:
		k1 ^= tail[0];

		k1 *= c1;
		k1 = (k1 << r1) | (k1 >> (32 - r1));
		k1 *= c2;
		hash ^= k1;
	}

	hash ^= len;
	hash ^= (hash >> 16);
	hash *= 0x85ebca6b;
	hash ^= (hash >> 13);
	hash *= 0xc2b2ae35;
	hash ^= (hash >> 16);

	return (hash & bitmask);
}

uint32_t jenkins_one_at_a_time_hash(char *key, size_t len)
{
    uint32_t hash, i;
    for(hash = i = 0; i < len; ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}


uint32_t iter_pos = 0;
uint32_t iter_j = 0;
struct hashkmer *iter_kmer;

void setIter(){
	iter_pos = 0;
	iter_kmer = NULL;
	iter_j = INITHASHSIZE(bitnum);
}

struct hashkmer* iterKmer(){
	if(iter_kmer == NULL){
		// search first element
		printf("Search first element\n");
		while(!dbHashg[iter_pos]){
			iter_pos++;
		}
		iter_kmer = dbHashg[iter_pos];
		return iter_kmer;
	}

	iter_kmer = iter_kmer->next;
	if(iter_kmer) return iter_kmer;
	else if(iter_pos+1<iter_j){
		iter_pos++;
		while(!dbHashg[iter_pos]){
			iter_pos++;
			if(iter_pos>=iter_j) return NULL;
		}
		iter_kmer = dbHashg[iter_pos];
		return iter_kmer;
	}
	else return NULL;
}



