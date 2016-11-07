/*
 * Misc_hash.c
 *
 *  Created on: Aug 11, 2016
 *      Author: kaempfpp
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "../include/misc_hash.h"

#define CHAIN_0 0
#define CHAIN_5 1
#define CHAIN_10 2
#define CHAIN_20 3
#define CHAIN_100 4
#define CHAIN_MAX 5

#define INITHASHSIZE(a) pow(2,a)

struct hashTable32* hashTable32_create(int size){
	printf("Create HashTable of size: %.0f (in bits: %i)\n",INITHASHSIZE(size),size);
	uint32_t i;
	struct hashTable32* ht = (struct hashTable32*)malloc(sizeof(struct hashTable32));
	ht->elements = 0;
	ht->size = size; // Size in Bits
	ht->bitmask = 0x0;
	ht->maxrp = 5;
	ht->bucket = (struct bucket**)malloc(sizeof(struct bucket*)*INITHASHSIZE(size));

	if(!ht->bucket){
		printf("Could not allocate the requested number of buckets. ABORT and EXIT!\n");
		exit(1);
	}

	for(i=0;i<size;i++){
		ht->bitmask = ht->bitmask << 1;
		ht->bitmask |= 1;
	}

	for(i=0;i<INITHASHSIZE(ht->size);i++){
		ht->bucket[i] = NULL;
	}

	return ht;
}

void hashTable32_resize(struct hashTable32* ht){
	printf("Resize HashTable\n");
	if(ht->size>=32){
		printf("HashTable reached its MaxSize already\n");
		printf("No Reallocation. Resize maximum bucket size instead\n");
		ht->maxrp *= 2;
		return;
	}
	ht->size++;
	ht->bitmask = ht->bitmask << 1;
	ht->bitmask |= 1;
	struct bucket** htb = (struct bucket**)realloc(ht->bucket,sizeof(struct bucket*)*INITHASHSIZE(ht->size));
	if(htb) ht->bucket = htb;
	else{
		printf("Could not reallocate the requested number of buckets. ABORT and EXIT!\n");
		exit(1);
	}

	printf("New hashSize: %i (in bits: %i)",(int)INITHASHSIZE(ht->size),(int)ht->size);
	uint32_t i;
	for(i=INITHASHSIZE((ht->size-1));i<INITHASHSIZE(ht->size);i++){
		ht->bucket[i] = NULL;
	}
	struct bucket* bucket_o;
	struct bucket* bucket_n;
	struct bucket* bucket_t;
	uint32_t bucket;
	for(i=0;i<=INITHASHSIZE((ht->size-1));i++){
		while(htb[i]){
			bucket_o = htb[i];
			bucket = misc_myhash(bucket_o->key) & ht->bitmask;
//			printf("Insert at backut: %u\n",bucket);
			if(bucket != i){
				bucket_n = bucket_o->next;
				bucket_o->next = htb[bucket];
				htb[bucket] = bucket_o;
				htb[i] = bucket_n;
			}
			else{
				bucket_n = bucket_o->next;
				while(bucket_n){
					bucket = misc_myhash(bucket_n->key) & ht->bitmask;
					if(bucket!=i){
						bucket_t = bucket_n->next;
						bucket_n->next = htb[bucket];
						htb[bucket] = bucket_n;
						bucket_o->next = bucket_t;
						bucket_n = bucket_t;
					}
					else{
						bucket_o = bucket_n;
						bucket_n = bucket_n->next;
					}
				}
				break;
			}
		}
	}
}

// clean
void hashTable32_clean(struct hashTable32* ht){
	uint32_t i;
	struct bucket** bucket = ht->bucket;
	struct bucket* bucket_t;
	for(i=0;i<INITHASHSIZE(ht->size);i++){
		if(bucket[i]){
			while(bucket[i]->next){
				bucket_t = bucket[i]->next->next;
				free(bucket[i]->next);
				bucket[i]->next = bucket_t;
			}
		}
		free(bucket[i]);
	}
	free(bucket);
	free(ht);
}

void hashTable32_stats(struct hashTable32* ht){
	printf("HashTableStats:\n");
	printf("Histogram\n");
	unsigned i, bkt_hist[CHAIN_MAX+1];
	uint32_t j = INITHASHSIZE(ht->size);
	int count, allcount = 0;
	double pct = 100.0/j;
	memset(bkt_hist,0,sizeof(bkt_hist));
	struct bucket* bucket;
	for(i=0; i < j; i++) {
		count=0;
		bucket = ht->bucket[i];
		while(bucket){
			count++;
			bucket = bucket->next;
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

// insert
void hashTale32_insert(struct hashTable32* ht, uint32_t key, uint32_t value){
	uint32_t bucket = misc_myhash(key) & ht->bitmask;
	uint8_t  elem = 0;
	if(!ht->bucket[bucket]){
		ht->bucket[bucket] = (struct bucket*)malloc(sizeof(struct bucket));
		ht->bucket[bucket]->key = key;
		ht->bucket[bucket]->value = value;
		ht->bucket[bucket]->next = NULL;
		ht->elements++;
		return;
	}
	else{
		elem++;
		struct bucket* bucket_n = ht->bucket[bucket];
		while(bucket_n){
			if(bucket_n->key == key) return;
			else if(!bucket_n->next){
				bucket_n->next = (struct bucket*)malloc(sizeof(struct bucket));
				bucket_n->next->key = key;
				bucket_n->next->value = value;
				bucket_n->next->next = NULL;
				ht->elements++;
				elem++;
				if(elem>ht->maxrp) hashTable32_resize(ht);
				return;
			}
			elem++;
			bucket_n = bucket_n->next;
		}
	}
}

// find -> return value 0 means element not found. Do NOT use if 0 is a proper value!!!
uint32_t hashTable32_find(struct hashTable32* ht, uint32_t key){
	uint32_t value = 0;
	uint32_t bucket = misc_myhash(key) & ht->bitmask;
	struct bucket* bucket_n = ht->bucket[bucket];
	while(bucket_n){
		if(bucket_n->key == key) return bucket_n->value;
		bucket_n = bucket_n->next;
	}
	return value;
}

// delete

// Derived from Jenkins's (one-at-a-time) HashFunction
uint32_t misc_myhash(uint32_t value){
	uint32_t hash, i;

    for(hash = i = 0; i < sizeof(uint32_t); ++i){
        hash += (char) value;
        value = value >> 8;
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return (hash);
}
