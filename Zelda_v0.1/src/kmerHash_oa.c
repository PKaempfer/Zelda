/*
 ============================================================================
 Name        : kmerHash_oa.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : This Data structure represents a open addressable hash table
 		 	   for concurrent, lock-free kmer hashing (Marçais and Kingsford,
 		 	   2011), exploiting the x86 architecture supported atomic CAS
 		 	   instruction
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#include "readDB.h"
#include "kmerHash_oa.h"
#include "kmerHash.h"
#include "FileReader.h"

volatile unsigned char resize_mutex = 0;
volatile unsigned char fin_mutex = 0;
unsigned char pthr_runN = 0;
// Set dependent on the given k size in mt_fileScheduler_DB
KmerBitBuffer NULL_KMER;
KmerBitBuffer NULL_KMER_LASTBITS;
KmerBitBuffer NULL_KMER_FIRSTBITS;

void createHashTable_oa(){
	printf("CHECKPOINT: Create HashTable\n");
	uint32_t i,j;
	NULL_KMER = 0;
	for(i=0;i<nK;i++){
		NULL_KMER = NULL_KMER << 2;
		NULL_KMER |= 3;
	}
	NULL_KMER_LASTBITS = NULL_KMER - 3;
	NULL_KMER_FIRSTBITS = 0;
	for(i=0;i<nK-1;i++){
		NULL_KMER_FIRSTBITS = NULL_KMER_FIRSTBITS << 2;
		NULL_KMER_FIRSTBITS |= 3;
	}

	empty = 0;
	for(i=0;i<sizeof(KmerBitBuffer)*4;i++){
		empty = empty << 2;
		empty |= 3;
	}

	j = (uint32_t)INITHASHSIZE(bitnum);
	dbHash_oa = (struct hashkmer_oa*)malloc(sizeof(struct hashkmer_oa)*j);
	for(i=0;i<j;i++){
		dbHash_oa[i].kmer = empty;
		dbHash_oa[i].count = 0;
		dbHash_oa[i].trans = 0;
		dbHash_oa[i].index = 0;
		dbHash_oa[i].ends = NULL;
	}
	bitmask = 0;
	for(i=0;i<bitnum;i++){
		bitmask = bitmask << 1;
		bitmask |= 1;
	}
	expansionThreshold = INITHASHSIZE(bitnum);
	printf("\tCreate HashTable of size: %i (expTreshold: %i)\n",j,expansionThreshold);
}

void freeHashTable_oa(){
	free((void*)dbHash_oa);
	dbHash_oa=NULL;
	bitnum -= 2;
}

void freeEnds_oa(){
	uint32_t j = INITHASHSIZE(bitnum);
	uint32_t i;
	volatile struct readEnd* readEnd;
	for(i=0;i<j;i++){
		readEnd = dbHash_oa[i].ends;
		while(readEnd){
			dbHash_oa[i].ends = readEnd->next;
			free((void*)readEnd);
			readEnd = dbHash_oa[i].ends;
		}
	}
}

volatile struct hashkmer_oa* resizeHashTable(){
	char verbose = 0;
	if(verbose) printf("Hash table full: RESIZE!\n");
	uint32_t i,j;
	volatile struct hashkmer_oa* tempHash_oa;
	if(verbose) hashStats_oa();
	bitnum++;
	j = (uint32_t)INITHASHSIZE(bitnum);
	if(verbose) printf("Set new bitnum: %i\n",bitnum);
	if(verbose) printf("New size of Hash table: %i\n",j);
	tempHash_oa = (struct hashkmer_oa*)malloc(sizeof(struct hashkmer_oa)*j);
	if(!tempHash_oa) return NULL;
	for(i=0;i<j;i++){
		tempHash_oa[i].kmer = empty;
		tempHash_oa[i].count = 0;
		tempHash_oa[i].trans = 0;
		tempHash_oa[i].index = 0;
		tempHash_oa[i].ends = NULL;
	}

	bitmask = 0;
	for(i=0;i<bitnum;i++){
		bitmask = bitmask << 1;
		bitmask |= 1;
	}

	j = (uint32_t)INITHASHSIZE((bitnum-1));
	uint32_t bucket;
	for(i=0;i<j;i++){
		if(dbHash_oa[i].count){
			bucket = my_hash(dbHash_oa[i].kmer);
			while(tempHash_oa[bucket].count)  bucket++;
			tempHash_oa[bucket] = dbHash_oa[i];
		}
	}
	free((void*)dbHash_oa);
	expansionThreshold = INITHASHSIZE(bitnum);
	if(verbose) printf("Create new HashTable of size: %i (expTreshold: %i)\n",j,expansionThreshold);
	return tempHash_oa;
}

char addReadEnd_oa(KmerBitBuffer current_new, readID read, char end, int32_t readlen){
	uint32_t bucket;
	struct readEnd* new_readEnd = (struct readEnd*)malloc(sizeof(struct readEnd));
	new_readEnd->next = NULL;
	new_readEnd->read = setID(read,end);
	if(end){
		readLenList[read] = readlen;
//		readStartList[read] = 0;
	}
//	new_readEnd->len = readlen;

	bucket = my_hash(current_new);
	int i=0;
	while(dbHash_oa[bucket].kmer != current_new){
		i++;
		bucket++;
	}

	struct readEnd** current_readEnd = (struct readEnd**)&dbHash_oa[bucket].ends;
	struct readEnd* temp;

	while((*current_readEnd)){
		current_readEnd = (struct readEnd**)&((*current_readEnd)->next);
	}

	do{
		if(*current_readEnd){
			current_readEnd = (struct readEnd**)&((*current_readEnd)->next);
		}
		temp = __sync_val_compare_and_swap(&(*current_readEnd),NULL,new_readEnd);
	} while(temp != NULL);

	return 1;
}

#ifdef TYPE128 // switch in kmer.h
// Support of Datatype depends on architecture -> Catch this cases (32 vs 64 bit, endianness)
// 128 bit version with linear probing (better balancing with double hashing, but higher number of cache misses)
// Two atomic 64 bit CAS instructions in a row guarantees data consistency
char addKmer128_oa(KmerBitBuffer current_new){
	uint32_t bucket;
	int i = 0;
	unsigned char counter;
	unsigned char new_counter;
	uint64_t current_key;

	if(resize_mutex){
		counter = resize_mutex;
		do{
			new_counter = counter;
			counter = __sync_val_compare_and_swap(&resize_mutex,new_counter,new_counter+1);
		} while(counter != new_counter);
//		printf("Thread Stop inserting, wait for free mutex -> oldHash: %p\n",dbHash_oa);
		// During waiting they could help rehashing:
		// -> Devide old hashTable into blocks:
		//    -> Each thread gets its individual part of the old hash table to rehash
		do{
			// Catch Thread in this while loop during resize!
		} while(resize_mutex);
//		printf("New HashTable found (%p), continue insert!\n",dbHash_oa);
	}

	bucket = my_hash(current_new)-1;
	do{
		do{
			if(i>max_reprobes){
				// There seems to be somehow a danger to end up in a deadlock in resize procedure!
				// Command all threads to Stop or return if another thread already is doing this
				if(!__sync_bool_compare_and_swap(&resize_mutex,0,1)){
//					printf("Already locked by other thread: Return and start new\n");
					return 0;
				}
				// Wait till all other threads report a stopped- or finished-state
				while(pthr_runN != resize_mutex + fin_mutex){
//					printf("pth: %i != %i + %i\n",pthr_runN, resize_mutex, fin_mutex);
//					sleep(1);
				}
//				printf("CHECKPOINT: Init new hash table (%i = %i + %i)\n",pthr_runN, resize_mutex, fin_mutex);
				// rehash all Keys -> do in parallel: rehash key -> set kmer -> memcpy
				volatile struct hashkmer_oa* temp = resizeHashTable();
//				printf("CHECKPOINT: Hash table finished\n");
				if(!temp){
					printf("Not able to reallocate memory for hash-table resize\nAbort!!!\n");
					exit(1);
				}
				dbHash_oa = temp;
//				hashStats_oa();
				resize_mutex = 0;
				i=0;
				bucket = my_hash(current_new)-1;
			}
			bucket++;
			current_key = __sync_val_compare_and_swap(&(((uint64_t*) &dbHash_oa[bucket].kmer)[1]),((uint64_t*)&empty)[1],((uint64_t*)&current_new)[1]);
			i++;
		} while(current_key != ((uint64_t*)&current_new)[1] && current_key != ((uint64_t*)&empty)[1]);
		current_key = __sync_val_compare_and_swap(&(((uint64_t*) &dbHash_oa[bucket].kmer)[0]),((uint64_t*)&empty)[0],((uint64_t*)&current_new)[0]);
	} while(current_key != ((uint64_t*)&current_new)[0] && current_key != ((uint64_t*)&empty)[0]);
//	dbHash_oa[bucket].count=1;
	counter = dbHash_oa[bucket].count;
	do{
		if(counter==255){
			break;
		}
		new_counter = counter;
		counter = __sync_val_compare_and_swap(&dbHash_oa[bucket].count,new_counter,new_counter+1);
	} while(counter != new_counter);
	return 1;
}

#else
// 64bit KmerBitBuffer data type
char addKmer_oa(KmerBitBuffer current_new){
	uint32_t bucket;
	int i = 0;
	int j;
	unsigned char counter;
	unsigned char new_counter;
	int reprobe_counter;
	int new_reprobe_counter;
	KmerBitBuffer current_key;

	// Find K-mer or set new if not existent
	// Atomic CAS with C11 standard Function in stdatomic.h (needs compilation with C11 standard)
	// 128bit dwCAS needs cmpxchg16b support by compiler and processor architecture
	// Does not work with current installed clang compiler version
//	do{
//		if(i > max_reprobe){
//			// means resize the hash table
//			// all threads have to stop for resizing
//			// Than double size the table and reorder all entries; time intensive (again parallelizable???)
//			printf("HashTable is full, resize!!!\n");
//			return FALSE;
//		}
//		bucket = my_hash(current_new) + i;
//		current_key = empty;
//		is_set = atomic_compare_exchange_weak(&dbHash_oa[bucket].kmer,&current_key,current_new);
//	} while(!is_set && current_key!=current_new);

	// Atomic CAS with C build_in Function
	do{

		if(i > max_reprobes){
			// means resize the hash table
			// all threads have to stop for resizing
			// Than double size the table and reorder all entries; time intensive (again parallelizable???)
			printf("Lock?: %i\n",lock.__align);
			pthread_mutex_lock(&lock);
			reprobe_counter = max_reprobes;
			do{
				new_reprobe_counter = reprobe_counter;
				reprobe_counter = __sync_val_compare_and_swap(&max_reprobes,new_reprobe_counter,new_reprobe_counter + 1);
			} while(reprobe_counter != new_reprobe_counter);
			j = INITHASHSIZE(bitnum);
			dbHash_oa[j-1].kmer = empty;
//			max_reprobes++;
			printf(" --> New max_reprobes number: %i\n", max_reprobes);
			sleep(1);
			printf("HashTable is full, resize!!!\n");
			pthread_mutex_unlock(&lock);
			return 0;
		}

		bucket = my_hash(current_new) + i;
//		current_key =  CAS_add(bucket,empty,current_new);
		current_key = __sync_val_compare_and_swap(&dbHash_oa[bucket].kmer,empty,current_new);
		i++;
	} while(current_key != current_new && current_key != empty);
	counter = dbHash_oa[bucket].count;
	do{
		if(counter==255){
//			printf("Counter Full\n");
			break;
		}
		new_counter = counter;
		counter = __sync_val_compare_and_swap(&dbHash_oa[bucket].count,new_counter,new_counter+1);
	} while(counter != new_counter);
	return 1;
}
#endif

uint32_t findKmer128_oa(KmerBitBuffer current_new){
	uint32_t bucket;

	bucket = my_hash(current_new);
//	printf("Bucket: %i\n",bucket);
	while(dbHash_oa[bucket].count && dbHash_oa[bucket].kmer != current_new){
		bucket++;
	}
	if(dbHash_oa[bucket].kmer == current_new) return bucket;

	return 0;
}


// Running Hash-Stats is essential to calc the total number of nodes in the following adjacency list
void hashStats_oa(){
	uint32_t dif_kmer = 0;
	uint64_t tot_num = 0;
	uint32_t j, i;
	j = INITHASHSIZE(bitnum);
	KmerBitBuffer old = empty;

	printf("\tNumber of Reads: %i\n",numreads);

	int unique = 0;
	int full = 0;
	int maxcount = 0;
	int distinct = 0;
	int totdist = 0;
	int setchains = 0;
	int usetchains = 0;
	int setchainlen = 0;
	int usetchainlen = 0;
	int chain=0;

	itemNum = 0;
	for(i=0;i<j;i++){
		if(dbHash_oa[i].kmer != empty){
			if(!chain){
				usetchains++;
				chain = 1;
			}
			setchainlen++;
			itemNum++;
			old = (KmerBitBuffer)dbHash_oa[i].kmer;
			dif_kmer++;
			tot_num += dbHash_oa[i].count;
			if(dbHash_oa[i].count==1) unique++;
			else{
				distinct ++;
				totdist += dbHash_oa[i].count;
			}
			if(dbHash_oa[i].count==255) full++;
			if(dbHash_oa[i].count > maxcount) maxcount = dbHash_oa[i].count;
			old = dbHash_oa[i].kmer;
		}
		else{
			if(chain){
				setchains++;
				chain = 0;
				usetchainlen++;
			}
			usetchainlen++;
		}
	}
	graphSize = itemNum+1;
	float coverage = (float)totdist/(float)distinct;
	float part = (float)maxReadLen/((float)(maxReadLen-nK)+1);
	float part2 = unique / (float)tot_num;
	printf("\tHashTable contains %i (%.2f%%) different k-mers with a total number of %lu\n",dif_kmer,(((float)dif_kmer/j))*100,tot_num);
	printf("\tUnique k-mer: (c=1)    %i\n",unique);
	printf("\tEstimated genome size: %i\n",distinct);
	printf("\tEstimated Coverage:    %.2f\n",(coverage * part) + (coverage * part * part2));
	printf("\tConunt >=255:          %i\n",full);
	printf("\tMaxCoverage:           %i\n",maxcount);
	printf("\tAL GraphSize:          %i\n",graphSize);
	printf("\tSetChains:             %i (avglen: %.2f)\n",setchains,(float)setchainlen/setchains);
}

void mt_createKmers(char* read, int readNum){
//	static clock_t times,timee=0,timeh=0;
//	if(readNum<16) printf("Read (%i): %s\n",readNum,read);
//	printf("Read (%i): %s\n",readNum,read);
	int len = strlen(read);
	int readLen;
	KmerBitBuffer temp;
	KmerBitBuffer revtemp;

	char* buffer;
	if(strchr(read,'0')){
		printf("Read contains Gaps\n");
		return;
	}

	readLen = strlen(read);
	if(readLen > maxRlen) maxRlen = len;
	buffer = read+(readLen-nK);
	len = readLen - nK;

	int testnum = 0;

	if(len<0){
		printf("return before ends were added (len: %i) -> read number inconsistency\n",len);
		return;
	}

	while(len!=-1){
		*(buffer+nK)='\0';
		temp = toBuffer(buffer);
//		printf("kmer: %s\n",toSeq(temp));
		revtemp = revKmer(temp);
		if(revtemp < temp){
			temp = revtemp;
		}
#ifdef TYPE128
		while(!addKmer128_oa(temp)){
			printf("Could not add Kmer, Try again\n");
		}
#else
		while(!addKmer_oa(temp)){
			printf("Could not add Kmer, Try again\n");
		}
#endif


		if(len == readLen - nK){
			addReadEnd_oa(temp,readNum,1,len);
			testnum++;
		}
		if(len == 0){
			addReadEnd_oa(temp,readNum,0,-1);
			testnum++;
		}
		len--;
		buffer--;
	}
	if(testnum != 2) printf("read number inconsistency\n");
}

/**
 * Function gets reads in binary representation and creats all its decompositions for insertion into the hash table.
 * DB-Support depends on architecture (Big-Endian vs. Little-Endian)
 * Catch the Case and build database regarding to the architecture
 * @param read		Represents the read in binary form
 * @param len		length of the read in decompressed format
 * @param readNum	ID of the read as in database
 */
void mt_createKmers_DB(char* readSeq, int len, int readNum){
	KmerBitBuffer kmercp;
	KmerBitBuffer temp;
	KmerBitBuffer revtemp;

	int readPos = 0;
	int bitpos = 0;
	int shift = 0;
	int shift2 = 0;
	char lastbyte;

	int byte = ((len+3)/4);

#ifdef TYPE128
	int cpByte = 16;
	if(byte < cpByte) cpByte = byte;
	if(len > maxRlen) maxRlen = len;
	byte-=(cpByte);
//	char* decompRead = decompressRead(readSeq,len);
//	printf("Read: %s\n",decompRead);
//	free(decompRead);

//	printf("readPos: %i, nK: %i, len: %i\n",readPos,nK,len);
	while(readPos + nK < len){
		temp = 0;
		if(bitpos==0){
//			printf("Copy bytes: %i (len: %i)\n",byte,cpByte);
			memcpy(&kmercp,&readSeq[byte],cpByte);
		}
		shift = ((8*cpByte) - (nK*2)) - (bitpos%8);
		if(shift>=0) temp = kmercp >> shift;
		else{
			shift2 = 8+shift;
			shift*=-1;
			lastbyte = readSeq[byte-1];
			temp = kmercp << shift;
			lastbyte = lastbyte >> shift2;
			lastbyte &= (char)pow(2,shift) - 1;
			temp |= lastbyte;
		}
		temp &= NULL_KMER;
		revtemp = revKmer(temp);
		if(revtemp < temp){
			temp = revtemp;
		}
		while(!addKmer128_oa(temp)){
			printf("Could not add Kmer, Try again\n");
		}
		if(readPos == 0){
			addReadEnd_oa(temp,readNum,0,-1);
		}
		bitpos += 2;
		readPos ++;
		if(bitpos==8){
			byte--;
			bitpos = 0;
		}
	}
	if(bitpos==0){
		memcpy(&kmercp,&readSeq[byte],cpByte);
	}
	shift = ((8*cpByte) - (nK*2)) - (bitpos%8);
	if(shift>=0) temp = kmercp >> shift;
	else{
		shift2 = 8+shift;
		shift*=-1;
		lastbyte = readSeq[byte-1];
		temp = kmercp << shift;
		lastbyte = lastbyte >> shift2;
		lastbyte &= (char)pow(2,shift) - 1;
		temp |= lastbyte;
	}
	temp &= NULL_KMER;
	revtemp = revKmer(temp);
	if(revtemp < temp){
		temp = revtemp;
	}
	while(!addKmer128_oa(temp)){
		printf("Could not add Kmer, Try again\n");
	}
	addReadEnd_oa(temp,readNum,1,len-nK);

#endif

}
