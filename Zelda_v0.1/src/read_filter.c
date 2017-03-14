/*
 * read_filter.c
 *
 *  Created on: Feb 22, 2017
 *      Author: lkaempfpp
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "read_filter.h"
#include "kmer.h"
#include "kmerHash.h"
#include "kmerHash_oa.h"

/**
 * Multi threaded reads filter function. Set length of erroneous reads to zero.
 * Is not used in next hash table construction.
 */
void* mt_filter_reads(void* filter_block){
	printf("Checkpoint: Create Mapping Thread\n");
	KmerBitBuffer kmercp;
	KmerBitBuffer temp;
	KmerBitBuffer revtemp;
	uint32_t bucket;

	struct filter_block block = *((struct filter_block*)filter_block);
	struct reads* reads = block.reads;
	long i = block.start;
	long end = block.end;

	int cov_tot;
	int cov_one;
	int cov;

	int len;

	int readPos;
	int bitpos;
	int shift;
	int shift2;
	char lastbyte;
	int byte;
	int cpByte;
	char* decomRead;
	char* readSeq;
	char* readSeqC=(char*)malloc(150000);
	int max_one;

	for(;i<=end;i++){
		len = reads[i].len;
		if(len >= nK){
//			decomRead = decompressRead(reads[i].seq,len);
//			printf("read: %s\n",decomRead);
//			printf("Thread: %i with read: %i (len: %i nK: %i)\n",block.pthr_id,reads[i].ID, len, nK);
//			free(decomRead);
			cov_tot = 0;
			cov_one = 0;
			readPos = 0;
			bitpos = 0;
			shift = 0;
			shift2 = 0;
			byte = ((len+3)/4);
			cpByte = 16;
			if(byte < cpByte) cpByte = byte;
//			printf("Init bytes: %i (len: %i)\n",byte,cpByte);
			byte-=(cpByte);
			memcpy(&readSeqC[16],reads[i].seq,(len+3)/4);
			readSeq=&readSeqC[16];
			while(readPos + nK < len){
				temp = 0;
				if(bitpos==0){
//					printf("Copy bytes: %i (len: %i)\n",byte,cpByte);
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
				bucket = findKmer128_oa(temp);
				cov = dbHash_oa[bucket].count;
				cov_tot += cov;
				if(cov==1) cov_one++;
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
			bucket = findKmer128_oa(temp);
			cov = dbHash_oa[bucket].count;
//			printf("Number of error Kmers: %i\n",cov_one);
			max_one =_min(nK,((len-nK)+1));
			if(cov_one >= max_one-15){
				reads[i].len = 0;
				free(reads[i].seq);
				reads[i].seq = NULL;
			}
		}
	}
	free(readSeqC);
	return NULL;
}

//void mt_correct_reads(void* filter_block){
//
//}

void filter_reads(struct reads* reads, const int pthr_num, pthread_t* threads){
	printf("TEst\n");
	printf("CHECKPOINT: Filter Reads\n");
	printf("TreadNumber: %i\n",pthr_num);
	struct filter_block* filter_block = (struct filter_block*)malloc(sizeof(struct filter_block)*pthr_num);
	printf("NumReads: %i\n",numreads);
	printf("TreadNumber: %i\n",pthr_num);

	int i;
	int part = numreads / pthr_num;

	long start = 0;
	printf("TreadNumber: %i\n",pthr_num);
	for(i=0;i<pthr_num;i++){
		printf("Start with thread: %i\n",i);
		filter_block[i].start = start;
		filter_block[i].end = (i+1)*part;
		start = filter_block[i].end + 1;
		filter_block[i].pthr_id = i;
		filter_block[i].pthr_num = pthr_num;
		filter_block[i].reads = reads;
		pthread_create(&threads[i],NULL,mt_filter_reads,(void*)&filter_block[i]);
	}

	void* status;
	for(i=0;i<pthr_num;i++){
		pthread_join(threads[i],&status);
	}
	free(filter_block);
}

