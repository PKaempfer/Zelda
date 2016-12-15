/*
 ============================================================================
 Name        : kmer.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : K-mer handling in character and binary based representation
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include "kmer.h"
#include "FileReader.h"
#include "kmerHash.h"
#include "DBGraph_oa.h"

volatile int max_reprobes = 50;
int maxRlen;
int nK = 47;
int dtSize;

int *readLenList;
int *readStartList;

KmerBitBuffer toBuffer(char* read){
	KmerBitBuffer kmer=0;
	int i;
	for (i = 0; i < nK; i++) {
		if(codes[(int)read[i]] == OB){
			printf("WARNING: Reading wrong base\n");
		}
		kmer |= ((KmerBitBuffer)codes[(int)read[i]] << ((nK-(i+1)))*2);
	}
	return kmer;
}

void toBufferPtr(char* read, KmerBitBuffer *kmer){
	static KmerBitBuffer test;
	test = 0;
	int i;
	for (i = 0; i < nK; i++) {
		test |= ((KmerBitBuffer)codes[(int)read[i]] << ((nK-(i+1))*2));
	}
	kmer=&test;
}

char* toSeq(KmerBitBuffer kmer){
	int i;
	static char* seq = NULL;
	if(!seq) seq=(char*)malloc(nK+1);
	for(i=0;i<nK;i++){
		seq[i] = rev_codes[ (kmer >> ((nK-(i+1))*2)) & 3];
	}
	seq[nK]='\0';
	return seq;
}

static KmerBitBuffer rev64(KmerBitBuffer kmer){
	kmer = ((kmer >> 2) & 0x3333333333333333)  | ((kmer << 2) & 0xcccccccccccccccc);
	kmer = ((kmer >> 4) & 0x0f0f0f0f0f0f0f0f)  | ((kmer << 4) & 0xf0f0f0f0f0f0f0f0);
	kmer = ((kmer >> 8) & 0x00ff00ff00ff00ff)  | ((kmer << 8) & 0xff00ff00ff00ff00);
	kmer = ((kmer >> 16) & 0x0000ffff0000ffff) | ((kmer << 16) & 0xffff0000ffff0000);
	kmer = ((kmer >> 32) & 0x00000000ffffffff) | ((kmer << 32) & 0xffffffff00000000);
	kmer=(~kmer)>>(64-(nK*2));
	return kmer;
}

static uint64_t rev64u(uint64_t kmer){
	kmer = ((kmer >> 2) & 0x3333333333333333)  | ((kmer << 2) & 0xcccccccccccccccc);
	kmer = ((kmer >> 4) & 0x0f0f0f0f0f0f0f0f)  | ((kmer << 4) & 0xf0f0f0f0f0f0f0f0);
	kmer = ((kmer >> 8) & 0x00ff00ff00ff00ff)  | ((kmer << 8) & 0xff00ff00ff00ff00);
	kmer = ((kmer >> 16) & 0x0000ffff0000ffff) | ((kmer << 16) & 0xffff0000ffff0000);
	kmer = ((kmer >> 32) & 0x00000000ffffffff) | ((kmer << 32) & 0xffffffff00000000);
	return kmer;
}

static KmerBitBuffer rev128(KmerBitBuffer kmer){
	uint64_t r1; // First 64 bit
	uint64_t r2; // Last 64 bit
	r1 = rev64u(kmer >> 64);
	r2 = rev64u(kmer);
	kmer = ((KmerBitBuffer)r2 << 64);
	kmer |= r1;
	kmer = ((~kmer) >> ((sizeof(KmerBitBuffer)*8)-(nK*2)));
	return kmer;
}

KmerBitBuffer revKmer(KmerBitBuffer kmer){
//	printf("dtSize: %i\n",dtSize);
	if(dtSize==64)	kmer = rev64(kmer);
	else kmer=rev128(kmer);
	return kmer;
}

char getTransBase(KmerBitBuffer *kmer, int dir){
	if(dir > 0) return rev_codes[((char)(*kmer)) & 0x03];
	else return rev_codes[( (char) ((~(*kmer)) >> ((nK-1)*2)) & 0x03)];
}

char getTransBaseDown(KmerBitBuffer *kmer, int dir){
	if(dir > 0) return rev_codes[(int)((*kmer) >> ((nK-1)*2))];
	else return rev_codes[((char)(~(*kmer))) & 0x03];
}

// can used only once for each vertex, else problems with indexing the nodes
struct hashkmer* hasChild_2(KmerBitBuffer *kmer, char base, char *dir){
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;
	temp = (*kmer) << 2;
	temp &= ~(((KmerBitBuffer)3) << (2*nK));
	temp |= (KmerBitBuffer)base;
	revtemp=revKmer(temp);
	if(revtemp<temp){
		temp = revtemp;
		(*dir) *= -1;
	}
	return getKmer(temp);
}

struct hashkmer* hasParent_2(KmerBitBuffer *kmer, char base, char *dir){
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;
	temp = (*kmer) >> 2;
	temp |= ((KmerBitBuffer)base << ((nK-1)*2));
	revtemp=revKmer(temp);
	if(revtemp<temp){
		temp = revtemp;
		(*dir) *= -1;
	}
	return getKmer(temp);
}

uint32_t hasChild_2_oa(KmerBitBuffer *kmer, char base, char *dir){
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;
	temp = (*kmer) << 2;
	temp &= ~(((KmerBitBuffer)3) << (2*nK));
	temp |= (KmerBitBuffer)base;
	revtemp=revKmer(temp);
	if(revtemp<temp){
		temp = revtemp;
		(*dir) *= -1;
	}
	return getKmer_oa(temp);
}

uint32_t hasParent_2_oa(KmerBitBuffer* kmer, char base, char* dir){
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;
	temp = (*kmer) >> 2;
	temp |= ((KmerBitBuffer)base << ((nK-1)*2));
	revtemp=revKmer(temp);
	if(revtemp<temp){
		temp = revtemp;
		(*dir) *= -1;
	}
	return getKmer_oa(temp);
}
