/*
 ============================================================================
 Name        : read_filter.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Prefilter to correct erroneous reads based on coverage drops
 	 	 	   in the initial build hash table.
 ============================================================================
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
void* mt_filter_reads_correction(void* filter_block){
	char verbose = 0;
	char verbose2 = 0;
	if(verbose) printf("Checkpoint: Create Mapping Thread\n");
	KmerBitBuffer kmercp;
	KmerBitBuffer tempCor;
	KmerBitBuffer tempCor2;
	KmerBitBuffer temp;
	KmerBitBuffer revtemp;
	KmerBitBuffer precurser = 0;
	uint32_t bucket;

	FILE* correctedR;
	if(verbose) correctedR = fopen("correctedReads","w");

	struct filter_block block = *((struct filter_block*)filter_block);
	struct reads* reads = block.reads;
	long i = block.start;
	long end = block.end;
//	if(block.pthr_id == 2) verbose = 1;

	int cov_tot;
	int cov_one;
	int cov_high;
	int cov;
	int pre_cov = 0;
	uint16_t best_len;
	uint16_t best_st;
	uint16_t best_end;
	int cutoff = 1;
	int j;

	int len;

	int readPos;
	int bitpos;
	int shift;
	int shift2;
	char lastbyte;
	int byte;
	int cpByte;
	char* decomRead = NULL;
	char* comRead = NULL;
	char* readSeq;
	char* readSeqC=(char*)malloc(150000);
	char* cutSeq=(char*)malloc(150000);
	int max_one;
	char* readSeqOrg = (char*)malloc((maxReadLen+3)/4);
	unsigned char* cArray = (unsigned char*)malloc(sizeof(unsigned char)*((maxReadLen-nK)+1));
	char redo = 0;
	char change_count = 0;

	if(verbose) printf("MaxReadLen: %i\n",maxReadLen);

	for(;i<=end;i++){
		if(redo == 0) change_count = 0;
		else change_count++;
		if(change_count > 3){
			redo = 0;
			continue;
		}
		if(verbose && (i-block.start)%100000==0) printf("Thread: %i: %li reads corrected\n",block.pthr_id,i-block.start);
		len = reads[i].len;
		if(len >= nK){
//			decomRead = decompressRead(reads[i].seq,len);
//			printf("read: %s\n",decomRead);
//			printf("Thread: %i with read: %i (len: %i nK: %i)\n",block.pthr_id,reads[i].ID, len, nK);
//			free(decomRead);
			cov_tot = 0;
			cov_one = 0;
			cov_high = 0;
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
				tempCor = temp;
				revtemp = revKmer(temp);
				if(revtemp < temp){
					temp = revtemp;
				}
				bucket = findKmer128_oa(temp);
				if(!bucket){
					if(temp==dbHash_oa[bucket].kmer && dbHash_oa[bucket].count){
						// Go ahead
					}
					else{
						if(verbose) printf("Read: %li -> Kmer is not existent. Undo Base substitution\n",i);
						decomRead = decompressRead(reads[i].seq,len);
						if(verbose) printf("old read: %s\n",decomRead);
						memcpy(reads[i].seq,readSeqOrg,(len+3)/4);
						if(redo == -2 && verbose){
							decomRead = decompressRead(reads[i].seq,len);
							printf("res read: %s\n",decomRead);
							exit(1);
						}
						redo = -1;
						i--;
						break;
					}

				}
				cov = dbHash_oa[bucket].count;
				cArray[readPos] = (unsigned char)cov;
				cov_tot += cov;
				if(cov<=cutoff) cov_one++;
				else cov_high++;
				bitpos += 2;

				// Read Correction
				if(readPos && redo != -2){
					// Case 1
					if(pre_cov < 3 && cov > 10){
						if(verbose2) printf("C1 -> ReadCorrection (Redo = %i)\n",redo);
						if(verbose) printf("END-Error - Read %li \n",i);
//						if(redo == 0) memcpy(readSeqOrg,reads[i].seq,(len+3)/4);
						char found = 0;
						precurser &= NULL_KMER_FIRSTBITS;
						char bestj = 0;
						unsigned char newCov = 0;
						for(char j=0;j<4;j++){
							KmerBitBuffer mask = (KmerBitBuffer)j << (nK-1)*2;
							tempCor2 = precurser | mask;

							revtemp = revKmer(tempCor2);
							if(revtemp < tempCor2){
								tempCor2 = revtemp;
							}
							bucket = findKmer128_oa(tempCor2);
							if(bucket && dbHash_oa[bucket].count > 5){
								newCov = dbHash_oa[bucket].count;
								found ++;
								bestj = j;
							}

							precurser &= NULL_KMER_FIRSTBITS;
						}
						if(found==1){
							if(verbose) printf("Found: %i -> Best: %c (old / new. %i/%i)\n",(int)found,rev_codes[(int)bestj],pre_cov,(int)newCov);
							if(redo == 0) memcpy(readSeqOrg,reads[i].seq,(len+3)/4);
//							memcpy(readSeqInter,reads[i].seq,(len+3/4));
							decomRead = decompressRead(reads[i].seq,len);
							if(verbose) printf("Old read: %s\n",decomRead);
							decomRead[readPos-1]=rev_codes[(int)bestj];
							if(verbose) printf("New read: %s\n",decomRead);
							comRead = compressRead(decomRead);
							memcpy(reads[i].seq,comRead,(len+3)/4);
							if(verbose2) printf("read: %s\n",decomRead);
							free(decomRead);
							free(comRead);
							redo = 1;
							if(verbose) printf("Thread: %i with read: %i (len: %i nK: %i)\n",block.pthr_id,reads[i].ID, len, nK);
							i--;
							break;

						}
						else{
							if(verbose) printf("No alternative base found\n");
						}
					}
					// Case 2
					else if(cov < 3 && pre_cov > 10){
						if(verbose2) printf("C2 -> ReadCorrection (Redo = %i)\n",redo);
						if(verbose) printf("BEGIN-Error - Read %li (redo=%i) pos: %i\n",i,(int)redo,readPos+(nK-1));
						char found = 0;
						char oldj = tempCor & 3;
						oldj = rev_codes[(int)oldj];
						tempCor &= NULL_KMER_LASTBITS;
						char bestj = 0;
						unsigned char newCov = 0;
						for(char j=0;j<4;j++){
							tempCor2 = tempCor;
							revtemp = revKmer(tempCor2);
							if(revtemp < tempCor2){
								tempCor2 = revtemp;
							}
							bucket = findKmer128_oa(tempCor2);
							if(bucket && dbHash_oa[bucket].count > 5){
								newCov = dbHash_oa[bucket].count;
								found ++;
								bestj = j;
							}
							tempCor ++;
						}
						if(found==1){
							if(verbose) printf("Found: %i (Base %c->%c (Pos: %i)) -> Best: %c (old / new. %i/%i)\n",(int)found,oldj,rev_codes[(int)bestj],readPos+(nK-1),rev_codes[(int)bestj],cov,(int)newCov);
							if(redo == 0) memcpy(readSeqOrg,reads[i].seq,(len+3)/4);
//							memcpy(readSeqInter,reads[i].seq,(len+3/4));
							decomRead = decompressRead(reads[i].seq,len);
							if(verbose) printf("Old read: %s\n",decomRead);
							decomRead[readPos+(nK-1)]=rev_codes[(int)bestj];
							if(verbose) printf("New read: %s\n",decomRead);
							comRead = compressRead(decomRead);
							memcpy(reads[i].seq,comRead,(len+3)/4);
							if(verbose2) printf("read: %s\n",decomRead);
							free(decomRead);
							free(comRead);
							redo = 1;
							if(verbose) printf("Thread: %i with read: %i (len: %i nK: %i)\n",block.pthr_id,reads[i].ID, len, nK);
							i--;
							break;
						}
						else{
							if(verbose) printf("No alternative base found\n");
						}

					}
				}

				pre_cov = cov;
				precurser = tempCor;

				readPos ++;
				if(bitpos==8){
					byte--;
					bitpos = 0;
				}
			}
			if(redo == 1 || redo == -1){
				if(redo == -1) redo = -2;
				else redo = 2;
				continue;
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
			if(cov<=cutoff) cov_one++;
			else cov_high++;
			cArray[readPos] = (unsigned char)cov;
//			printf("Number of error Kmers: %i\n",cov_one);
			max_one =_min(nK,((len-nK)+1));
			if(cov_one){
				if(cov_high > cov_one*3){
					redo = 0;
					continue;
					best_len = 0;
					best_st = 0;
					best_end = 0;
					for(j=0;j<=readPos;j++){
						if(cArray[j]>cutoff) best_len++;
						else{
							if(best_len > (best_end - best_st) + 1){
								best_end = j-1;
								best_st = best_end - (best_len - 1);
							}
							best_len = 0;
						}
						if(verbose) printf("%i  ",cArray[j]);
					}
					if(verbose)printf("\n");
					if(best_len > (best_end - best_st) + 1){
						best_end = j-1;
						best_st = best_end - (best_len - 1);
					}
					best_len = (best_end - best_st) + (nK);
					decomRead = decompressRead(reads[i].seq,len);
					if(verbose) printf("Cut from %i to %i\n",best_st,best_end+(nK-1));
					if(verbose) printf("Old Read: %s\n",decomRead);
					memcpy(cutSeq,&decomRead[best_st],best_len);
					cutSeq[best_len]='\0';
					if(verbose)printf("New Read: %s\n",cutSeq);
					comRead = compressRead(cutSeq);
					if(verbose) printf("Read %li -> Cut: %i bp -> %i bp\n",i,len,best_len);
					if(verbose) fprintf(correctedR,">%li\n",i);
					if(verbose) fprintf(correctedR,"%s\n",cutSeq);
					reads[i].len = best_len;
					memcpy(reads[i].seq,comRead,(best_len+3)/4);
					free(comRead);
					free(decomRead);
				}
				else if(!cov_high){
					if(verbose) printf("Read %li -> Deleted\n",i);
					reads[i].len = 0;
					free(reads[i].seq);
					reads[i].seq = NULL;
				}
			}

//			if(cov_one >= max_one-15){
//
//				if(verbose) printf("Read %li -> Deleted\n",i);
//				reads[i].len = 0;
//				free(reads[i].seq);
//				reads[i].seq = NULL;
//			}
			if(redo==2 && verbose2) printf("Read %li -> Corrected\n",i);
			redo = 0;
		}
	}
	if(verbose2) printf("Thread %i finished correction\n",block.pthr_id);
	free(readSeqC);
	free(cArray);
	if(verbose) fclose(correctedR);
	return NULL;
}

//void mt_correct_reads(void* filter_block){
//
//}

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
//	char* decomRead;
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

void filter_reads(struct reads* reads, const int pthr_num, pthread_t* threads){
	printf("CHECKPOINT: Filter Reads\n");
	struct filter_block* filter_block = (struct filter_block*)malloc(sizeof(struct filter_block)*pthr_num);

	int i;
	int part = numreads / pthr_num;

	long start = 0;
	for(i=0;i<pthr_num;i++){
//		printf("Start with thread: %i\n",i);
		filter_block[i].start = start;
		filter_block[i].end = (i+1)*part;
		start = filter_block[i].end + 1;
		filter_block[i].pthr_id = i;
		filter_block[i].pthr_num = pthr_num;
		filter_block[i].reads = reads;
		pthread_create(&threads[i],NULL,mt_filter_reads_correction,(void*)&filter_block[i]);
//		pthread_create(&threads[i],NULL,mt_filter_reads,(void*)&filter_block[i]);
	}

	void* status;
	for(i=0;i<pthr_num;i++){
		pthread_join(threads[i],&status);
	}
	free(filter_block);
}

