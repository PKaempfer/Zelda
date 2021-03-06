/*
 ============================================================================
 Name        : readDB.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Zelda Database Handling
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "readDB.h"
#include "FileReader.h"
#include "kmer.h"

/**
 * Defines size of the stream buffer size of reading blocks: Default setting of 4mb (may depend on architecture)
 */
#define BUFFER_SIZE (4 * 1024 * 1024)

int readTotNum = 0;
int maxreadTotNum = 16 * 1024;
int maxReadLen = 0;
struct reads* readsList;

static const char* peOri[3]={"FR","RF","FF"};

/**
 * MakeDB creates a database of all libraries under consideration in a single database. For this the reads are stored in DB in their final Bit-representation.
 * This may be critical because the interpretation of the bit order may depend on the architecture the tool is running on (little-endian vs. big-endian).
 * All functions depending on the database representation are only supported by little-endian machines.
 * The DB stores additionally meta information about the considered input as: library-names, library-type, insert size, as some statistics about the reads.
 * This should allow a certain approximation about the needed hash table size, to avoid (slow and sequential) hash resize function.
 * Todo: Check Endianess of the architecture: reorder bytes if architecture supports big-endianess.
 * @param outDB	is the path to meta-file of the final database. Also the small meta file is written in binary format and includes the final destination of the read-database.
 * @param blocks is the number of blocks the Database is separated in. This means the database is physically one block but indexed in the given number of blocks.
 * @param files contains the meta information about the given set of libraries
 * @return bool Value is 1 if database creation terminated without throwing any exception. 0 Otherwise
 */
int makeDB(char* outDB, int blocks, struct readFiles* files){
	printf("CHECKPOINT: Read raw libraries\n");
	// read fasta, fastq or list of files!!!
	if(files->libNum == 0) return 0;
	int i;
	int readNum=1;
	int readl, readr;
	char* format = (char*)malloc(sizeof(char)*50);

	readsList = (struct reads*)malloc(sizeof(struct reads)*maxreadTotNum);

	for(i=0;i<files->libNum;i++){
		// PE or MP
		if(files[i].rightReads){
			printf("\tRead Paired-End Library (%i)\n",i);
			printf("\tMinInsertSize: %i\n",files[i].minInsert);
			printf("\tMaxInsertSize: %i\n",files[i].maxInsert);
			files[i].avgInsert = (files[i].maxInsert + files[i].minInsert) / 2;
			printf("\tPairedEnd Orientation: %s \n",peOri[files[i].oriPE]);

			files[i].startId = readNum;
			// read left reads
			readl = readNum;
			readr = readNum+1;
			format = strrchr(files[i].leftReads,'.')+1;
			printf("\tLeft Reads: %s\n\tFile format: %s\n",files[i].leftReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readFastA_DB(files[i].leftReads,readl,2);
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readFastQ_DB(files[i].leftReads,readl,2);
			}
			else{
				printf("\tNeither fasta nor fastq. EXIT\n");
				return 0;
			}
			// read right reads
			format = strrchr(files[i].rightReads,'.')+1;
			printf("\tRight Reads: %s\n\tFile format: %s\n",files[i].rightReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA_DB(files[i].rightReads,readr,2) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ_DB(files[i].rightReads,readr,2) + 1;
			}
			else{
				printf("\tNeither fasta nor fastq. EXIT\n");
				return 0;
			}
			files[i].endId = readNum-1;
		}
		// SE
		else{
			printf("\tRead Single-End Library\n");
			format = strrchr(files[i].leftReads,'.')+1;
			printf("\tInFile: %s\n\tFile format: %s\n",files[i].leftReads,format);
			files[i].startId = readNum;
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA_DB(files[i].leftReads,readNum,1) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ_DB(files[i].leftReads,readNum,1) + 1;
			}
			else{
				printf("\tNeither fasta nor fastq. Exit !!!\n");
				return 0;
			}
			files[i].endId = readNum-1;
		}
	}

	printf("\tWrite database to: %s\n",outDB);
	writeDB(outDB,blocks,files);

	return 1;
}

/**
 * writeDB is the writing function for the database creation called directly by makeDB. It writes a meta information file and the real read database.
 * @param outDB	is the name of the meta file of the database.
 * @param blocks is the number of blocks the database in separated in by indexing.
 * @param files	is the internal used representation of the meta information about the libraries under consideration.
 */
void writeDB(char* outDB, int blocks, struct readFiles* files){

	int i = 0;
	int j = 0;
	int temp;
	int len;
	int maxlen = 0;
	int blsize = readTotNum / blocks;
	char* readDB = (char*)malloc(1000);
	strcpy(readDB,outDB);
	strcat(readDB,"_reads.db");

	FILE* db = fopen(readDB,"wb");
	uint64_t wPos=0;
	uint64_t** blocksPos = (uint64_t**)malloc(sizeof(uint64_t*)*2);
	blocksPos[0] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);
	blocksPos[1] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);

	blocksPos[0][0] = wPos;

	printf("CHECKPOINT: Write read DB\n");
	printf("\tNumber of all reads: %i\n",readTotNum);

	for(i=0;i<readTotNum;i++){
		if(i==blsize && j < blocks-1){
			blocksPos[1][j] = wPos;
			j++;
			blocksPos[0][j] = wPos;
			blsize += readTotNum / blocks;
		}
		len = readsList[i].len;
		if(len > maxlen) maxlen = len;
		fwrite(&len,sizeof(int),1,db);
		wPos += sizeof(int);
		if(len){
			if(len>100){
//				printf("Len: %i (i:%i)\n",len,i);
//				printf("ID: %i\n",readsList[i].ID);
			}
			fwrite(&readsList[i].ID,sizeof(int),1,db);
			fwrite(readsList[i].seq,sizeof(char),(len+3)/4,db);
			wPos += (sizeof(int)+((len+3)/4));
		}
	}
	blocksPos[1][j] = wPos;

	fclose(db);

	printf("CHECKPOINT: Write metaData to %s\n",outDB);

	db = fopen(outDB,"wb");

	// Write MetaInfo (Binary)
	maxReadLen = maxlen;
	fwrite(&maxlen,sizeof(int),1,db);
	fwrite(&files->libNum,sizeof(int),1,db);
	fwrite(files,sizeof(struct readFiles),files->libNum,db);
	for(i=0; i<files->libNum;i++){
		temp = strlen(files[i].leftReads);
		fwrite(&temp,sizeof(int),1,db);
		fwrite(files[i].leftReads,sizeof(char),temp,db);
		if(files[i].rightReads){
			temp = strlen(files[i].rightReads);
			fwrite(&temp,sizeof(int),1,db);
			fwrite(files[i].rightReads,sizeof(char),temp,db);
			fwrite(&files[i].minInsert,sizeof(int),1,db);
			fwrite(&files[i].maxInsert,sizeof(int),1,db);
			fwrite(&files[i].avgInsert,sizeof(int),1,db);
			fwrite(&files[i].oriPE,sizeof(int),1,db);
		}
		else{
			temp = 0;
			fwrite(&temp,sizeof(int),1,db);
		}
	}

	// Write Block information
	temp = strlen(readDB);
	fwrite(&temp,sizeof(int),1,db);
	fwrite(readDB,sizeof(char),temp,db);
	fwrite(&readTotNum,sizeof(int),1,db);
	fwrite(&blocks,sizeof(int),1,db);
	fwrite(blocksPos[0],sizeof(uint64_t),blocks,db);
	fwrite(blocksPos[1],sizeof(uint64_t),blocks,db);

	free(blocksPos[0]);
	free(blocksPos[1]);
	free(blocksPos);
	free(readDB);

	fclose(db);
}


void write_filteredFasta(struct readFiles* files, struct reads* reads){
	int i;
	int j;
	int lastD;
	int paired;
	int temp;

	FILE* fastaL;
	FILE* fastaR;
	FILE* fastaU;

	printf("MaxRreadLen: %i\n",maxReadLen);
	char* readseq = (char*)malloc(maxReadLen+1);
	char* fastanameL = (char*)malloc(1000);
	char* fastanameR = (char*)malloc(1000);
	char* fastanameU = (char*)malloc(1000);


	for(i=0; i<files->libNum;i++){
		if(files[i].rightReads) paired = 1;
		else paired = 0;
		temp = strlen(files[i].leftReads);
		strcpy(fastanameL,files[i].leftReads);
		lastD = 0;
		for(j=0;j<temp;j++){
			if(fastanameL[j]=='.') lastD = j;
		}
		if(lastD) fastanameL[lastD] = '\0';
		strcat(fastanameL,"_filter.fasta");
		fastaL = fopen(fastanameL,"w");
		if(paired){
			temp = strlen(files[i].rightReads);
			strcpy(fastanameR,files[i].rightReads);
			lastD = 0;
			for(j=0;j<temp;j++){
				if(fastanameR[j]=='.') lastD = j;
			}
			if(lastD) fastanameR[lastD] = '\0';
			strcat(fastanameR,"_filter.fasta");
			fastaR = fopen(fastanameR,"w");
			temp = strlen(files[i].leftReads);
			strcpy(fastanameU,files[i].leftReads);
			lastD = 0;
			for(j=0;j<temp;j++){
				if(fastanameU[j]=='.') lastD = j;
			}
			if(lastD) fastanameU[lastD] = '\0';
			strcat(fastanameU,"_U_filter.fasta");
			fastaU = fopen(fastanameU,"w");
		}
		for(j=files->startId;j<files->endId;){
			if(paired){
				if(reads[j].len && reads[j+1].len){
					decompressReadSt(reads[j].seq,readseq,reads[j].len);
					fprintf(fastaL,">%i\n",j);
					fprintf(fastaL,"%s\n",readseq);
					decompressReadSt(reads[j+1].seq,readseq,reads[j+1].len);
					fprintf(fastaR,">%i\n",j+1);
					fprintf(fastaR,"%s\n",readseq);
				}
				else{
					if(reads[j].len){
						decompressReadSt(reads[j].seq,readseq,reads[j].len);
						fprintf(fastaU,">%i\n",j);
						fprintf(fastaU,"%s\n",readseq);
					}
					if(reads[j+1].len){
						decompressReadSt(reads[j+1].seq,readseq,reads[j+1].len);
						fprintf(fastaU,">%i\n",j+1);
						fprintf(fastaU,"%s\n",readseq);
					}
				}
				j+=2;
			}
			else{
				if(reads[j].len){
					decompressReadSt(reads[j].seq,readseq,reads[j].len);
					fprintf(fastaL,">%i\n",j);
					fprintf(fastaL,"%s\n",readseq);
				}
				j++;
			}
		}
		fclose(fastaL);
		if(paired){
			fclose(fastaR);
			fclose(fastaU);
		}
	}
	free(readseq);
	free(fastanameL);
	free(fastanameR);
	free(fastanameU);
}

void write_filteredDB(char* outDB, int blocks, struct readFiles* files, struct reads* reads){

	int i = 0;
	int j = 0;
	int temp;
	int len;
	int maxlen = 0;
	int blsize = numreads / blocks;
	char* readDB = (char*)malloc(500);
	strcpy(readDB,outDB);
	strcat(readDB,"_reads.db");

	FILE* db = fopen(readDB,"wb");
	uint64_t wPos=0;
	uint64_t** blocksPos = (uint64_t**)malloc(sizeof(uint64_t*)*2);
	blocksPos[0] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);
	blocksPos[1] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);

	blocksPos[0][0] = wPos;

	printf("CHECKPOINT: Write read DB\n");
//	printf("Number of all reads: %i\n",numreads);

	for(i=0;i<numreads;i++){
		if(i==blsize && j < blocks-1){
			blocksPos[1][j] = wPos;
			j++;
			blocksPos[0][j] = wPos;
			blsize += numreads / blocks;
		}
		len = reads[i].len;
		if(len > maxlen) maxlen = len;
		fwrite(&len,sizeof(int),1,db);
		wPos += sizeof(int);
		if(len){
			fwrite(&reads[i].ID,sizeof(int),1,db);
			fwrite(reads[i].seq,sizeof(char),(len+3)/4,db);
			wPos += (sizeof(int)+((len+3)/4));
		}
	}
	blocksPos[1][j] = wPos;

	fclose(db);

	printf("CHECKPOINT: Write metaData to %s\n",outDB);

	db = fopen(outDB,"wb");

	// Write MetaInfo (Binary)
	maxReadLen = maxlen;
	fwrite(&maxlen,sizeof(int),1,db);
	fwrite(&files->libNum,sizeof(int),1,db);
	fwrite(files,sizeof(struct readFiles),files->libNum,db);
	for(i=0; i<files->libNum;i++){
		temp = strlen(files[i].leftReads);
		fwrite(&temp,sizeof(int),1,db);
		fwrite(files[i].leftReads,sizeof(char),temp,db);
		if(files[i].rightReads){
			temp = strlen(files[i].rightReads);
			fwrite(&temp,sizeof(int),1,db);
			fwrite(files[i].rightReads,sizeof(char),temp,db);
			fwrite(&files[i].minInsert,sizeof(int),1,db);
			fwrite(&files[i].maxInsert,sizeof(int),1,db);
			fwrite(&files[i].avgInsert,sizeof(int),1,db);
			fwrite(&files[i].oriPE,sizeof(int),1,db);
		}
		else{
			temp = 0;
			fwrite(&temp,sizeof(int),1,db);
		}
	}

	// Write Block information
	temp = strlen(readDB);
	fwrite(&temp,sizeof(int),1,db);
	fwrite(readDB,sizeof(char),temp,db);
	fwrite(&numreads,sizeof(int),1,db);
//	printf("readTotNum: %i\n",readTotNum);
//	fwrite(&readTotNum,sizeof(int),1,db);
	fwrite(&blocks,sizeof(int),1,db);
	fwrite(blocksPos[0],sizeof(uint64_t),blocks,db);
	fwrite(blocksPos[1],sizeof(uint64_t),blocks,db);

	free(blocksPos[0]);
	free(blocksPos[1]);
	free(blocksPos);
	free(readDB);

	fclose(db);
}

struct reads* readDB(char* outDB){
	char verbose = 0;
	FILE* metaDB = fopen(outDB,"rb");
	if(verbose) printf("Read Reads from: %s\n",outDB);
	struct reads* reads = NULL;

	int i;
	int temp;
	// Read MetaINFO
	fread(&maxReadLen,sizeof(int),1,metaDB);
	if(verbose) printf("MaxRead: %i\n",maxReadLen);
	fread(&temp,sizeof(int),1,metaDB);
	if(verbose) printf("Number of Libs: %i\n",temp);
	struct readFiles* files = (struct readFiles*)malloc(sizeof(struct readFiles)*temp);
	fread(files,sizeof(struct readFiles),temp,metaDB);
	for(i=0; i<files->libNum;i++){
		if(verbose) printf("EndId: %i\n",files[i].endId);
		numreads = files[i].endId;
		fread(&temp,sizeof(int),1,metaDB);
		files[i].leftReads = (char*)malloc(temp+1);
		fread(files[i].leftReads,sizeof(char),temp,metaDB);
		files[i].leftReads[temp] = '\0';

		fread(&temp,sizeof(int),1,metaDB);
		files[i].rightReads = NULL;
		if(temp){
			files[i].rightReads = (char*)malloc(temp+1);
			fread(files[i].rightReads,sizeof(char),temp,metaDB);
			files[i].rightReads[temp] = '\0';
			fread(&files[i].minInsert,sizeof(int),1,metaDB);
			fread(&files[i].maxInsert,sizeof(int),1,metaDB);
			fread(&files[i].avgInsert,sizeof(int),1,metaDB);
			fread(&files[i].oriPE,sizeof(int),1,metaDB);
			printf("\tMatePair Library -> InsertSize: %i\n",files[i].avgInsert);
			printf("\t\tLeftReads: \t\t%s\n",files[i].leftReads);
			printf("\t\tRightReads: \t\t%s\n",files[i].rightReads);
			printf("\t\tMinInsertSize: \t\t%i\n",files[i].minInsert);
			printf("\t\tMaxInsertSize: \t\t%i\n",files[i].maxInsert);
			printf("\t\tPE-read Orientation: \t%s\n",peOri[files[i].oriPE]);
		}
		else{
			printf("SingleEnd Library\n");
			printf("\tReads:  %s\n",files[i].leftReads);
		}
	}

	for(i=0;i<files->libNum;i++){
		free(files[i].leftReads);
		free(files[i].rightReads);
	}
	free(files);

	fread(&temp,sizeof(int),1,metaDB);
	char* readDBFile = (char*)malloc(temp+1);
	fread(readDBFile,sizeof(char),temp,metaDB);
	readDBFile[temp]='\0';
	if(verbose) printf("Path to readDB: %s\n",readDBFile);
	int readNumber;
	fread(&readNumber,sizeof(int),1,metaDB);
	numreads = readNumber;
	if(verbose) printf("readNumber: %i\n",readNumber);

	fclose(metaDB);

	if(readDBFile) reads = (struct reads*)malloc(sizeof(struct reads)*(readNumber+1));
	FILE* readDB = fopen(readDBFile,"rb");
	if(!readDB){
		printf("Can not open readDB!!!\nAbort\n");
		exit(1);
	}
	free(readDBFile);

	int ID;
	int len;

	for(i=0;i<readNumber;i++){
		reads[i].len = 0;
		reads[i].annotation = NULL;
		reads[i].seq = NULL;
	}

	for(i=0;i<readNumber;i++){
		fread(&len,sizeof(int),1,readDB);
		if(len>=31){
			if(len > maxReadLen) maxReadLen = len;
			fread(&ID,sizeof(int),1,readDB);
			reads[ID].seq = (char*)malloc((len+3)/4);
			fread(reads[ID].seq,sizeof(char),(len+3)/4,readDB);
			reads[ID].ID = ID;
			reads[ID].len = len;
		}
	}

	fclose(readDB);

	return reads;
}

void freeDB(struct reads* reads){
	char verbose = 0;
	if(verbose) printf("Numreads to free: %i\n",numreads);
	for(int i=0;i<=numreads;i++){
		if(reads[i].len){
			free(reads[i].seq);
			free(reads[i].annotation);

		}
	}
	free(reads);
}

char* compressRead(char* read){
	int len = strlen(read);
	if(len < 31) return NULL;
	int i = 0,j = ((len+3)/4)-1;
	int complen = (len+3)/4;
	char* compRead = (char*)malloc(complen);
	char current = 0;

	for(;i<len;i++){
		if(i!=0 && i%4==0){
			compRead[j] = current;
			current = 0;
			j--;
		}
		if(codes[(unsigned char)read[i]]==4){
			if(i>=31){
				read[i] = '\0';
				while(i%4!=0){
					current = current << 2;
					i++;
				}
				compRead[j] = current;
				return compRead;
			}
			else{
//				printf("return NULL, i: %i\n",i);
				free(compRead);
				return NULL;
			}
		}
		current = current << 2;
		current |= codes[(unsigned char)read[i]];
	}

	while(i%4!=0){
		current = current << 2;
		i++;
	}
//	if(complen)
		compRead[j] = current;

	return compRead;
}

char* decompressRead(char* compRead, int len){
	int i = 0 ,j = 0,k = ((len+3)/4)-1;
	char* read = (char*)malloc(len+1);
	if(len==0) return NULL;
	char current = 0;

	for(;i<len;i++){
		if(j%4==0){
			current = compRead[k--];
			j=0;
		}
		read[i] = rev_codes[(current >> (3-j)*2) & 3];
		j++;
	}
	read[len]='\0';
	return read;
}

void decompressReadSt(char* compRead, char* read, int len){
	int i = 0 ,j = 0,k = ((len+3)/4)-1;
	char current = 0;

	for(;i<len;i++){
		if(j%4==0){
			current = compRead[k--];
			j=0;
		}
		read[i] = rev_codes[(current >> (3-j)*2) & 3];
		j++;
	}
	read[len]='\0';
}

int readFastA_DB(char* inFile, int readNum, int jump){
	char verbose = 0;
	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	long filesize;
	long cursize=0;
	int i,n;
	int first = 0;

	FILE *fasta;
	if((fasta = fopen(inFile,"r")) == NULL){
		printf("%s can't be opened\n",inFile);
		exit(EXIT_FAILURE);
	}

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("\tFilesize: %li\n",filesize);

	char* buffer2;
	char* read = (char*)malloc(150000);
	char* comp;
	struct reads* temp;

	while((n = fread(filebuffer,sizeof(char),BUFFER_SIZE,fasta))){
		buffer2 = filebuffer;
		for(i=0;i<n;i++){
			if(filebuffer[i] == '\n'){
				filebuffer[i] = '\0';
				if((*buffer2)=='>'){
					if(first){
						// Save read
						if(verbose) printf("%i: len: %lu\n",readTotNum,strlen(read));
						comp = compressRead(read);
						if(comp){
							readsList[readTotNum].ID = readNum;
							readsList[readTotNum].len = strlen(read);
							readsList[readTotNum++].seq = comp;
						}
						else readsList[readTotNum++].len = 0;
					}
					read[0]='\0';
					if(first) readNum+=jump;
					else first=1;
					if(readTotNum >= maxreadTotNum){
						maxreadTotNum *= 2;
						temp = (struct reads*)realloc(readsList,sizeof(struct reads)*maxreadTotNum);
						if(temp) readsList = temp;
						else {
							printf("Can not reallocate memory\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				else{
					strcat(read,buffer2);
				}
				buffer2 = &filebuffer[i] + 1;
				if(n-i<150000 && n > 150000){
					cursize += i;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}
	// Add last read
	if(first){
		if(verbose) printf("%i: len: %lu\n",readTotNum,strlen(read));
		comp = compressRead(read);
		if(comp){
			readsList[readTotNum].ID = readNum;
			readsList[readTotNum].len = strlen(read);
			readsList[readTotNum++].seq = comp;
		}
		else readsList[readTotNum++].len = 0;
	}

	printf("\n");
	fflush(stdout);
	printf("Processed reads: %i\n",readNum);

	free(read);
	fclose(fasta);
	return readNum;


	return 1;
}

int readFastQ_DB(char* inFile, int readNum, int jump){
	printf("CHECKPOINT: ReadFastQ_DB\n");
	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	long filesize;
	long cursize=0;
	int i,n;
	int first = 0;

	FILE *fasta;
	if((fasta = fopen(inFile,"r")) == NULL){
		printf("Can't open inFile: %s\n",inFile);
		exit(EXIT_FAILURE);
	}

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("\tFilesize: %li\n",filesize);

	char* buffer1, *buffer2;
	char* read = (char*)malloc(150000);
	char* comp;
	struct reads* temp;


	buffer1 = filebuffer;

	first = -1;

	while((n = fread(filebuffer,sizeof(char),BUFFER_SIZE,fasta))){
		buffer2 = filebuffer;
		for(i=0;i<n;i++){
			if(filebuffer[i] == '\n'){
				filebuffer[i] = '\0';
				if(first < 1){
					// Name
					if(!first){
						// Save read (Quality Streams?)
//						printf("Before -> readlen: %i (pos: %p)\n",strlen(read),&read);
						comp = compressRead(read);
//						printf("After -> readlen: %i (pos: %p)\n",strlen(read),&read);
						if(comp){
							readsList[readTotNum].len = strlen(read);
							readsList[readTotNum].ID = readNum;
							readsList[readTotNum].seq = comp;
						}
						else{
							readsList[readTotNum].len = 0;
						}
						readTotNum++;
						readNum+=jump;
					}
					read[0]='\0';
					first=1;
					if(readTotNum >= maxreadTotNum-1){
						maxreadTotNum *= 2;
						temp = (struct reads*)realloc(readsList,sizeof(struct reads)*maxreadTotNum);
						if(temp) readsList = temp;
						else {
							printf("Can not reallocate memory\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				else if(first == 1){
					// Sequence
					strcat(read,buffer2);
					first++;
				}
				else if(first == 2){
					first++;
				}
				else{
					// Quality stream!!!
					first = 0;
				}
				buffer2 = &filebuffer[i] + 1;
				if(first == 0 && n-i < 150000 && n > 150000){
					cursize += i;
					first = 3;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}
	comp = compressRead(read);
	if(comp){
		readsList[readTotNum].len = strlen(read);
		readsList[readTotNum].ID = readNum;
		readsList[readTotNum].seq = comp;
	}
	else{
		readsList[readTotNum].len = 0;
	}
	readTotNum++;
	readNum+=jump;

	free(read);

	return readNum-jump;
}
