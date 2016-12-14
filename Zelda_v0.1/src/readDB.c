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
			printf("Read Paired-End Library (%i)\n",i);
			printf("MinInsertSize: %i\n",files[i].minInsert);
			printf("MaxInsertSize: %i\n",files[i].maxInsert);
			files[i].avgInsert = (files[i].maxInsert + files[i].minInsert) / 2;
			printf("MaxInsertSize: %i\n",files[i].maxInsert);
			printf("PairedEnd Orientation: %s \n",peOri[files[i].oriPE]);

			files[i].startId = readNum;
			// read left reads
			readl = readNum;
			readr = readNum+1;
			format = strrchr(files[i].leftReads,'.')+1;
			printf("Left Reads: %s\nFile format: %s\n",files[i].leftReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readFastA_DB(files[i].leftReads,readl,2);
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readFastQ_DB(files[i].leftReads,readl,2);
			}
			else{
				printf("Neither fasta nor fastq. EXIT\n");
				return 0;
			}
			// read right reads
			format = strrchr(files[i].rightReads,'.')+1;
			printf("Right Reads: %s\nFile format: %s\n",files[i].rightReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA_DB(files[i].rightReads,readr,2) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ_DB(files[i].rightReads,readr,2) + 1;
			}
			else{
				printf("Neither fasta nor fastq. EXIT\n");
				return 0;
			}
			files[i].endId = readNum-1;
		}
		// SE
		else{
			printf("Read Single-End Library\n");
			format = strrchr(files[i].leftReads,'.')+1;
			printf("InFile: %s\nFile format: %s\n",files[i].leftReads,format);
			files[i].startId = readNum;
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA_DB(files[i].leftReads,readNum,1) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ_DB(files[i].leftReads,readNum,1) + 1;
			}
			else{
				printf("Neither fasta nor fastq. Exit !!!\n");
				return 0;
			}
			files[i].endId = readNum-1;
		}
//		files[i].endId = readNum - 1;
	}

//	numreads = readNum - 1;

	printf("Write database to: %s\n",outDB);
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
	int blsize = readTotNum / blocks;
	char* readDB = (char*)malloc(500);
	strcpy(readDB,outDB);
	strcat(readDB,"_reads.db");

	FILE* db = fopen(readDB,"wb");
	uint64_t wPos=0;
	uint64_t** blocksPos = (uint64_t**)malloc(sizeof(uint64_t*)*2);
	blocksPos[0] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);
	blocksPos[1] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);

	blocksPos[0][0] = wPos;

	printf("CHECKPOINT write read DB\n");
	printf("Number of all reads: %i\n",readTotNum);

	for(i=0;i<readTotNum;i++){
		if(i==blsize && j < blocks-1){
			blocksPos[1][j] = wPos;
			j++;
			blocksPos[0][j] = wPos;
			blsize += readTotNum / blocks;
		}
		len = readsList[i].len;
		fwrite(&len,sizeof(int),1,db);
		wPos += sizeof(int);
		if(len){
			fwrite(&readsList[i].ID,sizeof(int),1,db);
			fwrite(readsList[i].seq,sizeof(char),(len+3)/4,db);
			wPos += (sizeof(int)+((len+3)/4));
		}
	}
	blocksPos[1][j] = wPos;

//	for(i = 0; i < blocks; i++){
//		printf("Block %i: From %lu to %lu\n",i,blocksPos[0][i],blocksPos[1][i]);
//	}

	fclose(db);

	printf("CHECKPOINT write metaData to %s\n",outDB);

	db = fopen(outDB,"wb");

	// Write MetaInfo (Do not write Binary)
	fwrite(&files->libNum,sizeof(int),1,db);
	fwrite(files,sizeof(struct readFiles),files->libNum,db);
	for(i=0; i<files->libNum;i++){
		temp = strlen(files[i].leftReads);
		fwrite(&temp,sizeof(int),1,db);
		fwrite(files[i].leftReads,sizeof(char),temp,db);
		if(files[i].rightReads){
//			printf("RightReads (lib %i): %s\n",i,files[i].rightReads);
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

struct reads* readDB(char* outDB){
	FILE* metaDB = fopen(outDB,"rb");
	struct reads* reads = NULL;

	int i;
	int temp;
	// Read MetaINFO
	fread(&temp,sizeof(int),1,metaDB);
	printf("Number of Libs: %i\n",temp);
	struct readFiles* files = (struct readFiles*)malloc(sizeof(struct readFiles)*temp);
	fread(files,sizeof(struct readFiles),temp,metaDB);
	for(i=0; i<files->libNum;i++){
		printf("EndId: %i\n",files[i].endId);
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
			printf("MatePair Library -> Insert: %i\n",files[i].avgInsert);
			printf("\tLeftReads:  %s\n",files[i].leftReads);
			printf("\tRightReads: %s\n",files[i].rightReads);
			printf("\tMinInsertSize: %i\n",files[i].minInsert);
			printf("\tMaxInsertSize: %i\n",files[i].maxInsert);
			printf("\tPE-read Orientation: %s\n",peOri[files[i].oriPE]);
		}
		else{
			printf("SingleEnd Library\n");
			printf("\tLeftReads:  %s\n",files[i].leftReads);
		}
	}

	for(i=0;i<files->libNum;i++){
		free(files[i].leftReads);
		free(files[i].rightReads);
	}
	free(files);

//	int blocks;
	fread(&temp,sizeof(int),1,metaDB);
	char* readDBFile = (char*)malloc(temp+1);
	fread(readDBFile,sizeof(char),temp,metaDB);
	readDBFile[temp]='\0';
	printf("Path to readDB: %s\n",readDBFile);
	int readNumber;
	fread(&readNumber,sizeof(int),1,metaDB);
	printf("numreads: %i\n",numreads);

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
//	char* decomp;

	for(i=0;i<readNumber;i++) reads[i].len = 0;

	for(i=0;i<readNumber;i++){
		fread(&len,sizeof(int),1,readDB);
		if(len){
			if(len > maxReadLen) maxReadLen = len;
			fread(&ID,sizeof(int),1,readDB);
			reads[ID].seq = (char*)malloc((len+3)/4);
			fread(reads[ID].seq,sizeof(char),(len+3)/4,readDB);
//			decomp = decompressRead(seq,len);
//			printf("ID: %i len: %i: Seq: %s\n",ID,len,decomp);
//			free(decomp);
			reads[ID].ID = ID;
			reads[ID].len = len;
		}
	}

	fclose(readDB);

	return reads;
}

void freeDB(struct reads* reads){
	for(int i=1;i<=numreads;i++){
		free(reads[i].seq);
		if(reads[i].annotation) free(reads[i].annotation);
	}
	free(reads);
}

char* compressRead(char* read){
//	printf("Read: %s\n",read);
	int len = strlen(read);
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
		if(codes[(int)read[i]]==-1){
			free(compRead);
			return NULL;
		}
		current = current << 2;
		current |= codes[(int)read[i]];
	}

	while(i%4!=0){
		current = current << 2;
		i++;
	}

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
	printf("Filesize: %li\n",filesize);

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
	printf("CHECKPOINT : ReadFastQ_DB\n");
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
	printf("Filesize: %li\n",filesize);

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
						comp = compressRead(read);
						if(comp){
							readsList[readTotNum].len = strlen(read);
							readsList[readTotNum].ID = readNum;
							readsList[readTotNum].seq = comp;
//							printf("Len: %i\n",readsList[readTotNum].len);
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
//					printf("Load new Block\n");
					cursize += i;
					first = 3;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}
	readsList[readTotNum].len = strlen(read);
	comp = compressRead(read);
	if(comp){
		readsList[readTotNum].ID = readNum;
		readsList[readTotNum].seq = comp;
	}
	readTotNum++;
	readNum+=jump;

	free(read);

	return readNum-jump;
}

void runTest(){
	char* test = (char*)malloc(15);
	strcpy(test,"ACGGTGCACTGGTT");
	char* comp = compressRead(test);
	printf("Test org: %s\n",test);
	char* decomp = decompressRead(comp,14);
	printf("Test com: %s\n",decomp);
}
