/*
 ============================================================================
 Name        : FileReader.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Zelda File- and Parameter-Handling
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "readDB.h"
#include "kmer.h"
#include "FileReader.h"
#include "kmerHash.h"
#include "kmerHash_oa.h"

int readListSize = 100000;
struct hashTable *kmere = NULL;
struct readLink *links = NULL;
int graphSize;
char **readList;
int numreads = 0;



void printUsageAll(){
	printf("Usage:\n");
	printf("\tZelda <Command> -@ <name> [options]\n");
	printf("Commands:\n");
	printf("\tmakedb\t\t Creates Read database in binary format\n");
	printf("\tdb\t\t Make Assembly on the basis of a given database\n");
	printf("\tassem\t\t Creates read database in binary format and makes assembly.\n");
	printf("\n");
}

void printUsageMakeDB(){
	printf("Usage: makedb creates read database in binary format\n");
	printf("\tZelda makedb -@ <name> [options]\n");
	printf("Options:\n");
	printf("\t-@       <char>\t\t\t: Assembly Name, Creates Folder and writes database in the <name>/DB/ subfolder\n");
	printf("\t-bn      <int>\t\t\t: Block number of the read database is split into. Should be > thread number [default: 64]\n");
	printf("\t-t       <int>\t\t\t: Number of threads [default: 1]\n");
	printf("Paired-end reads:\n");
	printf("\t-pe      <char char int int>\t: Paired-end library: [left-reads file, right reads file, minSize, maxSize]\n");
	printf("\t-mp      <char char int int>\t: Mate-pair library:  [left-reads file, right reads file, minSize, maxSize]\n");
	printf("Singe-end reads:\n");
	printf("\t-se      <char>\t\t\t: Single-end reads: [reads-file]\n");
	printf("\n");
}

void printUsageDB(){
	printf("Usage: db make assembly with existing database in existing assembly folder\n");
	printf("\tZelda db -@ <name> [options]\n");
	printf("Options:\n");
	printf("\t-@       <char>\t\t\t: Assembly Name, Creates Folder and writes database in the <name>/DB/ subfolder\n");
	printf("\t-t       <int>\t\t\t: Number of threads\n");
	printf("\t-k       <1..63>\t\t: k-mer size [default 47]\n");
	printf("\t-m       <int>\t\t\t: Min overlap length in DeBruijn graph >=k [default: k-mer size, k]\n");
	printf("\n");
}

void printUsageAssem(){
	printf("Usage: assem creates read database in binary format and makes assembly. Overwrites old database if already existing\n");
	printf("\tZelda assem -@ <name> [options]\n");
	printf("Options:\n");
	printf("\t-@       <char>\t\t\t: Assembly Name, Creates Folder and writes database in the <name>/DB/ subfolder\n");
	printf("\t-bn      <int>\t\t\t: Block number of the read database is split into. Should be > thread number [default: 64]\n");
	printf("\t-t       <int>\t\t\t: Number of threads [default: 1]\n");
	printf("\t-k       <1..63>\t\t: k-mer size [default 47]\n");

	printf("Paired-end reads:\n");
	printf("\t-pe      <char char int int>\t: Paired-end library: [left-reads file, right reads file, minSize, maxSize]\n");
	printf("\t-mp      <char char int int>\t: Mate-pair library:  [left-reads file, right reads file, minSize, maxSize]\n");
	printf("Singe-end reads:\n");
	printf("\t-se      <char>\t\t\t: Single-end reads: [reads-file]\n");
	printf("\n");
}

static inline void errorAbort(struct para* para){
	switch(para->run){
		case 1:  printUsageMakeDB(); break;
		case 2:  printUsageDB(); break;
		case 3:  printUsageAssem(); break;
		default: printUsageAll(); break;
	}
	if(para->assemblyName) free(para->assemblyName);
	if(para->files){
		for(int i=0;i<=para->files->libNum;i++){
			free(para->files[i].leftReads);
			free(para->files[i].rightReads);
		}
		free(para->files);
	}
	if(para->readDB) free(para->readDB);
	if(para->asemblyFolder) free(para->asemblyFolder);
	free(para);
	exit(1);
}

void finished(struct para* para){
	int i;
	if(para->assemblyName) free(para->assemblyName);
	if(para->files){
		for(i=0;i<para->files->libNum;i++){
			free(para->files[i].leftReads);
			free(para->files[i].rightReads);
		}
		free(para->files);
	}
	if(para->readDB) free(para->readDB);
	if(para->asemblyFolder) free(para->asemblyFolder);
	free(para);
	printf("Assembly Successful!!!\n");
	exit(EXIT_SUCCESS);
}

struct para* readCMDline(int argc, char *argv[]){
	struct para* para = (struct para*)malloc(sizeof(struct para));
	// defaults
	para->assemblyName = NULL;
	para->readDB = NULL;
	para->asemblyFolder = NULL;
	para->files = NULL;
	para->blocks = 64;
	para->kSize = 47;
	para->minOvlLen = 0;
	para->threads = 1;
	para->run = 0;
	dtSize = sizeof(KmerBitBuffer)*8;

	int i;
	for(i=0; i<argc; i++){
		if(strcmp(argv[i],"makedb")==0) para->run = 1;
		if(strcmp(argv[i],"db")==0) para->run = 2;
		if(strcmp(argv[i],"assem")==0) para->run = 3;
		if(strcmp(argv[i],"-@")==0 && i+1 < argc && argv[i+1][0] != '-'){
			para->assemblyName = (char*)malloc(strlen(argv[i+1])+1);
			strcpy(para->assemblyName,argv[i+1]);
		}
	}
	if(!para->run || !para->assemblyName){
		errorAbort(para);
	}
	if(para->run == 1 || para->run == 3){
		int libNum = 0;
		for(i=0; i<argc; i++){
			if(strcmp(argv[i],"-se")==0 || strcmp(argv[i],"-pe")==0 || strcmp(argv[i],"-mp")==0) libNum++;
		}
		if(!libNum){
			errorAbort(para);
		}
		para->files = readCMDmakeDB(argc,argv,libNum);
		if(!para->files){
			errorAbort(para);
		}
	}
	if(para->run == 2){
		struct stat st;
		char* temp = (char*)malloc(1000);
		if(stat(para->assemblyName, &st) == -1) {
		    printf("Can not read database. File not found!\n");
		    errorAbort(para);
		}
		else{
			sprintf(temp,"%s/DB",para->assemblyName);
			if(stat(temp,&st) == -1){
			    printf("Can not read database. File not found!\n");
			    errorAbort(para);
			}
			else{
				sprintf(temp,"%s/DB/%s.db",para->assemblyName,para->assemblyName);
				if(stat(temp,&st) == -1){
				    printf("Can not read database. File not found!\n");
				    errorAbort(para);
				}
				else{
					para->readDB = (char*)malloc(strlen(temp)+100);
					strcpy(para->readDB,temp);
				}
			}
		}
		free(temp);
	}
	for(i=0; i<argc; i++){
		if(strcmp(argv[i],"-k")==0 && i+1 < argc && argv[i+1][0] != '-'){
			para->kSize = atoi(argv[i+1]);
			if(para->kSize < 1 || para->kSize+1 > dtSize/2){
				printf("Illegal k-mer Size (0 < k < %i)\n",dtSize/2);
				errorAbort(para);
			}
		}
		if(strcmp(argv[i],"-t")==0 && i+1 < argc && argv[i+1][0] != '-'){
			para->threads = atoi(argv[i+1]);
			if(para->threads < 1){
				printf("Illegal number of Threads\n");
				errorAbort(para);
			}
		}
		if(strcmp(argv[i],"-bn")==0 && i+1 < argc && argv[i+1][0] != '-'){
			para->blocks = atoi(argv[i+1]);
			if(para->blocks < 1){
				printf("Illegal number of DB-Blocks\n");
				errorAbort(para);
			}
		}
		if(strcmp(argv[i],"-m")==0 && i+1 < argc && argv[i+1][0] != '-'){
			para->minOvlLen = atoi(argv[i+1]);
		}
	}
	if(para->minOvlLen < para->kSize) para->minOvlLen = para->kSize;

	if(para->run == 1 || para->run == 3){
		struct stat st;
		char* temp = (char*)malloc(1000);
		if(stat(para->assemblyName, &st) == -1) mkdir(para->assemblyName, 0700);
		sprintf(temp,"%s/DB",para->assemblyName);
		if(stat(temp, &st) == -1) mkdir(temp, 0700);
		sprintf(temp,"%s/DB/%s.db",para->assemblyName,para->assemblyName);
		para->readDB = (char*)malloc(strlen(temp)+100);
		strcpy(para->readDB,temp);
		free(temp);
	}

	if(para->run == 2 || para->run == 3){
		struct stat st;
		char* temp = (char*)malloc(1000);
		sprintf(temp,"%s/Assembly",para->assemblyName);
		if(stat(temp, &st) == -1) mkdir(temp, 0700);
		para->asemblyFolder = (char*)malloc(strlen(temp)+1);
		strcpy(para->asemblyFolder,temp);
		free(temp);
	}

	return para;
}

struct readFiles* readCMDmakeDB(int argc, char *argv[],int libnumTot){
	int i,j;
	int libnum = 0;
	int maxnum = 20;
	struct readFiles* files = (struct readFiles*)malloc(sizeof(struct readFiles) * libnumTot);

	// find all read files and set their parameters
	for(i=0; i<argc;i++){
		// Read Paired-end or mate-pair libs before Single-End // Otherwise its not clear which reads build pairs
		if(strcmp(argv[i],"-pe")==0 || strcmp(argv[i],"-mp")==0){
			for(j=1;j<5;j++){
				if(!argv[i+j] || argv[i+j][0]=='-'){
					printf("Missing PE/MP Parameters\n");
					return NULL;
				}
			}
			printf("Found Paired-End libraries\n");
			files[libnum].leftReads = (char*)malloc(sizeof(char)*100);
			files[libnum].rightReads = (char*)malloc(sizeof(char)*100);
			strcpy(files[libnum].leftReads,argv[i+1]);
			strcpy(files[libnum].rightReads,argv[i+2]);
			files[libnum].minInsert = atoi(argv[i+3]);
			files[libnum].maxInsert = atoi(argv[i+4]);
			i+=4;
			if(argv[i+1] && argv[i+1][0]!='-'){
				if(strcmp(argv[i+1],"fr")==0 || strcmp(argv[i+1],"FR")==0) files[libnum].oriPE = FR;
				else if(strcmp(argv[i+1],"ff")==0 || strcmp(argv[i+1],"FF")==0) files[libnum].oriPE = FF;
				else if(strcmp(argv[i+1],"rf")==0 || strcmp(argv[i+1],"RF")==0) files[libnum].oriPE = RF;
				else{
					printf("Missing PE/MP Parameters\nAbort\n");
					return NULL;
				}
				i++;
			}
			else{
				files[libnum].oriPE = FR;
			}
			libnum++;
		}
		// Read Single-end libs after paired ends --> important
		if(strcmp(argv[i],"-se")==0){
			if(!argv[i+1] || argv[i+1][0]=='-'){
				printf("Missing SE Parameters\n");
				return NULL;
			}
			printf("Found Single-End libraries\n");
			files[libnum].rightReads = NULL;
			files[libnum].leftReads = (char*)malloc(sizeof(char)*100);
			strcpy(files[libnum].leftReads,argv[i+1]);
			libnum++;
			i++;
		}
		if(libnum == maxnum){
			files = (struct readFiles*)realloc(files,sizeof(struct readFiles)*(maxnum*2));
			maxnum*=2;
		}
	}
	files->libNum = libnum;
	return files;
}

#define BUFFER_SIZE (4 * 1024 * 1024)

volatile int readID_global = 1;

void* mt_fileReader(void* block){
	struct hash_block hash_block = *((struct hash_block*)block);
	long start = hash_block.start;
	long end = hash_block.end;
	long len = end - start;
	int pthr_id = hash_block.pthr_id;
	FILE* fasta = fopen(hash_block.fasta,"r");
	fseek(fasta,start,SEEK_CUR);

	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	char* readname = NULL;
	char* readseq;
	size_t n;
	char* buffer2;
	long cursize=0;
	int i;
	int readID;
	int tempnum;

	// Read blocks of 4mb (macro BUFFER_SIZE) from FileStream to the End
	while((n = fread(filebuffer,sizeof(char),BUFFER_SIZE,fasta))){

		buffer2 = filebuffer;
		for(i=0;i<(int)n;i++){
			if(filebuffer[i] == '\n'){
				filebuffer[i] = '\0';
				if((*buffer2)=='>'){
					// Header
					if(!readname) readname = (char*)malloc(1000);
					strcpy(readname,buffer2);
				}
				else{
					// Sequence Data
					if(readname && strlen(buffer2)){
						// Provides consistent assignment of IDs, but assignment is not deterministic (best resolution is a pre-build database)
						tempnum = readID_global;
						do{
							readID = tempnum;
							tempnum = __sync_val_compare_and_swap(&readID_global,readID,readID+1);
						} while(tempnum != readID);
						readseq = buffer2;
						mt_createKmers(readseq,readID);
					}
					if(cursize + i > len){
						printf("Code 0: Thread %i has done a good Job!\n",pthr_id);
						fclose(fasta);
						return NULL;
					}
				}
				buffer2 = &filebuffer[i] + 1;
				if(n-i<5000 && n > 5000){
					cursize += i;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}

	printf("Code 1: Thread %i has done a good Job and reached the end of the file!\n",pthr_id);

	fclose(fasta);
	return NULL;
}

void* mt_fileReaderDB(void* block){
	struct hash_block hash_block = *((struct hash_block*)block);
	uint64_t start = (uint64_t)hash_block.start;
	uint64_t end = (uint64_t)hash_block.end;

	FILE* fasta = fopen(hash_block.fasta,"rb");
	fseek(fasta,start,SEEK_SET);

	char* readsequence = (char*)malloc(150000);
	int readID=0;
	int readLen=0;
	uint64_t wPos = start;

	while(wPos < end){
		fread(&readLen,sizeof(int),1,fasta);
		if(readLen>=31){
			fread(&readID,sizeof(int),1,fasta);
			fread(&readsequence[16],sizeof(char),(readLen+3)/4,fasta);
//			printf("ID: %i, len: %i\n",readID,readLen);
			if(readLen >= nK){
//				printf("ID: %i, len: %i\n",readID,readLen);
				mt_createKmers_DB(&readsequence[16],readLen,readID);
			}
			wPos += (sizeof(int) + ((readLen+3)/4));
		}
		wPos += sizeof(int);
//		printf("wpos: %i, end: %i\n",wPos, end);
	}

	unsigned char old_mutex = fin_mutex;
	unsigned char new_mutex;
	do{
		new_mutex = old_mutex;
		old_mutex = __sync_val_compare_and_swap(&fin_mutex,new_mutex,new_mutex+1);
	} while(old_mutex != new_mutex);

	fclose(fasta);
	free(readsequence);
	return NULL;
}


void fileScheduler(char* inFile, int pthr_num, pthread_t* threads){
	FILE* fasta = fopen(inFile,"r");
	long filesize;
	long start, end;
	long blocksize;
	int i;

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("Filesize: %li\n",filesize);
	// Creates HashTable
	createHashTable_oa();

	struct hash_block* hash_block = (struct hash_block*)malloc(sizeof(struct hash_block)*pthr_num);
	blocksize = filesize/pthr_num;
	start = 0;

	for(i=0;i<pthr_num-1;i++){
		end = (i+1)*blocksize;
		hash_block[i].pthr_id = i;
		hash_block[i].start = start;
		hash_block[i].end = end;
		hash_block[i].fasta = (char*)malloc(strlen(inFile)+1);
		strcpy(hash_block[i].fasta,inFile);
		hash_block[i].pthr_num = pthr_num;
	    pthread_create(&threads[i], NULL, mt_fileReader, (void*)&hash_block[i]);
	    start = end+1;
	}
	end = filesize;
	i = pthr_num-1;
	hash_block[i].pthr_id = i;
	hash_block[i].start = start;
	hash_block[i].end = end;
	hash_block[i].fasta = (char*)malloc(strlen(inFile)+1);
	strcpy(hash_block[i].fasta,inFile);
	hash_block[i].pthr_num = pthr_num;
	pthread_create(&threads[i], NULL, mt_fileReader, (void*)&hash_block[i]);

	fclose(fasta);
	void* status;
    for(i = 0; i < pthr_num; i++){
        pthread_join(threads[i], &status);
    }
    numreads = readID_global-1;
    hashStats_oa();
}
/**
 * Read raw-read libraries in binary database format. Initialize the hash table and start threads to fill
 * the table blockwise in parallel, join threads and return.
 * @param dbFile 	Database file name. Previously build by makeDB. Database contains metainformation about the libs and
 *  a db file containing the row reads
 * @param pthr_num	Number of threads to start in parallel
 * @param threads	Thread Identifier
 * @return			The library meta information (SE-,PE-,MP-lib, insert size, read number, ...)
 */
struct readFiles* fileScheduler_DB(char* dbFile, int pthr_num, pthread_t* threads){
	FILE* metaDB = fopen(dbFile,"rb");

	if(!metaDB){
		printf("Database not Found!\n");
	}

	int i;
	int temp;
	// Read MetaINFO
	fread(&maxReadLen,sizeof(int),1,metaDB);
	fread(&temp,sizeof(int),1,metaDB);
	printf("Number of Libs: %i\n",temp);
	printf("Max Read Len: %i\n",maxReadLen);
	struct readFiles* files = (struct readFiles*)malloc(sizeof(struct readFiles)*temp);
	fread(files,sizeof(struct readFiles),temp,metaDB);
	for(i=0; i<files->libNum;i++){
		printf("EndId: %i\n",files[i].endId);
		numreads = files[i].endId;
		fread(&temp,sizeof(int),1,metaDB);
		files[i].leftReads = (char*)malloc(temp+1);
		fread(files[i].leftReads,sizeof(char),temp,metaDB);
		files[i].leftReads[temp]='\0';
		temp = 0;
		fread(&temp,sizeof(int),1,metaDB);
		files[i].rightReads = NULL;
		if(temp){
			files[i].rightReads = (char*)malloc(temp+1);
			fread(files[i].rightReads,sizeof(char),temp,metaDB);
			files[i].rightReads[temp]='\0';
			fread(&files[i].minInsert,sizeof(int),1,metaDB);
			fread(&files[i].maxInsert,sizeof(int),1,metaDB);
			fread(&files[i].avgInsert,sizeof(int),1,metaDB);
			fread(&files[i].oriPE,sizeof(int),1,metaDB);
			printf("MP/PE Library -> Insert: %i - %i\n",files[i].minInsert,files[i].maxInsert);
			printf("\tLeftReads:  %s\n",files[i].leftReads);
			printf("\tRightReads: %s\n",files[i].rightReads);
		}
		else{
			printf("SingleEnd Library\n");
			printf("\tLeftReads:  %s\n",files[i].leftReads);
		}
	}

	printf("numreads: %i\n",numreads);
	readLenList = (int16_t*)malloc(sizeof(int16_t)*(numreads+1));
	readStartList = (int*)malloc(sizeof(int)*(numreads+1));

	// Read BlockINFO
	// Lock if ReadBlock distribution is consistent with pthr_num
	int blocks;
	fread(&temp,sizeof(int),1,metaDB);
	char* readDBFile = (char*)malloc(temp+1);
	fread(readDBFile,sizeof(char),temp,metaDB);
	readDBFile[temp]='\0';
	printf("Path to readDB: %s\n",readDBFile);
	int readNumber;
	fread(&readNumber,sizeof(int),1,metaDB);

	fread(&blocks,sizeof(int),1,metaDB);
	printf("Number of Block in this DB: %i\n",blocks);
	uint64_t** blocksPos = (uint64_t**)malloc(sizeof(uint64_t*)*2);
	blocksPos[0] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);
	blocksPos[1] = (uint64_t*)malloc(sizeof(uint64_t)*blocks);
	fread(blocksPos[0],sizeof(uint64_t),blocks,metaDB);
	fread(blocksPos[1],sizeof(uint64_t),blocks,metaDB);

	fclose(metaDB);

	createHashTable_oa();
	if(blocks < pthr_num){
		printf("#Threads > #Blocks: Set pthr_num to %i\n",blocks);
		pthr_num = blocks;
	}
	struct hash_block* hash_block = (struct hash_block*)malloc(sizeof(struct hash_block)*pthr_num);

	// Sync #blocks and #threads. Start Threads.
	if(!blocks || !pthr_num){
		if(!blocks) printf("No Blocks Defined\n Abort!!!\n");
		if(!pthr_num) printf("Number of threads set to 0\n Abort!!!\n");
		exit(EXIT_FAILURE);
	}
	int blocks_per_thread;
	int rest_blocks;
	int blockst = 0;
	int blockend;
	rest_blocks = blocks%pthr_num;
	blocks_per_thread = blocks/pthr_num;
	pthr_runN = pthr_num;

	printf("Init Number of Threads: %i in %i blocks\n",pthr_num,blocks);
	printf("Number of blocks per Thread: %i (rest: %i)\n",blocks_per_thread,rest_blocks);
	for(i=0;i<pthr_num;i++){
		hash_block[i].pthr_id = i;
		hash_block[i].pthr_num = pthr_num;
		hash_block[i].fasta = (char*)malloc(strlen(readDBFile)+1);
		strcpy(hash_block[i].fasta,readDBFile);
		hash_block[i].start = blocksPos[0][blockst];
		blockend = (blockst + blocks_per_thread) - 1;
		if(rest_blocks){
			blockend++;
			rest_blocks--;
		}
		hash_block[i].end = blocksPos[1][blockend];
	    pthread_create(&threads[i], NULL, mt_fileReaderDB, (void*)&hash_block[i]);
		blockst = blockend+1;
	}

	void* status;
    for(i = 0; i < pthr_num; i++){
        pthread_join(threads[i], &status);
    }

    hashStats_oa();

    for(i = 0; i < pthr_num; i++){
        free(hash_block[i].fasta);
    }
    free(hash_block);
    free(blocksPos[0]);
    free(blocksPos[1]);
    free(blocksPos);
    free(readDBFile);

	return files;
}

void freeFiles(struct para* para){
	if(!para->files) return;
	for(int i=0; i<para->files->libNum; i++){
		free(para->files[i].leftReads);
		free(para->files[i].rightReads);
	}
	free(para->files);
	para->files = NULL;
}
