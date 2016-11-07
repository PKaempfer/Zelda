/*
 * FileReader.c
 *
 *  Created on: Oct 13, 2014
 *      Author: kaempfpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "kmer.h"
#include "FileReader.h"
//#include "uthash.h"
#include "kmerHash.h"
//#include "jellyhash.h"
#include "kmerHash_oa.h"

// Use Jellyfish master (C-Implementation of concurrent k-mer hashing)
//#define JELLY 1

#ifdef JELLY
	JellyHash jhash;
	HKey hpos;
#endif

int readListSize = 100000;
struct hashTable *kmere = NULL;
struct readLink *links = NULL;
int graphSize;
char **readList;
int numreads = 0;

void hash_chain_len_histogram(UT_hash_table *tbl){
	printf("Histogram\n");
	unsigned i, bkt_hist[CHAIN_MAX+1];
	double pct = 100.0/tbl->num_buckets;
	memset(bkt_hist,0,sizeof(bkt_hist));
	for(i=0; i < tbl->num_buckets; i++) {
		unsigned count = tbl->buckets[i].count;
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
}



void setLink(char* kmer, int readNum, int end){
	struct readLink *s;
	readID *temps;
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;

	temp = toBuffer(kmer);
	revtemp = revKmer(temp);
	if(revtemp < temp) temp = revtemp;
	HASH_FIND(hhb, links, &temp, sizeof(KmerBitBuffer),s);
	if(s==NULL){
		s = (struct readLink*)malloc(sizeof(struct readLink));
		s->kmer = temp;
		s->read = (readID*)malloc(sizeof(readID)*17);
		s->read[0] = setID(readNum, end);
		s->read[1] = 0;
		HASH_ADD(hhb, links, kmer, sizeof(KmerBitBuffer), s);
	}
	else{
		int i=1;
		while(s->read[i]!=0){
			i++;
		}
		s->read[i] = setID(readNum,end);
		if(i%16==0){
//			printf("Realloc readLink: %i\n",i);
			temps = (readID*)realloc(s->read,sizeof(readID)*(i+17));
			s->read = temps;
//			printf("Realloc SetRead\n");
//			s->read[i+1] = 0;
//			printf("Realloc readLinkAFTER\n");
		}
//		else{
//			if(i>=16) printf("i: %i\n",i);
		s->read[i+1] = 0;
//		}
	}
}

void iterLink(){
    struct readLink *s;
    s = (struct readLink*)malloc(sizeof(struct readLink));
    int i;
    int j=0;
    int begin=0;
    int end = 0;
    for(s=links; s != NULL; s=(struct readLink*)(s->hhb.next)) {
    	i=0;
    	while(s->read[i]!=0){
    		if(__builtin_clz(s->read[i])){
    			begin++;
    		}
    		else{
    			end ++;
    		}
    		i++;
    		j++;
    	}
    }
    printf("num in LinkList: %i\n",j);
    printf("Beg in LinkList: %i\n",begin);
    printf("End in LinkList: %i\n",end);
}

void createKmers(char* read, int readNum){
//	static clock_t times,timee=0,timeh=0;
#ifdef UT_HASH
	struct hashTable *s;
#elif JELLY
//	printf("Use JellyFisch\n");
	char kmer[nK];
	char out[nK];
	KmerBitBuffer out_bit;
	size_t key;
	int inserted = 10;
#else

#endif

	static int len;
	static KmerBitBuffer temp;
	static KmerBitBuffer revtemp;
	static char* buffer;
	if(strchr(read,'N')){
//		printf("Read contains Gaps\n");
		return;
	}
	len = strlen(read);
	if(len > maxRlen) maxRlen = len;
	buffer = read+(len-nK);
	len-=nK;
	if(len<0) return;
//	printf("setLink (end): %s\n",buffer);
	setLink(buffer,readNum,1);
	while(len!=-1){
		*(buffer+nK)='\0';
//		times = clock();
		temp = toBuffer(buffer);
		revtemp = revKmer(temp);
		if(revtemp < temp) temp = revtemp;
#ifdef UT_HASH
		HASH_FIND( hhb ,kmere, &temp, sizeof(KmerBitBuffer), s);
		if(s==NULL){
			s = (struct hashTable*)malloc(sizeof(struct hashTable));
			s->kmer = temp;
			s->trans = 0;
			s->index = 0;
			HASH_ADD(hhb, kmere, kmer, sizeof(KmerBitBuffer), s);
		}
#elif JELLY
//		strcpy(kmer,toSeq(temp));
//		printf("Insert kmer: %s\n",kmer);
		hpos = jelly_hash_find(&jhash, &temp, 2, &inserted);
//		hpos = jelly_hash_find(&jhash, &temp, 0, &inserted);
		if(hpos == JHASH_NULL){
			printf("Resize Hash Table\n");
			printf("Entry not found!\n");
			bitnum++;
			jelly_hash_alloc(&jhash, bitnum, 20, 2*nK);
			hpos = jelly_hash_find(&jhash, &temp, 2, &inserted);
			if(hpos == JHASH_NULL){
				printf("Entry not found!\n");
			}
		}
		else{
//			jelly_hash_get_key(&jhash, hpos, &out_bit); // get value for `key`
//			printf("Found value at position x: %s\n", toSeq(out_bit));
		}
#else
		insertKmer(&temp);
#endif

		len--;
		buffer--;
	}
//	printf("setLink (start): %s\n",read);
	setLink(read,readNum,0);
//	if(readNum%100000==0) printf("Time for Kmer processing: %.3f sec -> HASH: %.3f sec\n",(float)timee/CLOCKS_PER_SEC,((float)timeh/CLOCKS_PER_SEC)-(float)timee/CLOCKS_PER_SEC);
}

void cleanGraph(){
    struct hashTable *s, *t;
    s = (struct hashTable*)malloc(sizeof(struct hashTable));
    HASH_ITER(hhb, kmere, s, t){
    	HASH_DELETE(hhb, kmere, s);
    	free(s);
    }
    free(t);
    free(kmere);

//    struct readLink *u,*v;
//    u = (struct readLink*)malloc(sizeof(struct readLink));
//    HASH_ITER(hhb,links,u,v){
//    	free(u->read);
//    	free(u);
//    }
//    free(v);
//    free(links);
}

void printUsage(){
	puts("Usage:");
	puts("./DBSGA [options] -o <output_dir>");
	puts("Options:");
	puts("Create Database:");
	puts("\t-makedb  <char>\t\t\t: Creates database. [Path to MetaDB]");
	puts("\t-bn      <int>\t\t\t: Block number of the read database is split into [default: 64]");
	puts("Make assembly:");
	puts("\t-db      <char>\t\t\t: Reads a compressed database with all given libraries");
	puts("\t-t       <int>\t\t\t: Number of threads for parallelized hash table construction");
	puts("\t-k       <1..63>\t\t: k-mer size [default 47]");
	puts("\t-m       <int>\t\t\t: Min overlap length in DeBruijn graph [k-mer size]");
	puts("Paired-end reads:");
	puts("\t-pe      <char char int int>\t: Paired-end library: [left-reads file, right reads file, minSize, maxSize]");
	puts("\t-mp      <char char int int>\t: Mate-pair library: [left-reads file, right reads file, minSize, maxSize]");
	puts("Singe-end reads:");
	puts("\t-se      <char>\t\t\t: Single-end reads: [reads-file]");
	puts("");

}

struct readFiles* readCMDline(int argc, char *argv[], struct readFiles *files){
	int i,j;
	int libnum = 0;
	int maxnum = 20;
	files = (struct readFiles*)malloc(sizeof(struct readFiles) * maxnum);

	// find all read files and set their parameters
	for(i=0; i<argc;i++){
		// Read Paired-end or mate-pair libs before Single-End // Otherwise its not clear which reads build pairs
		if(strcmp(argv[i],"-pe")==0 || strcmp(argv[i],"-mp")==0){
			for(j=1;j<5;j++){
				if(!argv[i+j] || argv[i+j][0]=='-'){
					printf("Missing PE/MP Parameters\n");
					printUsage();
					exit(1);
				}
			}
			printf("Found Paired-End libraries\n");
			files[libnum].leftReads = (char*)malloc(sizeof(char)*100);
			files[libnum].rightReads = (char*)malloc(sizeof(char)*100);
			strcpy(files[libnum].leftReads,argv[i+1]);
			strcpy(files[libnum].rightReads,argv[i+2]);
//			files[libnum].insertSize = atoi(argv[i+3]);
			files[libnum].minInsert = atoi(argv[i+3]);
			files[libnum].maxInsert = atoi(argv[i+4]);
			i+=4;
			if(argv[i+1] && argv[i+1][0]!='-'){
				if(strcmp(argv[i+1],"fr")==0 || strcmp(argv[i+1],"FR")==0) files[libnum].oriPE = FR;
				else if(strcmp(argv[i+1],"ff")==0 || strcmp(argv[i+1],"FF")==0) files[libnum].oriPE = FF;
				else if(strcmp(argv[i+1],"rf")==0 || strcmp(argv[i+1],"RF")==0) files[libnum].oriPE = RF;
				else{
					printf("Missing PE/MP Parameters\nAbort\n");
					printUsage();
					exit(1);
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
				printUsage();
				exit(1);
			}
			printf("Found Single-End libraries\n");
			files[libnum].rightReads = NULL;
			files[libnum].leftReads = (char*)malloc(sizeof(char)*100);
			strcpy(files[libnum].leftReads,argv[i+1]);
			files[libnum].insertSize = 0;
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

int readFile(struct readFiles* files){
	// read fasta, fastq or list of files!!!
	int i;
	int readNum=1;
	int readl, readr;
	char* format = (char*)malloc(sizeof(char)*50);

#ifdef JELLY
	jelly_hash_alloc(&jhash, bitnum, 20, 2*nK);
#else
	createHashTable();
#endif


	if(files->libNum == 0) return 0;

	for(i=0;i<files->libNum;i++){
		// PE or MP
		if(files[i].rightReads){
			printf("Read Paired-End Library\n");
			printf("InsertSize: %i\n",files->insertSize);
			files[i].startId = readNum;
			// read left reads
			readl = readNum;
			readr = readNum+1;
			format = strrchr(files[i].leftReads,'.')+1;
			printf("InFile: %s\nFile format: %s\n",files[i].leftReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readFastA(files[i].leftReads,readl,2);
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readFastQ(files[i].leftReads,readl,2);
			}
			else{
				printf("Neither fasta nor fastq. EXIT\n");
				return 0;
			}
			// read right reads
			format = strrchr(files[i].rightReads,'.')+1;
			printf("InFile: %s\nFile format: %s\n",files[i].rightReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA(files[i].rightReads,readr,2) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ(files[i].rightReads,readr,2) + 1;
			}
			else{
				printf("Neither fasta nor fastq. EXIT\n");
				return 0;
			}
		}
		// SE
		else{
			format = strrchr(files[i].leftReads,'.')+1;
			printf("InFile: %s\nFile format: %s\n",files[i].leftReads,format);
			if(strcmp(format,"fasta") == 0 || strcmp(format,"fa") == 0){
				readNum = readFastA(files[i].leftReads,readNum,1) + 1;
			}
			else if(strcmp(format,"fastq") == 0 || strcmp(format,"fq") == 0){
				readNum = readFastQ(files[i].leftReads,readNum,1) + 1;
			}
			else{
				printf("Neither fasta nor fastq. Exit !!!\n");
				return 0;
			}
		}
		files[i].endId = readNum - 1;
	}

	numreads = readNum - 1;

#ifdef UT_HASH
	unsigned key_count = HASH_CNT(hhb,kmere);
	fprintf(stderr,"Number of unique k-mers: %i from %u reads\n", key_count, readNum);
	graphSize = key_count+1;
	hash_chain_len_histogram(kmere->hhb.tbl);
#elif JELLY
	FILE* jellystat = fopen("jelly.stat","w");
	jelly_hash_print_stats(&jhash,jellystat);
	fclose(jellystat);
#else
	printf("Myhash:\n");
	hashStats();
	graphSize = itemNum+1;
	printf("GraphSize:%i\n",graphSize);
#endif
	return 1;
}

#define BUFFER_SIZE (4 * 1024 * 1024)

int readFastA(char* inFile, int readNum, int jump){
//	printf("hashTableSize: %i (hashHandle: %i)\n",(int)sizeof(struct hashTable),sizeof(UT_hash_handle));
	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	int start;
	long filesize;
	long cursize=0;
	int i,j,n;
	int first = 0;

	FILE *fasta;
	if((fasta = fopen(inFile,"r")) == NULL){
		printf("%s can't be opend\n",inFile);
		exit(EXIT_FAILURE);
	}

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("Filesize: %li\n",filesize);

	char *buffer1, *buffer2;
	char *read = (char*)malloc(100000);
	if (readLenList == NULL) readLenList = (int*)malloc(sizeof(int)*readListSize);
	if (readStartList == NULL) readStartList = (int*)malloc(sizeof(int)*readListSize);
	//	readList = (char**)malloc(sizeof(char*)*readListSize);

//	char **readListTemp;
	int* readLenListTemp;
	buffer1 = filebuffer;

	while((n = fread(filebuffer,sizeof(char),BUFFER_SIZE,fasta))){
		buffer2 = filebuffer;
		for(i=0;i<n;i++){
			if(filebuffer[i] == '\n'){
				filebuffer[i] = '\0';
				if((*buffer2)=='>'){
					if(first){
						if(readNum >= readListSize-10){
							readLenListTemp = (int*)realloc(readLenList,sizeof(int)*(readListSize*2));
							readLenList = readLenListTemp;
							readLenListTemp = (int*)realloc(readStartList,sizeof(int)*(readListSize*2));
							readListSize*=2;
							readStartList = readLenListTemp;
						}
						readLenList[readNum] = strlen(read)-nK;
						readStartList[readNum] = start;
						if(readLenList[readNum] < 0) continue;

						if(!((readNum/jump)%10000)){
							printf("Processed reads: %i\r",readNum/jump);
							fflush(stdout);
						}
						createKmers(read, readNum);
					}
					for(j=1;j<1000;j++){
						if (buffer2[j]=='(') break;
					}
					sscanf(&buffer2[j],"(Strand %*c Offset %i--%*i)",&start);
					read[0]='\0';
					if(first) readNum+=jump;
					else first=1;
				}
				else{
					strcat(read,buffer2);
				}
				buffer2 = &filebuffer[i] + 1;
				if(n-i<50000 && n > 50000){
					cursize += i;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}

	readLenList[readNum] = strlen(read)-nK;
	readStartList[readNum] = start;
	createKmers(read, readNum);
	printf("\n");
	fflush(stdout);
	printf("Processed reads: %i\n",readNum);

	free(read);
	fclose(fasta);
	return readNum;
}


int readFastQ(char* inFile, int readNum, int jump){
	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	int start;
	long filesize;
	long cursize=0;
	int i,j,n;
	int first = 0;

	FILE *fasta;
	if((fasta = fopen(inFile,"r")) == NULL){
		printf("%s can't be opend\n",inFile);
		exit(EXIT_FAILURE);
	}

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("Filesize: %li\n",filesize);

	char *buffer1, *buffer2;
	char *read = (char*)malloc(70000);
	if (readLenList == NULL) readLenList = (int*)malloc(sizeof(int)*readListSize);
	if (readStartList == NULL) readStartList = (int*)malloc(sizeof(int)*readListSize);
	//	readList = (char**)malloc(sizeof(char*)*readListSize);

//	char **readListTemp;
	int* readLenListTemp;
	buffer1 = filebuffer;


	while((n = fread(filebuffer,sizeof(char),BUFFER_SIZE,fasta))){
		buffer2 = filebuffer;
		for(i=0;i<n;i++){
			if(filebuffer[i] == '\n'){
				filebuffer[i] = '\0';
				if((*buffer2)=='@'){
					if(first){
						if(readNum >= readListSize-10){
							readLenListTemp = (int*)realloc(readLenList,sizeof(int)*(readListSize*2));
							readLenList = readLenListTemp;
							readLenListTemp = (int*)realloc(readStartList,sizeof(int)*(readListSize*2));
							readListSize*=2;
							readStartList = readLenListTemp;
						}
						if(!((readNum/jump)%10000)){
							printf("Processed reads: %i\r",readNum/jump);
							fflush(stdout);
						}
//						printf("Read: %s\n",read);
						readLenList[readNum] = strlen(read)-nK;
						readStartList[readNum] = start;
						first=1;
						if(readLenList[readNum] > 0) createKmers(read, readNum);
					}
					for(j=1;j<1000;j++){
						if (buffer2[j]=='(') break;
					}
					sscanf(&buffer2[j],"(Strand %*c Offset %i--%*i)",&start);
					read[0]='\0';
					if(first) readNum+=jump;
					else first=1;
				}
				else if((*buffer2)=='+'){
//					printf("buffer2: %s\n",buffer2);
					first = 2;
				}
				else if(first == 1){
					strcat(read,buffer2);
//					printf("CatRead: %s\n",read);
				}
				buffer2 = &filebuffer[i] + 1;
				if(n-i<1000 && n > 1000){
					cursize += i;
					fseek(fasta, -(n-i), SEEK_CUR);
					break;
				}
			}
		}
	}
	return readNum-jump;
}

volatile int readID_global = 1;

void* mt_fileReader(void* block){
	struct hash_block hash_block = *((struct hash_block*)block);
	long start = hash_block.start;
	long end = hash_block.end;
	long len = end - start;
	int pthr_id = hash_block.pthr_id;
//	int pthr_num =  hash_block.pthr_num;
	FILE* fasta = fopen(hash_block.fasta,"r");
	fseek(fasta,start,SEEK_CUR);
//	printf("Thread %i starts reading pos: %ld %ld\n",pthr_id,start, end);

	char filebuffer[BUFFER_SIZE]; // 4 MB buffer
	char* readname = NULL;
	char* readseq;
	size_t n;
	char* buffer2;
	long cursize=0;
	int i;
	int readID;
//	int bits = bitOffset(pthr_num);
	int tempnum;
//	int newnum;

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
					// Sequence Dataas
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
//	long len = end - start;
	int pthr_id = hash_block.pthr_id;
//	int pthr_num =  hash_block.pthr_num;

	FILE* fasta = fopen(hash_block.fasta,"rb");
	fseek(fasta,start,SEEK_SET);

	char* readsequence = (char*)malloc(150000);
//	char* decomp;
//	size_t n;
//	int i;
	int readID=0;
	int readLen=0;
	uint64_t wPos = start;

	printf("Thread: %i opens DB: %s (start: %ld / End: %ld)\n",pthr_id,hash_block.fasta,start, end);

	while(wPos < end){
		fread(&readLen,sizeof(int),1,fasta);
		if(readLen){
			fread(&readID,sizeof(int),1,fasta);
			fread(&readsequence[16],sizeof(char),(readLen+3)/4,fasta);
			mt_createKmers_DB(&readsequence[16],readLen,readID);
			wPos += (sizeof(int) + ((readLen+3)/4));
		}
		wPos += sizeof(int);
	}

	unsigned char old_mutex = fin_mutex;
	unsigned char new_mutex;
	do{
		new_mutex = old_mutex;
		old_mutex = __sync_val_compare_and_swap(&fin_mutex,new_mutex,new_mutex+1);
	} while(old_mutex != new_mutex);
	if(fin_mutex>=0 )printf("Thread finished job and reports Finished-Status\n");

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
//	int maxlayer = 20;

	fseek (fasta , 0 , SEEK_END);
	filesize = ftell(fasta);
	rewind(fasta);
	printf("Filesize: %li\n",filesize);
	// Creates HashTable
	// Max number of hash table extensions
	createHashTable_oa();
//	pthread_mutex_init(&lock,NULL);

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
//    pthread_mutex_destroy(&lock);
    // Merge hash Tables

    numreads = readID_global-1;
//    exit(1);
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
	fread(&temp,sizeof(int),1,metaDB);
	printf("Number of Libs: %i\n",temp);
	struct readFiles* files = (struct readFiles*)malloc(sizeof(struct readFiles)*temp);
	fread(files,sizeof(struct readFiles),temp,metaDB);
	for(i=0; i<files->libNum;i++){
		printf("EndId: %i\n",files[i].endId);
		numreads = files[i].endId;
		fread(&temp,sizeof(int),1,metaDB);
//		printf("Left-read malloc size: %i\n",temp);
		files[i].leftReads = (char*)malloc(temp+1);
		fread(files[i].leftReads,sizeof(char),temp,metaDB);
		files[i].leftReads[temp]='\0';
		temp = 0;
		fread(&temp,sizeof(int),1,metaDB);
		files[i].rightReads = NULL;
		if(temp){
//			printf("right-read malloc size: %i\n",temp);
			files[i].rightReads = (char*)malloc(temp+1);
			fread(files[i].rightReads,sizeof(char),temp,metaDB);
			files[i].rightReads[temp]='\0';
			fread(&files[i].minInsert,sizeof(int),1,metaDB);
			fread(&files[i].maxInsert,sizeof(int),1,metaDB);
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
//	fin_mutex = -1;
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
		printf("Thread %i reads blocks %i -> %i\n",i,blockst,blockend);
	    pthread_create(&threads[i], NULL, mt_fileReaderDB, (void*)&hash_block[i]);
		blockst = blockend+1;
	}

//    do{
//    	sleep(5);
//    	printf("Waited another 5 seconds -> That means a new STATS (YEAHH)\n");
//    	if(resize_mutex){
//    		while(resize_mutex);
//    	}
//    	hashStats_oa();
//    } while(fin_mutex != pthr_num);

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
