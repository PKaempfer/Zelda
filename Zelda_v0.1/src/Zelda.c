/*
 ============================================================================
 Name        : Zelda.c
 Author      : Kämpfer, Philipp
 Version     :
 Copyright   : GPLv3 (general public license)
 Description : The Zelda Genome Assembler in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include "kmer.h"
#include "kmerHash.h"
#include "kmerHash_oa.h"
#include "FileReader.h"
#include "DBGraph.h"
#include "DBGraph_oa.h"
#include "DBGraph_stringer.h"
#include "DBGraph_reduced.h"
#include "DBGraph_error.h"
#include "DBGraph_scaffold.h"
#include "readDB.h"
#include "ConsensusCaller.h"

int main(int argc, char* argv[]) {
	char findotdump = 1;
	char dotdump = 0;
	char sleeptime = 0;
	int num_threads = 1;
	time_t start,stop;

	if(argc<3){
		printUsage();
		return EXIT_SUCCESS;
	}

	int i;
	char* fastaFile = NULL;
	char* outDB = NULL;
	char* inDB = NULL;
	char* assemName = NULL;
	char* pathDB = (char*)malloc(1000);
	char* pathAssembly = (char*)malloc(1000);
	char* tempPath = (char*)malloc(1000);
	int blocks = 64;
	dtSize = sizeof(KmerBitBuffer)*8;
	minovlLen=0;

	for(i=0; i<argc ; i++){
		if(strcmp(argv[i],"-@")==0){
			assemName = (char*)malloc(1000);
			strcpy(assemName,argv[i+1]);
			struct stat st = {0};
			if (stat(assemName, &st) == -1) {
			    mkdir(assemName, 0700);
			}
			sprintf(pathDB,"%s/DB",assemName);
			sprintf(pathAssembly,"%s/assembly",assemName);
		}
		if(strcmp(argv[i],"-k")==0){
			nK = atoi(argv[i+1]);
			if(nK+1 > dtSize/2){
				printf("k can not be bigger than %i, set k to %i\n",(dtSize/2)-1,(dtSize/2)-1);
				nK = (dtSize/2) -1;

			}
			else printf("k = %i\n",nK);
		}
		if(strcmp(argv[i],"-in")==0){
			fastaFile=(char*)malloc(500);
			strcpy(fastaFile,argv[i+1]);
			printf("Input File: %s\n",fastaFile);
		}
		if(strcmp(argv[i],"-m")==0){
			minovlLen = atoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-t")==0){
			num_threads = atoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-bn")==0){
			blocks = atoi(argv[i+1]);
		}
		if(strcmp(argv[i],"-makedb")==0){
			outDB = (char*)malloc(1000);
			strcpy(outDB,argv[i+1]);
			if(assemName){
				struct stat st = {0};
				if (stat(pathDB, &st) == -1) {
				    mkdir(pathDB, 0700);
				}
			}
			else{
				printUsage();
				exit(1);
			}
		}
		if(strcmp(argv[i],"-db")==0){
			inDB = (char*)malloc(1000);
			strcpy(inDB,argv[i+1]);
			if(assemName){
				struct stat st = {0};
				if (stat(pathAssembly, &st) == -1) {
				    mkdir(pathAssembly, 0700);
				}
			}
			else{
				printUsage();
				exit(1);
			}
		}
	}

	const int NUM_THREADS = num_threads;
    pthread_t threads[NUM_THREADS];
    struct readFiles *files = NULL;

	if(inDB){
		if(blocks < 1){
			printf("Number of blocks have to be > 0\n");
			exit(1);
		}
		printf("CHECKPOINT: read Fasta File -> create Hash Table\n");
		time(&start);
		files = fileScheduler_DB(inDB,NUM_THREADS,threads);
		printf("CHECKPOINT: create Hash Table END\n");
		time(&stop);
		printf("DBG Time: %0.2f\n",difftime (stop,start));
		sleep(sleeptime);
	}

	if(!files){
		files = readCMDline(argc,argv,files);
	}


	// Generate database
	if(outDB){
		if(blocks < 1){
			printf("Number of blocks have to be > 0\n");
			exit(1);
		}
		makeDB(outDB, blocks, files);
		free(outDB);
		exit(1);
	}

	if(fastaFile){
		printf("CHECKPOINT: read Fasta File -> create Hash Table\n");
		time(&start);
		fileScheduler(fastaFile,NUM_THREADS,threads);
		printf("CHECKPOINT: create Hash Table END\n");
		time(&stop);
		printf("Time: %0.2f\n",difftime (stop,start));
//		exit(1);
	}


	printf("CHECKPOINT: create Graph\n");
	time(&stop);
	createGraph(graphSize);
	time(&start);
	printf("Time: %0.2f\n",difftime (start,stop));

	printf("CHECKPOINT: Fill Adjacency List\n");
	hashToTabDFS_oa();
	time(&stop);
	printf("Time: %0.2f\n",difftime (stop,start));
	printf("Wait after Hashing\n");
	sleep(sleeptime);
	printf("Continue\n");

	printf("HashTable Destructor\n");
	if(inDB){
		freeHashTable_oa();
	}
	printf("Wait after Freeing HashTbale\n");
	sleep(sleeptime);
	printf("Continue\n");

//	printf("CHECKPOINT: Print Graph\n");
//	printGraph();

	printf("CHECKPOINT: Write Dot-File\n");
//	writeDot(NULL);

//	printf("CHECKPOINT: Clean Hash Table\n");
//	cleanGraph();

	printf("CHECKPOINT: Error Correction Table\n");
	time(&start);
	perfectErrorCorrection();// errorCorrection();
//	writeDot("./output/test_aftererror.dot");
	time(&stop);
	printf("Error Correction Time: %0.2f\n",difftime (stop,start));
	printf("Wait after Error correction\n");
	sleep(sleeptime);
	printf("continue\n");

	printf("CHECKPOINT: Prepare Reduction\n");
	time(&start);
	prepareGraph();
	printf("CHECKPOINT: Reduce Graph (light)\n");
//	printGraph();
	travToRedOVL_v2();
	freeGraph();
	labelGraph();
	if(dotdump) printRedDot("output/redGraph.dot");
//	printRedGraph();
//	exit(1);
//	printRedGraph();
	reduceRedGraph();
//	printRedDot("output/RedredGraph.dot");
//	printRedGraph();
	verticalReduction();
//	printRedGraph();
//	printRedDot("output/redRedredGraph.dot");
	do{
		printf("Horizontal reduction\n");
		reduceRedGraph();
		printf("Vertical reduction\n");
	}while(verticalReduction());
	printf("Write dot-File: Reduced_DBG.dot");
	if(findotdump){
		sprintf(tempPath,"%s/Reduced_DBG.dot",pathAssembly);
		printRedDot(tempPath);
	}
	countRemainingNodes();
//	printRedGraph();
	printf("CHECKPOINT: Reduce Graph (strong)\n");
	if(dotdump) printRedGraphToFile("./output/redGraphBefore.list");
	do{
		printf("Graph reduction\n");
		reduceRedGraph_strong();
	}while(verticalReduction());
	if(findotdump){
		sprintf(tempPath,"%s/Reduced_DBG_strong.dot",pathAssembly);
		printRedDot(tempPath);
	}
	time(&stop);
	printf("Reducer Time: %0.2f\n",difftime (stop,start));
	printf("Wait after Reduction\n");
	sleep(sleeptime);
	printf("continue\n");
	countRemainingNodes();
	redGraphConnector();
	if(dotdump){
		printRedGraphToFile("./output/redGraphAfter.list");
	}

//	testFunct();

//	printRedGraph();
//	printRedGraphToFile("./output/redGraph.list");
//	exit(1);
	printf("CHECKPOINT: Stringer:\n");
	time(&start);
	printf("Numread: %i\n",numreads);
	struct myovlList* G = initOVLgraph(numreads);
	struct string_graph* S = initStringGraph(G);
//	exit(1);
	time(&stop);
	printf("Stringer: %0.2f\n",difftime (stop,start));
	printf("Wait after Stringer\n");
	sleep(sleeptime);
	printf("continue\n");
	// exit for fixing the stringgraph bug
	printf("CHECKPOINT: ContigWriter\n");
	if(inDB){

		printf("Re-read the input Database");
		struct reads* reads = readDB(inDB);
//		initScaff(G,reads);

		// Scaffolding
		printf("CHECKPOINT: Scaffolding\n");
		time(&start);
		initScaff(G,reads);
		readTouring(G,files,reads);
		time(&stop);
		printf("Scaffolding: %0.2f\n",difftime (stop,start));
		printf("Wait after Scaffolding\n");
		sleep(sleeptime);
		printf("continue\n");
//		exit(1);

		if(findotdump){
			sprintf(tempPath,"%s/scaffGraph.dot",pathAssembly);
			scaffGraphDot(G,reads,tempPath);
		}
		printf("Scaffolding finished, hopefully correctly ;-)\n");
		exit(1);
//		No Scaffolding
		time(&start);
		printf("CHECKPOINT: POA\n");
		struct POG* contigs_pog = make_poaScaff(G,reads,0);
		time(&stop);
		printf("POA: %0.2f\n",difftime (stop,start));
		printf("Wait after POA\n");
		sleep(sleeptime);
		printf("continue\n");
		char* contigPath = (char*)malloc(100);
		sprintf(contigPath,"%s/contigs.fasta",pathAssembly);
		printf("CHECKPOINT: FastaOut\n");
		time(&start);
		poa_printContigs(contigs_pog,contigPath);
		time(&stop);
		printf("FASTA Out: %0.2f\n",difftime (stop,start));
		printf("Wait after FastaOut\n");
		sleep(sleeptime);
		printf("continue\n");
		exit(1);
		// POG Correction (Partial-Order-Graph)
//		struct POG* contigs_pog = make_poa(G,reads);
//		struct POG* contigs_pog = make_poaScaff(G,reads,0);
//		char* contigPath = "./output/correct_contigs.fasta";
//		poa_printContigs(contigs_pog,contigPath);
//		poa_toDot("output/poa.dot");
		// Scaffolding
//		struct POG* contigs_pog = make_poaScaff(G,reads,1);
//		struct POG* contigs_pog = make_poa(G,reads);
//		char* contigPath = "./output/correct_Scaffolds.fasta";
//		poa_printContigs(contigs_pog,contigPath);
//		poa_toDot("output/poa.dot");
//		exit(1);
		if(contigs_pog){
			struct Letter_T* current;
			struct LetterEdge* edge;
			struct LetterEdge* nextedge;
			struct LetterSource_S* sedge;
			struct LetterSource_S* nextsedge;
			for(i=0;i<numNodes;i++){
//				printf("Letter NUmber %i\n",i);
				current = &Letters[i];
				edge = current->left;
				while(edge){
//					printf("Free left edge\n");
					nextedge = edge->next;
					free(edge);
					edge = nextedge;
				}
				edge = current->right;
				while(edge){
//					printf("Free right edge\n");
//					printf("Free right edge at pos: %i\n",i);
					nextedge = edge->next;
					free(edge);
					edge = nextedge;
				}
				sedge = &current->source;
				if(sedge){
					sedge = sedge->next;
					while(sedge){
//						printf("Free source\n");
						nextsedge = sedge->next;
						free(sedge);
						sedge = nextsedge;
					}
				}

			}
			for(i=0;i<contigs_pog->contigNum;i++){
				free(contigs_pog->contig[i].name);
				free(contigs_pog->contig[i].sequence);
			}

			free(Letters);
			free(contigs_pog->contig);
			free(contigs_pog);
		}


		for(i=1;i<=numreads;i++){
			free(reads[i].seq);
		}
		free(reads);
//		for(i=0;i<contigs->num;i++){
//			free(contigs->contig[i].seq);
//		}
//		free(contigs->contig);
//		free(contigs);
	}

	// Destruct all memory consuming objects

	struct bread* bread;
	struct bread* breadN;
	for(i=0;i<=G->V;i++){
		if(G->read[i]){
			bread = G->read[i]->first;
			if(bread){
//				printf("Free bread: %i of node: %i\n",bread->ID,i);
				while(bread->next){
					breadN = bread->next;
					free(bread->dest);
					free(bread);
					bread = breadN;
				}
				free(bread->dest);
				free(bread);
			}
		}
		free(G->read[i]);
	}
	free(G->read);
	free(G);

	free(S->ID);
	free(S->edge);
	free(S->length);
	free(S->side);
	free(S->status);
	free(S);

	free(inDB);

	for(i=0;i<files->libNum;i++){
		free(files[i].leftReads);
		free(files[i].rightReads);
	}
	free(files);

	struct edge* edge;
	struct edge* edgeN;
	struct ReadNode* read;
	struct ReadNode* readN;

	for(i=0;i<redGraph->V;i++){
		if(redGraph->array[i].head){
//			printf("Free Head-Edge\n");
			edge = redGraph->array[i].head;
			while(edge->next){
				edgeN = edge->next;
				free(edge);
				edge = edgeN;
			}
			free(edge);
		}
		if(redGraph->array[i].tail){
//			printf("Free Tail-Edge\n");
			edge = redGraph->array[i].tail;
			while(edge->next){
				edgeN = edge->next;
				free(edge);
				edge = edgeN;
			}
			free(edge);
		}
		if(redGraph->array[i].headread){
//			printf("Free Read-List\n");
			read = redGraph->array[i].headread;
			while(read->next){
				readN = read->next;
				free(read);
				read = readN;
			}
			free(read);
		}
	}
	printf("No program Faults occurred\nTool terminated regularly!!!\n");

	exit(1);
}
