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
	time_t start,stop;
	int i;

	struct para* para = readCMDline(argc, argv);

	if(para->run == 1 || para->run == 3){
		makeDB(para->readDB, para->blocks, para->files);
		if(para->run == 1) finished(para);
	}

	char* tempPath = (char*)malloc(1000);
	nK = para->kSize;
	minovlLen = para->minOvlLen;

	const int NUM_THREADS = para->threads;
    pthread_t threads[NUM_THREADS];
	para->files = fileScheduler_DB(para->readDB,NUM_THREADS,threads);

	printf("CHECKPOINT: 1. create Graph\n");
	time(&stop);
	createGraph(graphSize);
	time(&start);
	printf("Time: %0.2f\n",difftime (start,stop));

	printf("CHECKPOINT: 2. Fill Adjacency List\n");
	hashToTabDFS_oa();
	time(&stop);
	printf("Time: %0.2f\n",difftime (stop,start));
	printf("Wait after Hashing\n");
	sleep(sleeptime);
	printf("Continue\n");

	printf("CHECKPOINT: 3. HashTable Destructor\n");
	freeHashTable_oa();
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
		sprintf(tempPath,"%s/Reduced_DBG.dot",para->asemblyFolder);
		printRedDot(tempPath);
	}
	countRemainingNodes();
//	printRedGraph();
	printf("CHECKPOINT: Reduce Graph (strong)\n");
	if(dotdump){
		sprintf(tempPath,"%s/redGraphBefore.list",para->asemblyFolder);
		printRedGraphToFile(tempPath);
	}
	do{
		printf("Graph reduction\n");
		reduceRedGraph_strong();
	}while(verticalReduction());
	if(findotdump){
		sprintf(tempPath,"%s/Reduced_DBG_strong.dot",para->asemblyFolder);
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
		sprintf(tempPath,"%s/redGraphAfter.list",para->asemblyFolder);
		printRedGraphToFile(tempPath);
	}

//	testFunct();

//	printRedGraph();
//	printRedGraphToFile("./output/redGraph.list");
//	exit(1);
	printf("CHECKPOINT: Stringer:\n");
	time(&start);
	printf("Numread: %i\n",numreads);
	struct myovlList* G = initOVLgraph(numreads);
	struct string_graph* S = initStringGraph(G,para->asemblyFolder);
//	exit(1);
	time(&stop);
	printf("Stringer: %0.2f\n",difftime (stop,start));
	printf("Wait after Stringer\n");
	sleep(sleeptime);
	printf("continue\n");
	// exit for fixing the stringgraph bug
	printf("CHECKPOINT: ContigWriter\n");


	printf("Re-read the input Database");
	struct reads* reads = readDB(para->readDB);
//		initScaff(G,reads);

	// Scaffolding
	printf("CHECKPOINT: Scaffolding\n");
	time(&start);
	initScaff(G,reads);
	readTouring(G,para->files,reads);
	time(&stop);
	printf("Scaffolding: %0.2f\n",difftime (stop,start));
	printf("Wait after Scaffolding\n");
	sleep(sleeptime);
	printf("continue\n");
//		exit(1);

	if(findotdump){
		sprintf(tempPath,"%s/scaffGraph.dot",para->asemblyFolder);
		scaffGraphDot(G,reads,tempPath);
	}
	printf("Scaffolding finished, hopefully correctly ;-)\n");
//		exit(1);
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
	sprintf(contigPath,"%s/contigs.fasta",para->asemblyFolder);
	printf("CHECKPOINT: FastaOut\n");
	time(&start);
	poa_printContigs(contigs_pog,contigPath);
	time(&stop);
	printf("FASTA Out: %0.2f\n",difftime (stop,start));
	printf("Wait after FastaOut\n");
	sleep(sleeptime);
	printf("continue\n");
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
		free_POG(contigs_pog);
	}
	finished(para);

	for(i=1;i<=numreads;i++){
		free(reads[i].seq);
	}
	free(reads);
//		for(i=0;i<contigs->num;i++){
//			free(contigs->contig[i].seq);
//		}
//		free(contigs->contig);
//		free(contigs);


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
