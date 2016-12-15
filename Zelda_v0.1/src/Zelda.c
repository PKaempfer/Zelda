/*
 ============================================================================
 Name        : Zelda.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
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

	struct para* para = readCMDline(argc, argv);

	if(para->run == 1 || para->run == 3){
		printf("CHECKPOINT: 0. craete DB\n");
		time(&start);
		makeDB(para->readDB, para->blocks, para->files);
		time(&stop);
		printf("Time: %0.2f\n",difftime (stop,start));
		sleep(sleeptime);
		printf("Continue\n");
		freeFiles(para);
		if(para->run == 1) finished(para);
	}

	char* tempPath = (char*)malloc(1000);
	nK = para->kSize;
	minovlLen = para->minOvlLen;

	const int NUM_THREADS = para->threads;
    pthread_t threads[NUM_THREADS];
    printf("CHECKPOINT: 1. create Graph\n");
	time(&stop);
    para->files = fileScheduler_DB(para->readDB,NUM_THREADS,threads);
	createGraph(graphSize);
	time(&start);
	printf("Time: %0.2f\n",difftime (start,stop));
	sleep(sleeptime);
	printf("Continue\n");

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

	printf("CHECKPOINT: 4 Error Correction\n");
	time(&start);
	perfectErrorCorrection();// errorCorrection();
//	writeDot("./output/test_aftererror.dot");
	time(&stop);
	printf("Error Correction Time: %0.2f\n",difftime (stop,start));
	printf("Wait after Error correction\n");
	sleep(sleeptime);
	printf("continue\n");

	printf("CHECKPOINT: 5. Graph Reduction\n");
	time(&start);
	prepareGraph();
	printf("CHECKPOINT: Reduce Graph (light)\n");
//	printGraph();
	travToRedOVL_v2();
	freeGraph();
	labelGraph();
//	countRemainingNodes();
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
//		printf("Horizontal reduction\n");
		reduceRedGraph();
//		printf("Vertical reduction\n");
	}while(verticalReduction());
	printf("Write dot-File: Reduced_DBG.dot");
	if(findotdump){
		sprintf(tempPath,"%s/Reduced_DBG.dot",para->asemblyFolder);
		printRedDot(tempPath);
	}
//	printRedGraph();
	countRemainingNodes();

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

	printf("CHECKPOINT: 6. Stringer (Overlaps to String Graph):\n");
	time(&start);
	struct myovlList* G = initOVLgraph(numreads);
	struct string_graph* S = initStringGraph(G,para->asemblyFolder,findotdump);
//	exit(1);
	time(&stop);
	printf("Stringer: %0.2f\n",difftime (stop,start));
	printf("Wait after Stringer\n");
	sleep(sleeptime);
	printf("CHECKPOINT: ContigWriter\n");


	// Scaffolding
	printf("CHECKPOINT: 7. Scaffolding\n");
	printf("Re-read the input Database\n");
	struct reads* reads = readDB(para->readDB);
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
//		exit(1);
//		No Scaffolding
	time(&start);

	printf("CHECKPOINT: 8. POA (Layout-Consensus)\n");
	struct POG* contigs_pog = make_poaScaff(G,reads,0,para); // 3. Argument, scaffolding 1 - yes, 0 - no
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
	sprintf(tempPath,"%s/contigs.vcf",para->asemblyFolder);
	printf("CHECKPOINT: Variation Calling\n");
	poa_reportVariant(contigs_pog,tempPath,contigPath);
	time(&stop);
	printf("FASTA Out: %0.2f\n",difftime (stop,start));
	printf("Wait after FastaOut\n");
	sleep(sleeptime);
	printf("continue\n");
	poa_deleteVariant(contigs_pog);
	if(contigs_pog) free_POG(contigs_pog);
	freeDB(reads);
	freeMyOvlList(G,S);
	free(contigPath);
	free(tempPath);
	finished(para);
}
