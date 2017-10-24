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
#include "CC.h"
#include "read_filter.h"

int main(int argc, char* argv[]) {

	char scaffolding = 0;
	char heuristic = 1;
	char prefilter;
	char findotdump = 0;
	char dotdump = 0;
	char sleeptime = 0;
	char prefilteronly = 0;
	time_t start,stop;

	struct para* para = readCMDline(argc, argv);
	prefilter = para->prefilter;

	if(para->run == 1 || para->run == 3){
		printf("\nStep 0: Craete DB\n");
		printf("###################################################\n");
		time(&start);
		makeDB(para->readDB, para->blocks, para->files);
		time(&stop);
		printf("Time: %0.2f\n",difftime (stop,start));
		if(sleeptime){
			sleep(sleeptime);
			printf("Continue\n");
		}
		freeFiles(para);
		if(para->run == 1) finished(para);
	}
	printf("\n");



	char* tempPath = (char*)malloc(1000);
	nK = para->kSize;
	minovlLen = para->minOvlLen;

	const int NUM_THREADS = para->threads;
    pthread_t threads[NUM_THREADS];


    if(prefilter == 1){
		// preliminary filter
//    	nK = 13;
	    printf("Step 1: Create Graph\n");
	    printf("###################################################\n");
		time(&stop);
		para->files = fileScheduler_DB(para->readDB,NUM_THREADS,threads);
		struct reads* reads1 = readDB(para->readDB);
		filter_reads(reads1,NUM_THREADS,threads);
		strcat(para->readDB,"filter");
		write_filteredDB(para->readDB,para->blocks,para->files,reads1);
		if(prefilteronly){
			write_filteredFasta(para->files,reads1);
			exit(1);
		}
		freeEnds_oa();
		freeHashTable_oa();
		printf("CHECKPOINT: Free Reads\n");
		freeDB(reads1);
		freeFiles(para);
		nK = para->kSize;
    }
    else if(prefilter == 2){
    	// Hashing
	    printf("Step 1: Create Graph\n");
	    printf("###################################################\n");
    	printf("CHECKPOINT: Start New Hashing\n");
    	time(&stop);
    	strcat(para->readDB,"filter");
    }

    printf("CHECKPOINT: Build DeBruijn Graph\n");
	para->files = fileScheduler_DB(para->readDB,NUM_THREADS,threads);

	createGraph(graphSize);
	time(&start);
	printf("Time: %0.2f\n",difftime (start,stop));
	if(sleeptime){
		sleep(sleeptime);
		printf("Continue\n");
	}
	printf("\n");


	printf("Step 2: Fill Adjacency List\n");
	printf("###################################################\n");
	hashToTabDFS_oa();
	time(&stop);
	printf("Time: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after Hashing\n");
		sleep(sleeptime);
		printf("Continue\n");
	}
	printf("\n");

//	printf("Step 3: Free HashTable\n");
	freeHashTable_oa();
//	if(sleeptime){
//		printf("Wait after Freeing HashTbale\n");
//		sleep(sleeptime);
//		printf("Continue\n");
//	}
//	printf("\n");

	printf("Step 3: Error Correction\n");
	printf("###################################################\n");
	time(&start);
	perfectErrorCorrection();// errorCorrection();
//	writeDot("./output/test_aftererror.dot");
	time(&stop);
	printf("Time: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after Error correction\n");
		sleep(sleeptime);
		printf("continue\n");
	}
	printf("\n");


	printf("Step 4: Graph Reduction\n");
	printf("###################################################\n");
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
		reduceRedGraph();
	}while(verticalReduction());
	if(findotdump){
		printf("CHECKPOINT: Write dot-File: Reduced_DBG.dot\n");
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
//	if(findotdump){
//		sprintf(tempPath,"%s/Reduced_DBG_strong.dot",para->asemblyFolder);
//		printRedDot(tempPath);
//	}
	time(&stop);
	printf("Reducer Time: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after Reduction\n");
		sleep(sleeptime);
		printf("continue\n");
	}

	countRemainingNodes();
	redGraphConnector();
	if(dotdump){
		sprintf(tempPath,"%s/redGraphAfter.list",para->asemblyFolder);
		printRedGraphToFile(tempPath);
	}

	printf("Step 5: Stringer (Overlaps to String Graph):\n");
	printf("###################################################\n");
	time(&start);
	printf("CHECKPOINT: Re-read the input Database\n");
	struct reads* reads = readDB(para->readDB);
	struct myovlList* G = initOVLgraph(numreads);
	struct string_graph* S = initStringGraph(G,para->asemblyFolder,findotdump,reads);
	time(&stop);
	printf("Stringer: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after Stringer\n");
		sleep(sleeptime);
	}
	printf("\n");

	// Scaffolding
	printf("Step 6: Scaffolding\n");
	printf("###################################################\n");
	time(&start);
	initScaff(G,reads);
	readTouring(G,para->files,reads);

	findotdump = 1;

	contig_repeatFinder();
	balancePaths(G,reads);
//	hierholzerTourAll(G,reads);

	if(findotdump){
		sprintf(tempPath,"%s/scaffGraph.dot",para->asemblyFolder);
		scaffGraphDot(G,reads,tempPath);
	}
	time(&stop);
	printf("Scaffolding: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after Scaffolding\n");
		sleep(sleeptime);
		printf("continue\n");
	}
	printf("\n");

	printf("Step 7: POA (Layout-Consensus)\n");
	printf("###################################################\n");
	time(&start);
	struct POG* contigs_pog;
	contigs_pog = OLC(G,reads,scaffolding,heuristic, para);
	time(&stop);
	printf("POA: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after POA\n");
		sleep(sleeptime);
		printf("continue\n");
	}
	printf("\n");

	printf("Step 8: Write GenomeAssembly Files\n");
	printf("###################################################\n");
	if(scaffolding){
//		exit(1);
//		contigs_pog = make_poaScaff(G,reads,1,para,heuristic);
		char* contigPath = (char*)malloc(100);
		sprintf(contigPath,"%s/scaff.fasta",para->asemblyFolder);
		printf("CHECKPOINT: Write Fasta File\n");
		time(&start);
		poa_printContigs(contigs_pog,contigPath);
		sprintf(tempPath,"%s/scaff.vcf",para->asemblyFolder);
		printf("CHECKPOINT: Write VCF File\n");
		poa_reportVariant(contigs_pog,tempPath,contigPath);
		free(contigPath);
	}
	else{
		char* contigPath = (char*)malloc(100);
		sprintf(contigPath,"%s/contigs.fasta",para->asemblyFolder);
		printf("CHECKPOINT: Write Fasta File\n");
		time(&start);
		poa_printContigs(contigs_pog,contigPath);
		sprintf(tempPath,"%s/contigs.vcf",para->asemblyFolder);
		printf("CHECKPOINT: Variation Calling\n");
		poa_reportVariant(contigs_pog,tempPath,contigPath);
		free(contigPath);
	}
	time(&stop);
	printf("FASTA Out: %0.2f\n",difftime (stop,start));
	if(sleeptime){
		printf("Wait after FastaOut\n");
		sleep(sleeptime);
		printf("continue\n");
	}
	printf("\n");

	poa_deleteVariant(contigs_pog);
	if(contigs_pog) free_POG(contigs_pog);
	freeDB(reads);
	freeMyOvlList(G,S);
	free(tempPath);
	finished(para);
	free(Letters);
}
