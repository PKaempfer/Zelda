/*
 * CC.c
 *
 *  Created on: Mar 13, 2017
 *      Author: lkaempfpp
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "CC.h"
#include "ConsensusCaller.h"
#include "kmer.h"
#include "DBGraph_scaffold.h"

#define MIN_CONTIG_LEN 100
#define MIN_SCAFF_LEN 	200
#define H_RANGE 		6

static inline int max_func(int a,int b,int c,int d){
	return _max((_max(a,b)),(_max(c,d)));
}

static inline int max_func2(int a,int b,int c){
	return _max((_max(a,b)),(c));
}

static inline void poa_resetMatrix(int line, int len){
	// free the letters and set matrix to 0
	int i,j;
	static int16_t* nullLine = NULL;
	if(!nullLine){
		nullLine = (int16_t*)malloc(sizeof(int16_t)*10000);
		for(j=0;j<10000;j++){
			nullLine[j] = (j+1) * GAP_PENALTY;;
		}
	}
	if(!alMatrix_Letter) printf("No alMatrix_Letter\n");
	for(i=1;i<line;i++){
		alMatrix_Letter[i]->score = 0;
		alMatrix_Letter[i]->ml = NULL;
		if(alMatrix_Letter[i]->junction!=0){
			printf("Junction WRONGGGG in line: %i: %i\n",i,alMatrix_Letter[i]->junction);
			alMatrix_Letter[i]->junction = 0;
		}
		memcpy(&alMatrix[i][1],nullLine,sizeof(int16_t)*len);
	}
}

/**
 * Function Prints the alignment in the classical blast design
 *
 * @param readseq	Sequence of the read including gaps introduces by the back tracing through the SW Matrix
 * @param refseq	Sequence of the local part of the reference the read was globally aligned to. Including gaps introduces by the back tracing through the SW Matrix
 * @param readlen	Lenght of the alignment
 */
static inline void poa_showAlignment(char* readseq, char* refseq, int readlen){

	int i;

	printf("Alignment: \n");
	printf("Read: ");
	for(i=readlen-1;i>=0;i--){
		printf("%c",readseq[i]);
	}
	printf("\n");
	printf("      ");
	for(i=readlen-1;i>=0;i--){
		if(readseq[i] == refseq[i]){
			printf("|");
		}
		else{
			printf(" ");
		}
	}
	printf("\n");
	printf("Ref:  ");
	for(i=readlen-1;i>=0;i--){
		printf("%c",refseq[i]);
	}
	printf("\n\n");
}

static inline void resetLetters(struct Letter_T* Letters){
	int i;
	struct Letter_T* current;
	struct LetterEdge* edge;
	struct LetterEdge* nextedge;
	for(i=0;i<numNodes;i++){
		current = &Letters[i];
		edge = current->left;
		while(edge){
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
		edge = current->right;
		while(edge){
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
		current->vFlag = 0;
		current->counter = 0;
		current->align_ring = NULL;
		current->junction = 0;
		current->score = 0;
		current->left = NULL;
		current->right = NULL;
	}
}

void POG_doubletest(struct POGseq* contig){
	int i;
	int end = contig->length;
	printf(">%s\n",contig->name);
	for(i=0;i<end;i++){
		if(i%50==0 && i) printf("\n");
		if(i%50==0) printf("%i:  ",i);
		printf(" %i",Letters[i].counter);
		if(!Letters[i].right) printf("No right edge\n");

	}
}

void POG_writeBackbone(struct POGseq* contig,char* fastaFile){
	FILE* fasta = fopen(fastaFile,"a");
	char* seq = (char*)malloc(sizeof(char)*1000000);
	struct Letter_T current = Letters[contig->startLetter.dest];
	int i = 0;
	seq[i++] = current.letter;
	fprintf(fasta,">%s\n",contig->name);
	while(current.right){
		current = Letters[current.right->dest];
		seq[i++] = current.letter;
	}
	seq[i]='\0';

	int len = contig->length;
	int len2 = strlen(seq);
	printf("Length Contig: %i -> Length Sequence: %i\n",len,len2);
	for(i=0;i<len;i+=80){
		fprintf(fasta,"%.80s\n",&seq[i]);
	}

	fclose(fasta);
}

static inline void poa_LetterSizeCheck(){
	uint32_t j;
	if(numNodes == maxNumNodes){
		maxNumNodes *=2;
		struct Letter_T* temp = (struct Letter_T*)realloc(Letters,sizeof(struct Letter_T)*maxNumNodes);
		if(temp){
			printf("Letter resize\n");
			Letters = temp;
			for(j=numNodes;j<maxNumNodes;j++){
				Letters[j].left = NULL;
				Letters[j].right = NULL;
				Letters[j].junction = 0;
				Letters[j].score = 0;
				Letters[j].vFlag = 0;
			}
		}
		else{
			printf("Could not resize Letters-Array\nAbort\n");
			exit(1);
		}
	}
}

static inline void push_pogread(int readID, struct POGreadsSet* pogreadset, int seqlen,char ori){
	uint32_t preadID = pogreadset->number;
	struct POGreads* pogreads = pogreadset->pogreads;
	if(preadID == pogreadset->size){
		pogreadset->size *= 2;
		pogreadset->pogreads = (struct POGreads*)realloc(pogreadset->pogreads,(sizeof(struct POGreads))*pogreadset->size);
		if(!pogreadset->pogreads){
			printf("Error in POG-Read container reallocation\n");
			exit(1);
		}
		pogreads = pogreadset->pogreads;
	}
	pogreads[preadID].ID = readID;
	pogreads[preadID].end = numNodes-1;
	pogreads[preadID].start = numNodes - seqlen;
	pogreads[preadID].ori = ori;
	pogreadset->number++;
}

void POG_initbackbone(struct POGseq* contig, char* seq){ //struct reads* read
	int i;
	uint32_t leftID = 0;

	contig->length = strlen(seq);

	struct Letter_T* left = NULL;
	struct Letter_T* current;

	contig->startLetter.dest = numNodes;
	contig->startLetter.next = NULL;

	for(i=0;i<contig->length;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->counter = 0;
		if(left){
			current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
			current->left->counter = 0;
			current->left->dest = leftID;
			current->left->next = NULL;
			left = &Letters[leftID];
			left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
			left->right->counter = 0;
			left->right->dest = numNodes;
			left->right->next = NULL;
		}

		left = current;
		leftID = numNodes;
		numNodes++;
		poa_LetterSizeCheck();
	}
}

void POG_appendbackbone(struct POGseq* contig, char* seq, int overhang){
	int i;
	int len = strlen(seq);
	uint32_t leftID = numNodes-1;

	contig->length += overhang;

	struct Letter_T* left;
	struct Letter_T* current;

	for(i=len-overhang;i<len;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->counter = 0;
		// left-right connection
		current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		current->left->counter = 0;
		current->left->dest = leftID;
		current->left->next = NULL;
		left = &Letters[leftID];
		left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		left->right->counter = 0;
		left->right->dest = numNodes;
		left->right->next = NULL;
		left = current;
		leftID = numNodes;
		numNodes++;
		poa_LetterSizeCheck();
	}
}

struct POGreadsSet* OLC_backbone(struct POGseq* contig, struct reads* reads, struct myovlList* G, struct scaffold_set* aS, int scaffID){
    char verbose = 0;
    char verbose2 = 0;

    static struct POGreadsSet* pogreads = NULL;
    if(!pogreads){
    	pogreads = (struct POGreadsSet*)malloc(sizeof(struct POGreadsSet));
    	pogreads->size = pow(2,16);
    	pogreads->pogreads = (struct POGreads*)malloc(sizeof(struct POGreads)*pogreads->size);
    }
	pogreads->number = 0;

	struct scaffEdge* scaffEdge;
    int startJunction;

    int dir;
    int bdir;
    int finJunction;
    int breadID;
//    int oldbreadID;
    int multidir;
    char ori;
    int overhang;
    int inserts = 0;

    struct bread* bread;
    struct bread* internb;
    char* readseq = (char*)malloc(sizeof(char)*(maxReadLen+1));
    char* revreadseq = (char*)malloc(maxReadLen+1);


    scaffEdge = aS->scaff[scaffID].first;
	startJunction = aS->scaff[scaffID].startJunction;
	dir = G->read[startJunction]->dir;
	bread = G->read[startJunction]->first;
	if(aS->scaff[scaffID].type == 1){
		finJunction = aS->scaff[scaffID].endJunction;
		while(bread){
			if(verbose && bread->dest) printf("destID: %i == endJunction: %i\n",bread->dest->ID, aS->scaff[scaffID].endJunction);
			if(bread->dest && bread->dest->ID == aS->scaff[scaffID].endJunction) break;
			else bread = bread->next;
		}

	}
	else{
		while(scaffEdge){
			finJunction = scaffEdge->targetJunction;
			scaffEdge = scaffEdge->next;
			printf("FinJunction: %i\n",finJunction);
		}
		scaffEdge = aS->scaff[scaffID].first;
		while(bread){
			printf("destPathID: %i == scaffEdgeID: %i\n",bread->dest->pathID, scaffEdge->ID);
			if(bread->dest && bread->dest->pathID == scaffEdge->ID) break;
			bread = bread->next;
		}
	}

	if(bread){
		breadID = bread->ID;
		if(verbose) printf("breadID  : %i \n",breadID);
		if(dir){
			if(bread->sideflag) 	multidir = 3;
			else					multidir = 0;
		}
		else{
			if(bread->sideflag)		multidir = 1;
			else					multidir = 2;
		}
		bdir = bread->sideflag;

		// Insert First Junction
		decompressReadSt(reads[aS->scaff[scaffID].startJunction].seq,readseq,reads[aS->scaff[scaffID].startJunction].len);
		if(multidir>1){
			revReadSt(readseq,revreadseq);
			strcpy(readseq,revreadseq);
			ori = 1;
		}
		else ori = 0;
		POG_initbackbone(contig,readseq);
		push_pogread(aS->scaff[scaffID].startJunction,pogreads,reads[aS->scaff[scaffID].startJunction].len,ori);

		// Insert Contained reads in first Junction
		internb = G->read[aS->scaff[scaffID].startJunction]->first;
		while(internb){
			if(G->read[internb->ID]->flag != CONTAINED && verbose) printf("Found bread: %i (FLAG: %i) dest: %i (scaffedgeID: %i)\n",internb->ID,G->read[internb->ID]->flag,internb->dest->pathID,scaffEdge->ID);
			if(G->read[internb->ID]->flag == CONTAINED){
    			ori = 0;
    			if(multidir%2==1 && !G->read[internb->ID]->dir)	ori = 1;
    			if(multidir%2==0 && G->read[internb->ID]->dir)	ori = 1;
    			push_pogread(internb->ID, pogreads, reads[internb->ID].len,ori);
    			if(verbose) printf("ALINGING CONTAINED READ IN JUNCTION\n");
    			if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
    			if(verbose) printf("Read: %s\n",readseq);
    			inserts++;
			}
			internb = internb->next;
		}

		// Insert first proper bread of the Junction
		overhang = bread->overhang;
		decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
		ori = 0;
		if(multidir%2==1 && G->read[breadID]->dir){
			revReadSt(readseq,revreadseq);
			strcpy(readseq,revreadseq);
			ori = 1;
		}
		if(multidir%2==0 && !G->read[breadID]->dir){
			revReadSt(readseq,revreadseq);
			strcpy(readseq,revreadseq);
			ori = 1;
		}
		POG_appendbackbone(contig,readseq,overhang);
		push_pogread(breadID,pogreads,strlen(readseq),ori);

		inserts++;
		internb = bread;
		// Go on, insert all the rest
		while(breadID != finJunction){
			if(G->read[breadID]->flag == JUNCTION){
				printf("bread (%i) is a JUNCTION (finjunction: %i)\n",breadID,finJunction);
				scaffEdge = scaffEdge->next;
				if(scaffEdge->bridge){

				}
				if(!scaffEdge){
					printf("Catch: no more Edge on Junction : %i , but not endJunction reached\nAbort\n",breadID);
					exit(1);
				}
				// Insert Contained reads of the intermediate Junction
				internb = G->read[breadID]->first;
				while(internb){
					if(G->read[internb->ID]->flag != CONTAINED && verbose) printf("Found bread: %i (FLAG: %i) dest: %i (scaffedgeID: %i)\n",internb->ID,G->read[internb->ID]->flag,internb->dest->pathID,scaffEdge->ID);
					if(G->read[internb->ID]->flag == CONTAINED){
	        			ori = 0;
	        			if(multidir%2==1 && !G->read[internb->ID]->dir)	ori = 1;
	        			if(multidir%2==0 && G->read[internb->ID]->dir)	ori = 1;
	        			overhang = internb->overhang;
	        			push_pogread(internb->ID, pogreads, reads[internb->ID].len,ori);
	        			if(verbose) printf("ALINGING CONTAINED READ IN JUNCTION\n");
	        			if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
	        			if(verbose) decompressReadSt(reads[internb->ID].seq,readseq,reads[internb->ID].len);
	        			if(verbose) printf("Read: %s\n",readseq);
	        			inserts++;
					}
					internb = internb->next;
				}
				// Insert the proper read intermediate junction along the scaffold
				internb = G->read[breadID]->first;
				while(internb){
					if(internb->dest && internb->dest->pathID == scaffEdge->ID){
		    			overhang = internb->overhang;
		    			breadID = internb->ID;
	        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
	        			ori = 0;
	        			if(multidir%2==1 && G->read[breadID]->dir){
	        				revReadSt(readseq,revreadseq);
	        				strcpy(readseq,revreadseq);
	        				ori = 1;
	        			}
	        			if(multidir%2==0 && !G->read[breadID]->dir){
	            			revReadSt(readseq,revreadseq);
	            			strcpy(readseq,revreadseq);
	            			ori = 1;
	        			}
	        			POG_appendbackbone(contig,readseq,overhang);
	        			push_pogread(breadID,pogreads,strlen(readseq),ori);

						if(verbose) printf("ALINING PROPER READ\n");
						if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
	        			if(verbose) printf("Read: %s\n",readseq);
						inserts++;
						break;
					}
					internb = internb->next;
				}
			}
			else{
				// Insert the Contained reads of the proper bread
				internb = G->read[breadID]->first;
				while(internb){
					if(G->read[internb->ID]->flag == CONTAINED){
	        			ori = 0;
	        			if(multidir%2==1 && !G->read[internb->ID]->dir){
	        				revReadSt(readseq,revreadseq);
	        				strcpy(readseq,revreadseq);
	        				ori = 1;
	        			}
	        			if(multidir%2==0 && G->read[internb->ID]->dir){
	            			revReadSt(readseq,revreadseq);
	            			strcpy(readseq,revreadseq);
	            			ori = 1;
	        			}
	        			push_pogread(internb->ID,pogreads,reads[internb->ID].len,ori);
	        			if(verbose2) printf("ALIGING CONTAINED READ\n");
	        			if(verbose2) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
	        			if(verbose2) decompressReadSt(reads[internb->ID].seq,readseq,reads[internb->ID].len);
	        			if(verbose2) printf("Read: %s\n",readseq);
	        			inserts++;
					}
					internb = internb->next;
				}
				// Insert next proper read of the proper bread
				internb = G->read[breadID]->first;
				while(internb){
					if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
		    			overhang = internb->overhang;
		    			breadID = internb->ID;
	        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
	        			ori = 0;
	        			if(multidir%2==1 && G->read[breadID]->dir){
	        				revReadSt(readseq,revreadseq);
	        				strcpy(readseq,revreadseq);
	        				ori = 1;
	        			}
	        			if(multidir%2==0 && !G->read[breadID]->dir){
	            			revReadSt(readseq,revreadseq);
	            			strcpy(readseq,revreadseq);
	            			ori = 1;
	        			}
	        			POG_appendbackbone(contig,readseq,overhang);
	        			push_pogread(breadID,pogreads,strlen(readseq),ori);
						if(verbose2) printf("ALIGING PROPER READ\n");
						if(verbose2) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
	        			if(verbose2) printf("Read: %s\n",readseq);
						inserts++;
						break;
					}
					internb = internb->next;
				}
			}
		}
	}

	free(readseq);
	free(revreadseq);

	return pogreads;
}

/**
 * TODO: All Fine, May have a look at the real end node!!! Seems to cause a problem at some point.
 * 1. Alignment Function
 * Initializes the Alignment Matrix. Each Node in the part of the actual poa graph gets a matrix line. The order of the lines is then defined by the graph order.
 * Nodes in alignment rings are on the same level in the Matrix. For Detail see the POA Paper. Count Junction in the graph to not calculate paths after junctions multiple times.
 *
 * @param contig 	The Contig builds a unity of the sequence form a start to an end. Is read out by the Consensus Function.
 * @param len		The lengt of the read, that is aligned to the poa graph. This length defines the length of the matrix rows.
 * @return			Returns the number of rows in the Matrix from the beginning of the left overlapping read to the end of the right one.
 */
static inline int POG_alignPrepro(int len, uint32_t st_pos, uint32_t end_pos){
	char verbose = 0;
	char strictverbose = 0;
	if(verbose) printf("CHECKPOINT: Graph_prepro (Start: %i - End: %i)\n",st_pos,end_pos);

	struct Letter_T* current = &Letters[st_pos];
	struct Letter_T* start_node = current;

	static struct Letter_T** new_letters = NULL;
	int new_num = 0;
	static struct Letter_T** old_letters = NULL;
	int old_num = 0;
	struct Letter_T** temp_letters;

	if(!new_letters){
		new_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000); // Max breadth of graph = 100
		old_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000);
	}

	struct Letter_T* end_node = &Letters[end_pos];
	int line = 1;

	do{
		alMatrix_Letter[line] = current;
		current->ml = alMatrix[line++];
		new_letters[new_num++] = current;
		if(current->align_ring && current->align_ring != start_node){
			current = current->align_ring;
		}
		else break;
	} while(1);

//	printf("Alternative Starts: %i\n",line-1);

	struct LetterEdge* edge;
	int depth = 1;

//	if(overhang >= 30) printf("Wide range overhang could cause a problem\n");

	while(new_num && depth <= len + 200){
		if(verbose && new_num >= 90) printf("Graph breadth > 100\n");
		temp_letters = old_letters;
		old_letters = new_letters;
		new_letters = temp_letters;
		old_num = new_num;
		new_num = 0;

		while(old_num){
			old_num--;
			if(old_letters[old_num] != end_node){
				edge = old_letters[old_num]->right;
				while(edge){
					if(Letters[edge->dest].junction == 0){
						Letters[edge->dest].junction = 1;
//						if(line>240)
						if(strictverbose) printf("1 PrePro -> l: %i (d: %i) dest: %i / ori: %li\n",line,depth,edge->dest,old_letters[old_num]-&Letters[0]);
						new_letters[new_num++] = &Letters[edge->dest];
						alMatrix_Letter[line] = &Letters[edge->dest];
						Letters[edge->dest].ml = alMatrix[line++];
					}
					else{
						Letters[edge->dest].junction++;
//						if(line>240)
						if(strictverbose) printf("2 PrePro -> l: %i (d: %i) dest: %i / ori: %li\n",line,depth,edge->dest,old_letters[old_num]-&Letters[0]);
					}
//					if(insNum == 4077)	printf("Letter: %c (id: %i) num: %i\n",Letters[edge->dest].letter,edge->dest,(int)Letters[edge->dest].junction);
					edge = edge->next;
				}
			}
			else{
//				printf("EndNode found\n");
			}
		}
		// Limit number of fields in matrix to compute. Give a maximum distance from the diagonal (e.g. 5bp)
		depth++;
	}
	if(strictverbose && depth<90){
		printf("LineNumber: %i / %i\n",line,depth);
		gdepth = depth;
	}

	return line;
}

/**
 * Set first non-gap matrix line for all nodes in the alignment-ring with the contig->readleft
 * @param current 		Letter of the last included read, which is the start point for the next read alignment
 * @param new_letters	List of all possible start nodes of the alignment ring
 * @param seq			Sequence of the read which is actually aligned to the graph
 * @return				returns the number of the alignment starts; Nodes in the alignment ring
 */
static inline int POG_alignInitMatrix(struct Letter_T* current, struct Letter_T** new_letters, unsigned char* seq, char fullMatrix, int len){
	char verbose = 0;
	int j,k;
	int mat_end;
	int best_sc;
	int new_num = 0;
	struct Letter_T* start_node = current;

	int range;
	if(!fullMatrix) range = H_RANGE;
	else range = len;

	do{
		// TODO Think about possibility of non-having a right edge of the contig->readleft
		new_letters[new_num++] = current;
//		new_letters[new_num++]->junction--;

		best_sc = 0;
		if(!fullMatrix) mat_end = _min(len,range);
		else mat_end = len;
//		printf("MatEnd: %i (len: %i)\n",mat_end,len);
//		printf("curren init: %p\n",current);
		for(j=1;j<=mat_end;j++){
			// k is pos in seq, j-1, because j==0 is first gap position;
			k = j-1;
//			current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(alMatrix[0][j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(alMatrix[0][j]+GAP_PENALTY));
//			printf("%i -> %i (j: %i)/n",current->ml[j],current->ml[j-1],j);
//			current->ml[j] = _max((_max((current->ml[j-1]+GAP_PENALTY),(alMatrix[0][j-1] + SM1[codes[current->letter]][codes[seq[k]]]))),(alMatrix[0][j]+GAP_PENALTY));
			current->ml[j] = max_func2(current->ml[j-1]+GAP_PENALTY,alMatrix[0][j-1] + SM1[codes[current->letter]][codes[seq[k]]],alMatrix[0][j]+GAP_PENALTY);
//			current->ml[j] = _max((current->ml[j]),(j*GAP_PENALTY + SM1[codes[current->letter]][codes[seq[k]]]));
			if(current->ml[best_sc] < current->ml[j]) best_sc = j;

		}
		current->score = best_sc;

		if(current->align_ring && current->align_ring != start_node){
			if(verbose) printf("Found alternative start point in alignment ring:\n");
			current = current->align_ring;
		}
		else{
			break;
		}
	} while(1);

//	printf("\n");
	return new_num;
}

static inline int POG_alignFillMatrix(int* new_numG, struct Letter_T** new_letters, unsigned char* seq,struct Letter_T* end_node , struct Letter_T** end_letters, char fullMatrix, int len){
	char strictverbose = 0;
	int j,k;
	int depth = 1;
	int id;
	char rightbool;
	int best_sc;
	int mat_st,mat_end;
	int end_num = 0;

	struct LetterEdge* edge;
	struct Letter_T* left;
	struct Letter_T* current;
	struct Letter_T** temp_letters;
	int old_num = 0;
	static struct Letter_T** old_letters = NULL;
	if(!old_letters) old_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*10000);

	int range;
	int new_num = (*new_numG);
	if(!fullMatrix) range = H_RANGE;
	else range = len;

	while(new_num && depth <= len + 200){ //maxReadLen
		if(depth > len+5 && strictverbose) printf("Go deeper: %i\n",depth);
//		memcpy(old_letters,new_letters,sizeof(struct Letter_T*)*new_num);
		temp_letters = old_letters;
		old_letters = new_letters;
		new_letters = temp_letters;
		old_num = new_num;
		new_num = 0;

		while(old_num){
			old_num--;
			if(old_letters[old_num] != end_node){
				edge = old_letters[old_num]->right;
				left = old_letters[old_num];
				id = old_letters[old_num] - Letters;
				if(id < 0) id *=-1;
//				printf("Letter: %c (id: %i) num: %i\n",old_letters[old_num]->letter,id,(int)old_letters[old_num]->junction);
				while(edge){
//					new_letters[new_num++] = &Letters[edge->dest];
					if(strictverbose && gdepth < 90)
					printf("edge: %i (depth: %i)\n",edge->dest,depth);
					if(Letters[edge->dest].junction == 1){
//						if((&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1) >240)
//							printf("FILL (%lu , %p - %p) -> l: %lu (d: %i)\n",&Letters[edge->dest].ml[0] - &alMatrix[0][0],&Letters[edge->dest].ml[0],&alMatrix[0][0],(&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1),depth);
						Letters[edge->dest].junction--;
						rightbool = 1;
					}
					else if(Letters[edge->dest].junction > 1){
//						if((&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1) >240)
//							printf("FILL (%lu , %p - %p) -> l: %lu (d: %i)\n",&Letters[edge->dest].ml[0] - &alMatrix[0][0],&Letters[edge->dest].ml[0],&alMatrix[0][0],(&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1),depth);
						Letters[edge->dest].junction--;
						Letters[edge->dest].junction *= -1;
						rightbool = 3;
					}
					else if(Letters[edge->dest].junction == -1){
//						if((&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1) >240)
//							printf("FILL (%lu , %p - %p) -> l: %lu (d: %i)\n",&Letters[edge->dest].ml[0] - &alMatrix[0][0],&Letters[edge->dest].ml[0],&alMatrix[0][0],(&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1),depth);
						Letters[edge->dest].junction++;
						rightbool = 0;
					}
					else if(Letters[edge->dest].junction < -1){
//						if((&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1) >240)
//							printf("FILL (%lu , %p - %p) -> l: %lu (d: %i)\n",&Letters[edge->dest].ml[0] - &alMatrix[0][0],&Letters[edge->dest].ml[0],&alMatrix[0][0],(&Letters[edge->dest].ml[0] - &alMatrix[0][0]) / (maxReadLen+1),depth);
						Letters[edge->dest].junction++;
						rightbool = 2;
					}
					else{
						printf("This case should not happen (Junction Number < 1: %i -> (oldnum: %i / j: %i))\n",Letters[edge->dest].junction,old_num,j);
						return -1;
	//					exit(1);
					}
					if(rightbool<2){
//						printf("edge new_num: %c\n",Letters[edge->dest].letter);
						new_letters[new_num++] = &Letters[edge->dest];
					}
					current = &Letters[edge->dest];
					if(!fullMatrix){
						mat_st = _max(1,(left->score-(range-1)));
						mat_end = _min(len,(left->score+(range+1)));
					}
					else{
						mat_st = 1;
						mat_end = len;
					}

//					printf("Best Scoring field of precursor node: %i (Score: %i)\n",left->score,left->ml[left->score]);
					best_sc = 0;
					for(j=mat_st;j<=mat_end;j++){
						// k is pos in seq, j-1, because j==0 is first gap position;
						k = j-1;
						// Smith-Waterman Scoring function: Best of itself, left, diagonal, top
//						current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(left->ml[j]+GAP_PENALTY));
						if(rightbool%2==0){
							current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(left->ml[j]+GAP_PENALTY));

						}
						else{
//							if(current->ml[j] <= left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]) current->ml[j] = left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]];
//							else{
//								printf("%i > %i (at i: %i ,j: %i)\n",current->ml[j],left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]],depth,j);
//							}
//							current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(left->ml[j]+GAP_PENALTY));
							current->ml[j] = max_func2((current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(left->ml[j]+GAP_PENALTY));
						}
//						current->ml[j] = _max((_max((current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]))),(left->ml[j]+GAP_PENALTY));
//						current->ml[j] = max_func2((current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]),(left->ml[j]+GAP_PENALTY));
//						if(current->ml[j] <= left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]) current->ml[j] = left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]];
//						else{
//							printf("%i > %i\n",current->ml[j],left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]);
//						}
//						if(current->ml[j] < current->ml[j-1]+GAP_PENALTY) current->ml[j] = current->ml[j-1]+GAP_PENALTY;
//						if(current->ml[j] < left->ml[j]+GAP_PENALTY) current->ml[j] = left->ml[j]+GAP_PENALTY;
//
						if(current->ml[best_sc] < current->ml[j]) best_sc = j;
					}
					if(current->ml[current->score] < current->ml[best_sc]) current->score = best_sc;
					edge = edge->next;
				}
			}
			else{
//				printf("ENDNODE FOUND\n");
				// Critical
//				old_letters[old_num]->junction--;
//				if(old_letters[old_num]->junction == 0) end_letters[end_num++] = end_node;

				old_letters[old_num]->junction = 0;
//				if(end_node->junction) end_node->junction = 0;
				end_letters[end_num++] = end_node;
//				printf("EndNode found\n");
			}
		}
		// Limit number of fields in matrix to compute. Give a maximum distance from the diagonal (e.g. 5bp)
		depth++;
	}

	(*new_numG) = new_num;
	if(depth%2==0){
		memcpy(old_letters,new_letters,sizeof(struct Letter_T)*new_num);
		temp_letters = old_letters;
		old_letters = new_letters;
		new_letters = temp_letters;
	}
//	if(gdepth<90) exit(1);

	return end_num;

}

//static inline int poa_searchEndPoint(int line, unsigned char* seq, int insNum, char backbone, char print_align, int overhang, struct Sequence* contig, char fullMatrix,int len){
static inline int POG_alignEndpoint(int line, char fullMatrix,int len){
	char verbose = 0;
	int i;
	int best_Letter = 0;
	int best_Score = 0;
//	for(i=0;i<end_num;i++){
//		if(end_letters[i]->ml[len] > best_Score){
//			best_Letter = i;
//			best_Score = end_letters[i]->ml[len];
//		}
//	}
	for(i=1;i<line;i++){
		if(alMatrix_Letter[i]->ml[len] > best_Score){
			best_Letter = i;
			best_Score = alMatrix_Letter[i]->ml[len];
//			printf("BestLetter: %i - Best Score: %i\n",best_Letter,best_Score);
		}
	}
//	poa_showMatrix(len,line,seq);
//	printf("EndPoint: %c\n",alMatrix_Letter[best_Letter]->letter);

//	printf("Return BestLetter : %i (score: %i)\n",best_Letter,best_Score);

	if(best_Score<0){
		printf("Matrix End was not calculated\n");
		exit(1);
	}

	if(!fullMatrix && best_Score < (len*SM1[0][0])*0.98){
//		printf("BestScore: %i\n",best_Score);
		return -1;
	}
	return best_Letter;
}

static inline struct pairAlign* POG_alignBacktrack(unsigned char* seq, struct Letter_T* current,char print_Message, int j, int line, int readID){ // Parameter 4:  int readID,
	char verbose = 0;
	char verbose2 = 0;
//	static int alignmentcounter = 0;
//	alignmentcounter++;
//	if(alignmentcounter > 2262) verbose = 1;


	if(print_Message) printf("Start back tracing (j:%i)\n",j);
	char suspectVerbose = 1;
	int k;
	int seqlen = j;
	struct LetterEdge* edge;

	struct pairAlign* align = NULL;
	char* readseq = (char*)malloc(j*2);
	char* refseq = (char*)malloc(j*2);
	int len=0;
	char leftbool = 0;
	char breakF = 0;

//	struct LetterEdge* counteredge;
	struct Letter_T* ringletter;
	struct Letter_T* newLetter;
	struct Letter_T* left;
	struct Letter_T* newLetterRight = NULL;
	int32_t newLetterRightID = -1;
	struct LetterEdge* newEdge;
	struct Letter_T* current_Right = NULL;

	// Find correct letter for local end point of contained sequence alignments
//	if(!backbone){
//		int bestline=line-1;
//		for(i=line-2;i>=0;i--){
//			if(alMatrix[i][j] > alMatrix[bestline][j]) bestline = i;
//		}
//		current = alMatrix_Letter[bestline];
//	}


	int nextbool = 0;
	while(j){
		if(verbose) printf("j: %i\n",j);
		edge = current->left;
		nextbool = 0;
		if(!edge && verbose) printf("NoEdge -> Element: %lu\n",(int)(current-Letters)/sizeof(struct Letter_T));
		while(edge){
			if(Letters[edge->dest].ml){
				if(verbose) printf("edge->dest: %i (%li)\n",edge->dest,Letters[edge->dest].ml-alMatrix[0]);
				// not sure
				leftbool = 1;
				left = &Letters[edge->dest];
				k = j-1;

				if(j>0 && current->ml[j] == left->ml[j-1] + SM1[codes[current->letter]][codes[seq[k]]]){
					nextbool = 1;
					// Entry from Diagonal
					readseq[len] = (char)seq[k];
					refseq[len++] = (char)current->letter;

					// Update POG
					if(current->letter == seq[k]){
						if(verbose)	printf("Match\n");
						// Same Letter -> Match: increase counter of existing letter
						if(current->counter<255) current->counter++;
						// Position Calculations depends on endianess.
						// Check if the right letter was already a match or if it comes from a new path
						if(newLetterRight){
							if(current_Right == newLetterRight){
								// Update edge counts
								newEdge = current->right;
								while(newEdge){
									if(newLetterRightID == newEdge->dest){
										if(newEdge->counter<255) newEdge->counter++;
										break;
									}
									newEdge = newEdge->next;
								}
								newEdge = current_Right->left;
								newLetterRightID = current-Letters;
								if(newLetterRightID < 0) newLetterRightID *= -1;
								while(newEdge){
									if(newLetterRightID == newEdge->dest){
										if(newEdge->counter<255) newEdge->counter++;
										break;
									}
									newEdge = newEdge->next;
								}
								// last letter matched already
								if(verbose) printf("\t->Last letter was a ref-match\n");
							}
							else{
								// Last Letter was a mismatch and newly created
								if(verbose) printf("\t->Last letter was NO ref-match\n");
								newEdge = current->right;
								current->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
								current->right->dest = newLetterRightID;
								current->right->counter = 1;
								current->right->vFlag = 0;
								current->right->next = newEdge;
//								printf("Letters: %p current: %p\n",Letters,current);
								newLetterRightID = current - Letters;
								if(newLetterRightID < 0) newLetterRightID *= -1;
								if(verbose) printf("Old ID for next connection = %i\n",newLetterRightID);
								newEdge = newLetterRight->left;
								newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
								newLetterRight->left->dest = newLetterRightID;
								newLetterRight->left->counter = 1;
								newLetterRight->left->next = newEdge;
							}
						}
//						else{
//							printf("Diagonal and no precurser\n");
//						}
						current_Right = current;
						newLetterRight = current;
						newLetterRightID = current-Letters;
						if(newLetterRightID < 0) newLetterRightID *= -1;
						if(verbose) printf("Letters: %p current: %p structSize: %lu -> dif: %i\n",Letters,current,sizeof(struct Letter_T),current-Letters);
						if(verbose) printf("Old ID for next connection = %i\n",newLetterRightID);

					}
					else{
						// Differnt letters -> Mismatch: Create new Letter and new edges
						if(verbose) printf("Mismatch\n");
						newLetter = &Letters[numNodes];
						newLetter->counter = 1;
						newLetter->letter = seq[k];
//						newLetter->source.ipos = k;
//						newLetter->source.iseq = readID;
//						newLetter->source.next = NULL;
						// Make, Close alignment ring
						// ToDo: look if the align-ring already contains the correct letter
						ringletter = current->align_ring;
						while(ringletter && ringletter != current){
							if(ringletter->letter == seq[k]){
								printf("CASE 1: New letter is created, although the letter is already in the align ring\n");
								breakF = 1;
							}
							ringletter = ringletter->align_ring;
						}
						if(current->align_ring){
							newLetter->align_ring = current->align_ring;
							current->align_ring = newLetter;
						}
						else{
							current->align_ring = newLetter;
							newLetter->align_ring = current;
						}
						// Connect new letter and take the new latter as new current point -> Close path if next step is match
						if(newLetterRight){
							newEdge = newLetterRight->left;
							newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
							newLetterRight->left->dest = numNodes;
							newLetterRight->left->counter = 1;
							newLetterRight->left->next = newEdge;

							newLetter->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
							newLetter->right->dest = newLetterRightID;
							if(verbose) printf("dest: %i\n",newLetter->right->dest);
							newLetter->right->counter = 1;
							newLetter->right->vFlag = 0;
							newLetter->right->next = NULL;
						}
						newLetterRightID = numNodes;
						newLetterRight = newLetter;

						// ToDo: Does this make sense, or any difference at all
						current_Right = current; // ToDo: Does this make sense, or any difference at all

						numNodes++;
						if(breakF) return NULL;
						poa_LetterSizeCheck();
					}

					current = left;
					j--;
					break;
				}
				else if(j>0 && current->ml[j] == current->ml[j-1] + GAP_PENALTY){
					nextbool = 1;
					// Entry from left -> Gap in ref -> Stay in current matrix line;
					if(print_Message) printf("Gap in REF\n");
					readseq[len] = (char)seq[k];
					refseq[len++] = '-';

					// Update POG
					newLetter = &Letters[numNodes];
					newLetter->counter = 1;
					newLetter->letter = seq[k];
//					newLetter->source.ipos = k;
//					newLetter->source.iseq = readID;
//					newLetter->source.next = NULL;
					if(newLetterRight){
						newEdge = newLetterRight->left;
						newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
						newLetterRight->left->dest = numNodes;
						newLetterRight->left->counter = 1;
						newLetterRight->left->next = newEdge;

						newLetter->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
						newLetter->right->dest = newLetterRightID;
						newLetter->right->counter = 1;
						newLetter->right->vFlag = 0;
						newLetter->right->next = NULL;
					}
					newLetterRightID = numNodes;
					newLetterRight = newLetter;
					numNodes++;
					poa_LetterSizeCheck();

					j--;
					continue;
				}
				else if(current->ml[j] == left->ml[j] + GAP_PENALTY){
					nextbool = 1;
					// Entry from above -> Gap in seq
					if(print_Message) printf("Gap in Seq\n");
					readseq[len] = '-';
					refseq[len++] = (char)current->letter;

					// ToDo: Does this make sense, or any difference at all
					current_Right = current;
//					current_Right = NULL;

					current = left;
					break;
				}
				else{
					if(verbose2) printf("Noting is true, there must be an alternative path!\n");
					if(verbose2) printf("j: %i",j);
				}
			}
			edge = edge->next;
		}
		if(!leftbool) break;
		else leftbool = 0;
		if(!nextbool){
			// TODO: Circumvent in another way!!! Not just return NULL and break the contig
			if(suspectVerbose) printf("Read %i: Captured in infinite loop (j=%i), line: %lu, Abort\n",readID,j,(int)(current->ml-alMatrix[0])/((maxReadLen+1)*sizeof(int16_t)));
			POG_showMatrix(seqlen,line,seq);
//			free(refseq);
//			free(readseq);
//			return align;
			align = (struct pairAlign*)malloc(sizeof(struct pairAlign));
			readseq[len] = '\0';
			refseq[len] = '\0';
			align->len = strlen(readseq);
			printf("Alignment Length: %i\n",align->len);
			align->readSeq = (char*)malloc((align->len*2)+1);
			strcpy(align->readSeq,readseq);
			align->refSeq = (char*)malloc((align->len*2)+1);
			strcpy(align->refSeq,refseq);
			align->j = j;
			align->current = current;
			poa_showAlignment(readseq,refseq,align->len);
			free(refseq);
			free(readseq);
//			align.current = NULL;
			exit(1);
			return align;

			if(j==99) exit(1);
		}
	}

	readseq[len] = '\0';
	refseq[len] = '\0';

	align = (struct pairAlign*)malloc(sizeof(struct pairAlign));
	align->len = strlen(readseq);
	align->readSeq = (char*)malloc((align->len*2)+1);
	strcpy(align->readSeq,readseq);
	align->refSeq = (char*)malloc((align->len*2)+1);
	strcpy(align->refSeq,refseq);
	align->j = j;
	align->current = current;
	printf("J at the end of Backtracking: %i\n",j);
//	printf("Alignment Length: %i\n",align->len);
//	printf("refSeq: %s\n",refseq);
//	printf("seqSeq: %s\n",readseq);

	free(refseq);
	free(readseq);

	return align;
}

static inline char POG_alignUpdateGraph(unsigned char* seq, struct pairAlign* align, char print_Message, int len){
	// --> 6. BEGINN Connect to matrix origin
	// Connect to matrix origin
	char verbose = 0;
	int j = align->j;
	int i;
	struct Letter_T* current = align->current;
	struct Letter_T* align_ring;
	int length = align->len;
	char* readseq = align->readSeq;
	char* refseq = align->refSeq;
	int k = j-1;
	if(print_Message) printf("Current Alignment Length: %i;  pos of the seq: %i\n",length,k);

	if(print_Message){
		if(j>0) printf("Connect to matrix origin j (score: %i) (%c/%c): %i\n",current->ml[j],current->letter,seq[k],j);
	}
	if(verbose) poa_showAlignment(readseq,refseq,length);
	while(j!=0){
		k = j-1;
		if(print_Message){
			if(j>0) printf("%li -> %i = %i + %i\n",current->ml - alMatrix[0],current->ml[j], alMatrix[0][j-1],SM1[codes[current->letter]][codes[seq[k]]]);
		}
		if(j>0 && current->ml[j] == alMatrix[0][j-1] + SM1[codes[current->letter]][codes[seq[k]]]){
			// Entry from Diagonal
//			if(print_Message) printf("Origin Diagonal\n");
			readseq[length] = seq[k];
			refseq[length++] = current->letter;
			if(current->letter == seq[k]){
				current->counter++;
			}
			else{
				printf("MISMATCH Origin of the Problem????????????\n");
				align_ring = current;
				// make new letter
				// how to connect???
			}
			j--;
			for(;j>0;j--){
				k = j-1;
				readseq[length] = seq[k];
				refseq[length++] = '-';
				if(print_Message) printf("go left till origin is reached\n");
			}
			break;
		}
		else if(j>0 && current->ml[j] == current->ml[j-1] + GAP_PENALTY){
			// Entry from left -> Gap in ref -> Stay in current matrix line;
			if(print_Message) printf("Origin left\n");
			readseq[length] = seq[k];
			refseq[length++] = '-';
			j--;
			continue;
		}
		else if(current->ml[j] == alMatrix[0][j] + GAP_PENALTY){
			// Entry from above -> Gap in seq
			if(print_Message) printf("Origin top\n");
			readseq[length] = '-';
			refseq[length++] = current->letter;
			for(;j>0;j--){
				k = j-1;
				readseq[length] = seq[k];
				refseq[length++] = '-';
				if(print_Message) printf("go left till origin is reached\n");
			}
			break;
		}
		else{
			if(verbose){
				printf("%c (j: %i from %i)\n",seq[k],j,align->j);
				printf("-\t");
				for(i=0;i<=len;i++){
					printf(" %i",alMatrix[0][i]);
				}
				printf("\n");
				printf("%c\t",current->letter);
				for(i=0;i<=len;i++){
					printf(" %i",current->ml[i]);
				}
				printf("\n");
//				printf("Number of ends: %i\n",end_num);
//				printf("Depth: %i\n", depth);
				printf("Matrixline: %li\n",current->ml - alMatrix[0]);
				printf("Nothing is true: This Case should not happen\n");
				poa_showAlignment(readseq,refseq,length);
			}
			return 0;
//			exit(1);

		}
	}
	if(print_Message && j!= 0) printf("j: %i\n",j);
	return 1;
}

static inline char POG_alignUpdateGraph2(unsigned char* seq, struct pairAlign* align, char print_Message, int len){
	// --> 6. BEGINN Connect to matrix origin
	// Connect to matrix origin
	char verbose = 1;
	int j = align->j;
	int i;
	struct Letter_T* current = align->current;
	int length = align->len;
	char* readseq = align->readSeq;
	char* refseq = align->refSeq;
	int k = j-1;
	if(print_Message) printf("Current Alignment Length: %i;  pos of the seq: %i\n",length,k);

	if(print_Message){
		if(j>0) printf("Connect to matrix origin j (score: %i) (%c/%c): %i\n",current->ml[j],current->letter,seq[k],j);
	}
	while(j!=0){
		k = j-1;
		if(print_Message){
			if(j>0) printf("%li -> %i = %i + %i\n",current->ml - alMatrix[0],current->ml[j], alMatrix[0][j-1],SM1[codes[current->letter]][codes[seq[k]]]);
		}
		if(j>0 && current->ml[j] == alMatrix[0][j-1] + SM1[codes[current->letter]][codes[seq[k]]]){
			// Entry from Diagonal
//			if(print_Message) printf("Origin Diagonal\n");
			readseq[length] = seq[k];
			refseq[length++] = current->letter;
			if(current->letter == seq[k]){
				current->counter++;
			}
			else{
				// make new letter
				// how to connect???
			}
			j--;
			for(;j>0;j--){
				k = j-1;
				readseq[length] = seq[k];
				refseq[length++] = '-';
				if(print_Message) printf("go left till origin is reached\n");
			}
			break;
		}
		else if(j>0 && current->ml[j] == current->ml[j-1] + GAP_PENALTY){
			// Entry from left -> Gap in ref -> Stay in current matrix line;
			if(print_Message) printf("Origin left\n");
			readseq[length] = seq[k];
			refseq[length++] = '-';
			j--;
			continue;
		}
		else if(current->ml[j] == alMatrix[0][j] + GAP_PENALTY){
			// Entry from above -> Gap in seq
			if(print_Message) printf("Origin top\n");
			readseq[length] = '-';
			refseq[length++] = current->letter;
			for(;j>0;j--){
				k = j-1;
				readseq[length] = seq[k];
				refseq[length++] = '-';
				if(print_Message) printf("go left till origin is reached\n");
			}
			break;
		}
		else{
			if(verbose){
				printf("%c (j: %i from %i)\n",seq[k],j,align->j);
				printf("-\t");
				for(i=0;i<=len;i++){
					printf(" %i",alMatrix[0][i]);
				}
				printf("\n");
				printf("%c\t",current->letter);
				for(i=0;i<=len;i++){
					printf(" %i",current->ml[i]);
				}
				printf("\n");
//				printf("Number of ends: %i\n",end_num);
//				printf("Depth: %i\n", depth);
				printf("Matrixline: %li\n",current->ml - alMatrix[0]);
				printf("Nothing is true: This Case should not happen\n");
				poa_showAlignment(readseq,refseq,length);
			}
			return 0;
//			exit(1);

		}
	}
	if(print_Message && j!= 0) printf("j: %i\n",j);
	return 1;
}

char POG_readAlign(unsigned char* seq, int seqlen, char heuristic, uint32_t st_pos, uint32_t end_pos, int readID){
	static struct timespec alignmentSt;
	static struct timespec alignmentEnd;
	char print_Message = 0;
	char print_align = 0;
	static int line = 0;
	struct Letter_T* current;
	struct Letter_T* end_node;
	struct pairAlign* align;
	static struct Letter_T** new_letters = NULL;
	static struct Letter_T** end_letters = NULL;
	if(!new_letters){
		new_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*10000); // Max breadth of graph = 100
		end_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*10000);
	}

	char fin = 0;
	int i;
	int new_num = 0;
	int end_num = 0;
	int best_Letter;
	static int numFull = 0;
	static int numPart = 0;

	char fullMatrix = !heuristic;

	while(1){
		new_num = 0;
		end_num = 0;
#ifdef TIMEM
		clock_gettime(CLOCK_MONOTONIC, &alignmentSt);
#endif

		// 1. Alignemnt matrix size definition and line assignment to the poa Nodes
//		printf("SeqLen: %i\n",seqlen);
		poa_resetMatrix(line,seqlen);
		line = POG_alignPrepro(seqlen,st_pos,end_pos);

#ifdef TIMEM
		clock_gettime(CLOCK_MONOTONIC, &alignmentEnd);
		alignmentTime += (((alignmentEnd.tv_sec * 1000000000) + alignmentEnd.tv_nsec) - ((alignmentSt.tv_sec * 1000000000) + alignmentSt.tv_nsec));
#endif

#ifdef TIMEM
		clock_gettime(CLOCK_MONOTONIC, &ts_start);
#endif

	// 2. Initialize matrix Rows the form the first line(s) of the initial matrix before filling
		current = &Letters[st_pos];
		new_num = POG_alignInitMatrix(current, new_letters, seq,fullMatrix,seqlen);

#ifdef TIMEM
		clock_gettime(CLOCK_MONOTONIC, &ts_finish);
		sumMatrix += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));
#endif

	// 3. Fill the Alignment Matrix (Most Time Consuming Part of the entire tool)
		end_node = &Letters[end_pos];
		end_num = POG_alignFillMatrix(&new_num, new_letters, seq, end_node, end_letters, fullMatrix, seqlen);
		if(end_num <= 0){
			printf("EndNum<0 . RETURN\n");
			if(!fullMatrix){
				fullMatrix = 1;
				numFull++;
				continue;
			}
			else{
				POG_showMatrix(seqlen,line,(char*)seq);
//				printf("Ref: ");
//				for(i=0;i<line;i++){
//					printf("%c",alMatrix_Letter[i]->letter);
//				}
//				printf("seq: %s\n",(char*)seq);
//				printf("\n");
//				printf("Error in Alignment Matrix -> Abort Contig\n");
				printf("\tFILL MATRIX BREAK\n");
				return 0;
			}

		}
		// ????
//		if(new_num){
//			printf("\tAdd new Letters\n");
//			memcpy(&end_letters[end_num],new_letters,sizeof(struct Letter_T*)*new_num);
//			end_num += new_num;
//		}

	// 4. Search best Alignmet end Point
		best_Letter = POG_alignEndpoint(line,fullMatrix,seqlen);

		//	printf("BestLetter 123: %i\n",best_Letter);
//		POG_showMatrix(seqlen,line,(char*)seq);
		if(best_Letter>=0) current = alMatrix_Letter[best_Letter];
		else{
			if(!fullMatrix){
				fullMatrix = 1;
				numFull++;
				continue;
			}
			else{
				printf("\tEND MATRIX BREAK\n");
				return 0;
			}
		}

		#ifdef TIMEM
			clock_gettime(CLOCK_MONOTONIC, &ts_start);
		#endif


//		if(!align && fullMatrix) printf("\tBACKTRACK BREAK\n");
		if(fullMatrix || best_Letter>=0){
			break;
//			if(!align && !fullMatrix){
//				fullMatrix=1;
//				numFull++;
//			}
//			else{
//				numPart++;
//				break;
//			}
		}
		else{
			fullMatrix=1;
			numFull++;
//			if(verbose)
			printf("NumFull: %i (Part: %i)\n",numFull,numPart);
			if(numFull%1000 == 0){
				printf("NumFull: %i (Part: %i)\n",numFull,numPart);
			}
		}
	}

	// 5. Make backtrace
	align = POG_alignBacktrack(seq,current,print_Message,seqlen, line, readID);

	// 6. Connect to Matrix origin
	if(align){
		fin = POG_alignUpdateGraph(seq,align,print_Message,seqlen);
		if(!fin) printf("\tUPDATE GRAPH BREAK\n");
	}
	else return 0;

	char* refseq = align->refSeq;
	char* readseq = align->readSeq;

#ifdef TIMEM
	clock_gettime(CLOCK_MONOTONIC, &ts_finish);
	sumTrace += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));
#endif

	if(print_align){
		int length = align->len;
		current = align->current;
		poa_showAlignment(readseq,refseq,length);
	}
	free(refseq);
	free(readseq);
	free(align);
	return fin;
}


char POG_align(struct reads* reads, struct POGreadsSet* pogreadsSet, char heuristic, uint32_t contigLen){
	printf("CHECKPOINT: Start Alignments\n");
	char verbose = 0;
	char verbose2 = 1;
	char fin;
	int i;
	int readID;
	int readLen;
	uint32_t st_pos;
	uint32_t end_pos;
	struct POGreads* pogreads = pogreadsSet->pogreads;

	char* readseq = (char*)malloc(sizeof(char)*maxReadLen+1);
	char* revreadseq = (char*)malloc(sizeof(char)*maxReadLen+1);

	for(i=0;i<pogreadsSet->number;i++){
		readID = pogreads[i].ID;
		readLen = reads[readID].len;
		st_pos = pogreads[i].start;
		end_pos = pogreads[i].end;
		decompressReadSt(reads[readID].seq,readseq,readLen);
		if(pogreads[i].ori){
			revReadSt(readseq,revreadseq);
			strcpy(readseq,revreadseq);
		}
		if(verbose) printf("Aligning read: %i (%i -> %i)\n",i,st_pos,end_pos);
		if(st_pos < 5) st_pos = 0;
		else st_pos -= 5;
		if(end_pos+5 > contigLen) end_pos = contigLen;
		else end_pos += 5;
		fin = POG_readAlign((unsigned char*)readseq,readLen,heuristic,st_pos,end_pos,i);
		if(verbose2 && !fin){
			printf("\tAlignment Denied\n");
			return 0;
		}

	}
	free(readseq);
	free(revreadseq);
	return 1;
}


struct POG* OLC(struct myovlList* G, struct reads* reads, char scaffolding, char heuristic, struct para* para){
	char verbose = 1;

#ifdef TIMEM
	struct timespec consenusSt;
	struct timespec consenusEnd;
	static long consensusTime = 0;
#endif

	printf("Checkpoint: ConsensusCaller (CC)\n");
//	para = NULL;
	struct scaffold_set* aS;
	if(scaffolding){
		printf("Checkpoint: Init Scaffold Correction\n");
		aS = scaffold_init2();
	}
	else{
		printf("Checkpoint: Init Contig Correction\n");
		aS = contigs_init(G); // ,reads
	}

	aS = scaffold_stats(aS);


	int i,j;
    printf("MaxReadLen: %i\n",maxReadLen);
    char* name = (char*)malloc(1000);
    char* dotPath = (char*)malloc(1000);

    // Init Graph structure
    struct POGset* pog = (struct POGset*)malloc(sizeof(struct POGset));
    pog->contigNum = 0;
    pog->maxNum = 1000;
    pog->contig = (struct POGseq*)malloc(sizeof(struct POGseq)*pog->maxNum);
    if(!Letters){
    	printf("Malloc new Letters\n");
    	Letters = (struct Letter_T*)malloc(sizeof(struct Letter_T)*maxNumNodes);
        for(i=0;i<maxNumNodes;i++){
        	Letters[i].left = NULL;
        	Letters[i].right = NULL;
        	Letters[i].junction = 0;
        	Letters[i].score = 0;
        	Letters[i].vFlag = 0;
        }
    }

    if(!alMatrix_Letter){
    	printf("Init MatrixLetter\n");
    	int16_t data[maxReadLen*MATRIX_MAX_BR+1][maxReadLen+1];
    	alMatrix = (int16_t**)malloc(sizeof(int16_t*)*(maxReadLen*MATRIX_MAX_BR+1)); // Convention that the aligning part of the graph do not contain more than 5*maxReadLen nodes
    	alMatrix_Letter = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*(maxReadLen*MATRIX_MAX_BR+1));
        for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
        	alMatrix[i] = &data[i][0];
        	alMatrix_Letter[i] = NULL;
        	for(j=0;j<=maxReadLen;j++){
        		alMatrix[i][j] = j * GAP_PENALTY;
        	}
        }
    }

    struct POGreadsSet* pogreadsset;
    char dotverbose;

    for(i=0;i<aS->numbridge;i++){
    	if(aS->scaff[i].len > MIN_SCAFF_LEN || i >= aS->num){
    		if(i>=aS->num) printf("\t\tGebridgetes Scaffold (%i)\n",i);
    		pogreadsset = OLC_backbone(&pog->contig[pog->contigNum],reads,G,aS,i);
    		printf("Contig_%i:%i_%i_len:%i\n",pog->contigNum,pogreadsset->pogreads[0].ID,pogreadsset->pogreads[pogreadsset->number-1].ID,pog->contig[pog->contigNum].length);
			if(scaffolding) sprintf(name,"Scaff_%i:%i_%i_len:",pog->contigNum,pogreadsset->pogreads[0].ID,pogreadsset->pogreads[pogreadsset->number-1].ID);
			else sprintf(name,"Contig_%i:%i_%i_len:",pog->contigNum,pogreadsset->pogreads[0].ID,pogreadsset->pogreads[pogreadsset->number-1].ID);
			pog->contig[pog->contigNum].name = (char*)malloc(strlen(name)+100);
			strcpy(pog->contig[pog->contigNum].name,name);
//    		POG_writeBackbone(&pog->contig[pog->contigNum],"test.fasta");
    		dotverbose = POG_align(reads,pogreadsset, heuristic,pog->contig[pog->contigNum].length-1);

    		// debug
    		if(verbose && !dotverbose){
    			sprintf(dotPath,"%s/%s.dot",para->asemblyFolder,pog->contig[pog->contigNum].name);
    			poa_toDot(dotPath);
    		}

//    		POG_doubletest(&pog->contig[pog->contigNum]);
			if((float)numNodes/(float)aS->scaff[i].len < 2){
#ifdef TIMEM
				clock_gettime(CLOCK_MONOTONIC, &consenusSt);
#endif
				POG_alignConsensus(&pog->contig[pog->contigNum]);
#ifdef TIMEM
	    		clock_gettime(CLOCK_MONOTONIC, &consenusEnd);
	    		consensusTime += (((consenusEnd.tv_sec * 1000000000) + consenusEnd.tv_nsec) - ((consenusSt.tv_sec * 1000000000) + consenusSt.tv_nsec));
#endif
	    		if(verbose){
	    			sprintf(dotPath,"%s/%s.dot",para->asemblyFolder,pog->contig[pog->contigNum].name);
	    			poa_toDot(dotPath);
	    		}

	    		aS->scaff[i].scaffoldID = pog->contigNum;
	    		printf("i: %i , scaffID: %i\n",i,aS->scaff[i].scaffoldID);
	    		if(aS->scaff[i].next>=0){
	    			printf("Scaffold %i has a connection\n",pog->contigNum);
	    			pog->contig[pog->contigNum].seqEdge = (struct sequenceEdge*)malloc(sizeof(struct sequenceEdge));
	//    			pog->contig[pog->contigNum].seqEdge->insertLen = aS->scaff[i].next->first->bridge->estLen;
	    			pog->contig[pog->contigNum].seqEdge->insertLen = aS->scaff[aS->scaff[i].next].first->bridge->estLen;
	    			pog->contig[pog->contigNum].seqEdge->ori = 0;
	    		}
	    		else pog->contig[pog->contigNum].seqEdge = NULL;
	    		if(i >= aS->num) pog->contig[pog->contigNum].vflag = 1;
	    		else pog->contig[pog->contigNum].vflag = 0;
	    		pog->contigNum++;
			}
			else if(i>=aS->num){
				aS->scaff[i].scaffoldID = -1;
			}
    		resetLetters(Letters);
    		numNodes = 0;
    		printf("Aligned Reads: %i\n",alignedReads);
#ifdef TIMEM
    		printf("Init time:      %.3f s\n",(float)alignmentTime/1000000000);
    		printf("Alignment time: %.3f s\n",(float)sumMatrix/1000000000);
    		printf("Backtrack time: %.3f s\n",(float)sumTrace/1000000000);
    		printf("Consensus time: %.3f s\n",(float)consensusTime/1000000000);

#endif
    	}
		if(pog->contigNum == pog->maxNum){
			printf("ContigList full, Realloc memory\n");
			pog->contig = (struct POGseq*)realloc(pog->contig,sizeof(struct POGseq)*(pog->maxNum*2));
			pog->maxNum *= 2;
		}
    }
    printf("POG Finisched\n");
    verbose = 0;

    for(i=0;i<aS->numbridge;i++){
    	if(aS->scaff[i].len > MIN_SCAFF_LEN || i >= aS->num){
    		if(aS->scaff[i].next >= 0){
        		if(aS->scaff[aS->scaff[i].next].scaffoldID<0){
        			free(pog->contig[aS->scaff[i].scaffoldID].seqEdge);
        			pog->contig[aS->scaff[i].scaffoldID].seqEdge = NULL;
        			continue;
        		}
    			printf("Bild Bridge to %i\n",aS->scaff[i].next);
//    			pog->contig[aS->scaff[i].scaffoldID].seqEdge->nextScaff = aS->scaff[i].next->scaffoldID;
    			pog->contig[aS->scaff[i].scaffoldID].seqEdge->nextScaff = aS->scaff[aS->scaff[i].next].scaffoldID;
//    			printf("Connect Scaffold %i with %i bp to scaffold %i\n",aS->scaff[i].scaffoldID,pog->contig[aS->scaff[i].scaffoldID].seqEdge->insertLen,aS->scaff[i].next->scaffoldID);
    			if(verbose) printf("Connect Scaffold %i with %i bp to scaffold %i\n",aS->scaff[i].scaffoldID,pog->contig[aS->scaff[i].scaffoldID].seqEdge->insertLen,aS->scaff[aS->scaff[i].next].scaffoldID);
    		}
    	}
    }

    // _________ Scaffold Consensus Calling
    if(verbose)printf("Free ScaffoldsSet\n");
    free_schaffoldSet(aS);
    if(verbose)printf("Free ScaffoldsSeqs\n");

    struct POG* realpog = (struct POG*)malloc(sizeof(struct POG));
    realpog->contigNum = pog->contigNum;
    realpog->maxNum = pog->maxNum;
    realpog->contig = (struct Sequence*)malloc(sizeof(struct Sequence)*pog->contigNum);
    for(i=0;i<pog->contigNum;i++){
    	realpog->contig[i].name = pog->contig[i].name;
    	realpog->contig[i].length = pog->contig[i].length;
    	realpog->contig[i].seqEdge = pog->contig[i].seqEdge;
    	realpog->contig[i].startLetter = pog->contig[i].startLetter;
    	realpog->contig[i].vflag = pog->contig[i].vflag;
    	realpog->contig[i].sequence = pog->contig[i].sequence;
    	realpog->contig[i].var = pog->contig[i].var;

    }
    free(pog->contig);
    free(pog);
    free(name);
	return realpog;
}
