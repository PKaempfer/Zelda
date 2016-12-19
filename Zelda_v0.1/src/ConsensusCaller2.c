/*
 ============================================================================
 Name        : ConsensusCaller2.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Consensus Caller utilizes a re-implementation of POA - Algorithms
 	 	 	   (Lee et al.,2002) for multiple sequence alignments to create a
 	 	 	   Layout
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ConsensusCaller.h"
#include "kmer.h"


#define H_RANGE 200

/**
 * Function returns the maximum of all given integer values
 */
static inline int max_func(int a,int b,int c,int d){
	return _max((_max(a,b)),(_max(c,d)));
}

/**
 * Check size of Letters-array and resize if full
 */
static inline void poa_LetterSizeCheck(){
	uint32_t j;
	if(numNodes == maxNumNodes){
		maxNumNodes *=2;
		struct Letter_T* temp = (struct Letter_T*)realloc(Letters,sizeof(struct Letter_T)*maxNumNodes);
		if(temp){
			Letters = temp;
			for(j=numNodes;j<maxNumNodes;j++){
				Letters[j].left = NULL;
				Letters[j].right = NULL;
				Letters[j].junction = 0;
//				Letters[j].source.next = NULL;
				Letters[j].vFlag = 0;
			}
		}
		else{
			printf("Could not resize Letters-Array\nAbort\n");
			exit(1);
		}
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

/**
 * Function prints the alignment matrix according to the poa graph. Following lines could form alignment rings and may not be successors.
 * @param row		Number of the letters in the read
 * @param column	Number of the letters in the part of the reference
 * @param seq		Sequence of the read in correct orientation
 */
static inline void poa_showMatrix(int row, int column, char* seq){
	printf("CHECKPOINT: Print Matrix (row: %i, col: %i)\n",row,column);
	int i,j;
	printf("\t");
	for(i=0;i<strlen(seq);i++){
		printf("\t%c",seq[i]);
	}
	printf("\n");

	for(i=0;i<=row;i++){
		printf("\t%i",alMatrix[0][i]);
	}
	printf("\n");
	for(j=1;j<column;j++){
		printf("%c",alMatrix_Letter[j]->letter);
		for(i=0;i<=row;i++){
			printf("\t%i",alMatrix[j][i]);
		}
		printf("\n");
	}
}

/**
 * Function called for all proper reads. The overhanging bases of the next proper reads are concatenated on the end of the last proper read.
 * Then the read will be aligned to the concatenated graph
 *
 * @param contig	Meta-struct of the contig including the poa graph
 * @param G			The String graph containing all overlap information
 * @param read		Pointer to the proper concatenated read in the database including the sequence and PE information
 * @param seq		Sequence of the read already in the correct orientation
 * @param leftID	ID if the precursor read
 * @param rightID	ID of the successor read
 */
void poa_catBackbone2(struct Sequence* contig, struct myovlList *G, char* seq, int leftID, int rightID){ // Parameter 3: struct reads* read
	char verbose = 0;

	struct bread* leftB;
	struct bread* rightB;

	leftB = G->read[rightID]->first;
	while(leftB){
		if(leftB->ID == leftID) break;
		leftB = leftB->next;
	}
	rightB = G->read[leftID]->first;
	while(rightB){
		if(rightB->ID == rightID) break;
		rightB = rightB->next;
	}
	if(!leftB || !rightB){
		printf("Breads not found! Abort");
		exit(1);
	}

	int right_ovh;
//	left_ovh = leftB->overhang;
	right_ovh = rightB->overhang;
	contig->length += right_ovh;

	int len = strlen(seq);
	int i =  len - right_ovh;
	leftID = contig->readright;
	struct Letter_T* left;
	struct Letter_T* current;
	struct LetterEdge* oldEdge;

	if(verbose) printf("Cat Backbone: ");
	for(;i<len;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->counter=0;
//		// _________
//		// This things not until alignment
//		current->source.ipos = i;
//		current->source.iseq = read->ID;
//		current->source.next = NULL;
//		// _________
//		printf("NewLetter: %i (%c)\n",numNodes,current->letter);
		current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		current->left->counter = 0;
		current->left->dest = leftID;
		current->left->next = NULL;
//		printf("Set Right edge to leftID: %i\n",leftID);
		left = &Letters[leftID];
		oldEdge = left->right;
		left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		left->right->counter = 0;
		left->right->dest = numNodes;
		left->right->vFlag = 0;
		left->right->next = oldEdge;

		leftID = numNodes;
		numNodes++;
		poa_LetterSizeCheck();
		if(verbose) printf("%c",seq[i]);
	}
	if(verbose) printf("\n");
	contig->readright = numNodes-1;
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
int poa_align_prepro(struct Sequence* contig, int len, int overhang){
//	printf("CHECKPOINT: Graph_prepro\n");

	struct Letter_T* current = &Letters[contig->readleft];

	static struct Letter_T** new_letters = NULL;
	int new_num = 0;
	static struct Letter_T** old_letters = NULL;
	int old_num = 0;

	if(!new_letters){
		new_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000); // Max breadth of graph = 100
		old_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000);
	}

	struct Letter_T* end_node = &Letters[contig->readright];
	int line = 1;

	do{
		if(!current->ml){
			alMatrix_Letter[line] = current;
			current->ml = alMatrix[line++];
			new_letters[new_num++] = current;
			current->junction = 1;
		}
		else break;
		if(current->align_ring){
			current = current->align_ring;
		}
		else break;
	} while(1);

//	printf("Alternative Starts: %i\n",line-1);

	struct LetterEdge* edge;
	int depth = 0;

//	if(overhang >= 30) printf("Wide range overhang could cause a problem\n");

	while(new_num && depth <= len + overhang + 50){
//		printf("Go deeper: %i\n",depth);
		if(new_num >= 90) printf("Graph breadth > 100\n");
		memcpy(old_letters,new_letters,sizeof(struct Letter_T)*new_num);
		old_num = new_num;
		new_num = 0;

		while(old_num){
			old_num--;
			if(old_letters[old_num] != end_node){
				edge = old_letters[old_num]->right;
				while(edge){
					if(Letters[edge->dest].junction == 0){
						Letters[edge->dest].junction = 1;
						new_letters[new_num++] = &Letters[edge->dest];
						alMatrix_Letter[line] = &Letters[edge->dest];
						Letters[edge->dest].ml = alMatrix[line++];
					}
					else{
						Letters[edge->dest].junction++;
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
//	printf("LineNumber: %i\n",line);

	return line;
}


/**
 * Set first non-gap matrix line for all nodes in the alignment-ring with the contig->readleft
 * @param current 		Letter of the last included read, which is the start point for the next read alignment
 * @param new_letters	List of all possible start nodes of the alignment ring
 * @param seq			Sequence of the read which is actually aligned to the graph
 * @return				returns the number of the alignment starts; Nodes in the alignment ring
 */
int poa_initMatrix(struct Letter_T* current, struct Letter_T** new_letters, char* seq){
	char verbose = 0;
	int j,k;
	int mat_end;
	int best_sc;
	int len = strlen(seq);
	int new_num = 0;
	struct Letter_T* start_node = current;

	do{
		// TODO Think about possibility of non-having a right edge of the contig->readleft
		new_letters[new_num++] = current;

		best_sc = 0;
		mat_end = _min(len,H_RANGE);
		for(j=1;j<=mat_end;j++){
			// k is pos in seq, j-1, because j==0 is first gap position;
			k = j-1;
			current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(alMatrix[0][j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]),(alMatrix[0][j]+GAP_PENALTY));
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

	return new_num;
}

int poa_fillMatrix(int new_num, struct Letter_T** new_letters, char* seq,struct Letter_T* end_node , struct Letter_T** end_letters, int overhang){
	int j,k;
	int depth = 1;
	int id;
	char rightbool;
	int best_sc;
	int mat_st,mat_end;
	int len = strlen(seq);
	int end_num = 0;

	struct LetterEdge* edge;
	struct Letter_T* left;
	struct Letter_T* current;
	int old_num = 0;
	static struct Letter_T** old_letters = NULL;
	if(!old_letters) old_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T)*10000);

	while(new_num && depth <= maxReadLen + overhang + 50){
//		printf("Go deeper: %i\n",depth);
		memcpy(old_letters,new_letters,sizeof(struct Letter_T)*new_num);
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
				if(old_letters[old_num]->junction == 1){
					old_letters[old_num]->junction--;
					rightbool = 1;
				}
				else if(old_letters[old_num]->junction > 1){
					old_letters[old_num]->junction--;
					rightbool = 0;
				}
				else{
					printf("This case should not happen (Junction Number < 1: %i -> (oldnum: %i / j: %i))\n",old_letters[old_num]->junction,old_num,j);
					return -1;
//					exit(1);
				}

				while(edge){
//					new_letters[new_num++] = &Letters[edge->dest];
//					printf("edge: %c (depth: %i)\n",Letters[edge->dest].letter,depth);
					if(rightbool){
//						printf("edge new_num: %c\n",Letters[edge->dest].letter);
						new_letters[new_num++] = &Letters[edge->dest];
					}
					current = &Letters[edge->dest];

					mat_st = _max(1,(left->score-(H_RANGE-1)));
					mat_end = _min(len,(left->score+(H_RANGE+1)));
					best_sc = 0;
					for(j=mat_st;j<=mat_end;j++){
						// k is pos in seq, j-1, because j==0 is first gap position;
						k = j-1;
						// Smith-Waterman Scoring function: Best of itself, left, diagonal, top
						current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]),(left->ml[j]+GAP_PENALTY));
						if(current->ml[best_sc] < current->ml[j]) best_sc = j;
					}
					if(current->ml[current->score] < current->ml[best_sc]) current->score = best_sc;
					edge = edge->next;
				}
			}
			else{
				old_letters[old_num]->junction = 0;
				if(end_node->junction) end_node->junction = 0;
				end_letters[end_num++] = end_node;
//				printf("EndNode found\n");
			}
		}
		// Limit number of fields in matrix to compute. Give a maximum distance from the diagonal (e.g. 5bp)
		depth++;
	}

	return end_num;

}

int poa_searchEndPoint(int line, char* seq, int insNum, char backbone, char print_align, int overhang, struct Sequence* contig){
	int i,j;

	int len = strlen(seq);
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
		}
	}
//	poa_showMatrix(len,line,seq);
//	printf("EndPoint: %c\n",alMatrix_Letter[best_Letter]->letter);


	if(best_Score<0){
		printf("Matrix End was not calculated\n");
		exit(1);
	}

	if(best_Score < (len*SM1[0][0])*0.10){ //!backbone &&
		if(!backbone) printf("Alignment %i -> Best Score of Matrix below threshold for Containment -> Alignment denied go to next (best Score: %i/%i)  \n",insNum,best_Score,len*SM1[0][0]);
		else printf("Alignment %i -> Best Score of Matrix below threshold for PROPER -> Alignment denied go to next (best Score: %i/%i)  \n",insNum,best_Score,len*SM1[0][0]);
//		poa_part_toDot("output/error.dot",contig);
		print_align = 1;
		if(!backbone){
			for(i=1;i<line;i++){
				alMatrix_Letter[i]->ml = NULL;
				alMatrix_Letter[i]->score = 0;
				alMatrix_Letter[i]->junction = 0;
				for(j=1;j<=len;j++){
					alMatrix[i][j] = j * GAP_PENALTY;
				}
			}
			return 0;
		}
		if(!print_align){
			for(i=1;i<line;i++){
				alMatrix_Letter[i]->ml = NULL;
				alMatrix_Letter[i]->score = 0;
				alMatrix_Letter[i]->junction = 0;
				for(j=1;j<=len;j++){
					alMatrix[i][j] = j * GAP_PENALTY;
				}
			}
			return 0;
		}
		else{
			printf("CHECKPOINT: Print Matrix (row: %i, col: %i,overhang: %i)\n",len,line,overhang);
			static int part = 0;
//			poa_showMatrix(len,line,seq);
			char* dotFile = (char*)malloc(100);
			sprintf(dotFile,"part_%i.dot",part);
			printf("POA: %s\n",dotFile);
			printf("Sequence: %s\n",seq);
			poa_part_toDot(dotFile,contig);
			part++;
		}
	}
	return best_Letter;
}

struct pairAlign poa_backtrace(struct Sequence* contig, char* seq, struct Letter_T* current,char print_Message, char backbone){ // Parameter 4:  int readID,
	if(print_Message) printf("Start back tracing\n");
	int j = strlen(seq);
	int k;
	struct LetterEdge* edge;

	struct pairAlign align;
	char* readseq = (char*)malloc(j*2);
	char* refseq = (char*)malloc(j*2);
	int len=0;
	char leftbool = 0;

//	struct LetterEdge* counteredge;
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
//		printf("j: %i\n",j);
		edge = current->left;
		nextbool = 0;
		while(edge){
			if(Letters[edge->dest].ml){
//				printf("edge->dest: %i (%i)\n",edge->dest,Letters[edge->dest].ml-alMatrix[0]);
				// not sure
				leftbool = 1;
				left = &Letters[edge->dest];
				k = j-1;

				if(j>0 && current->ml[j] == left->ml[j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]){
					nextbool = 1;
					// Entry from Diagonal
					readseq[len] = seq[k];
					refseq[len++] = current->letter;

					// Update POG
					if(current->letter == seq[k]){
//						printf("Match\n");
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
//								printf("\t->Last letter was a ref-match\n");
							}
							else{
								// Last Letter was a mismatch and newly created
//								printf("\t->Last letter was NO ref-match\n");
								newEdge = current->right;
								current->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
								current->right->dest = newLetterRightID;
								current->right->counter = 1;
								current->right->vFlag = 0;
								current->right->next = newEdge;
//								printf("Letters: %p current: %p\n",Letters,current);
								newLetterRightID = current - Letters;
								if(newLetterRightID < 0) newLetterRightID *= -1;
//								printf("Old ID for next connection = %i\n",newLetterRightID);
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
//						printf("Letters: %p current: %p structSize: %i -> dif: %i\n",Letters,current,sizeof(struct Letter_T),current-Letters);
//						printf("Old ID for next connection = %i\n",newLetterRightID);

					}
					else{
						// Differnt letters -> Mismatch: Create new Letter and new edges
//						printf("Mismatch\n");
						newLetter = &Letters[numNodes];
						newLetter->counter = 1;
						newLetter->letter = seq[k];
//						newLetter->source.ipos = k;
//						newLetter->source.iseq = readID;
//						newLetter->source.next = NULL;
						// Make, Close alignment ring
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
//							printf("dest: %i\n",newLetter->right->dest);
							newLetter->right->counter = 1;
							newLetter->right->vFlag = 0;
							newLetter->right->next = NULL;
						}
						newLetterRightID = numNodes;
						newLetterRight = newLetter;
						numNodes++;
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
					readseq[len] = seq[k];
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
					refseq[len++] = current->letter;
					current_Right = NULL;
					current = left;
					break;
				}
				else{
//					printf("Noting is true, there must be an alternative path!\n");
//					printf("j: %i",j);
				}
			}
			edge = edge->next;
		}
		if(!leftbool) break;
		else leftbool = 0;
		if(!nextbool){
			// TODO: Circumvent in another way!!! Not just return NULL and break the contig
			printf("Captured in infinite loop (j=%i), Abort\n",j);
			free(refseq);
			free(readseq);
			align.current = NULL;
			return align;
			exit(1);
			if(j==99) exit(1);
		}
	}

	if(backbone){
		if(print_Message)printf("New LeftStart: %i (%c)\n",newLetterRightID,Letters[newLetterRightID].letter);
		contig->readleft = newLetterRightID;
	}

	readseq[len] = '\0';
	refseq[len] = '\0';

	align.len = strlen(readseq);
	align.readSeq = (char*)malloc(strlen(seq)*2+1);
	strcpy(align.readSeq,readseq);
	align.refSeq = (char*)malloc(strlen(seq)*2+1);
	strcpy(align.refSeq,refseq);
	align.j = j;
	align.current = current;

	free(refseq);
	free(readseq);

	return align;
}

void poa_updateGraph(char* seq, struct pairAlign* align, char print_Message){
	// --> 6. BEGINN Connect to matrix origin
	// Connect to matrix origin
	int j = align->j;
	int i;
	struct Letter_T* current = align->current;
	int length = align->len;
	char* readseq = align->readSeq;
	char* refseq = align->refSeq;
	int len = strlen(seq);
	int k = j-1;
	if(print_Message) printf("Current Alignment Length: %i;  pos of the seq: %i\n",length,k);

	if(print_Message){
		if(j>0) printf("Connect to matrix origin j (score: %i) (%c/%c): %i\n",current->ml[j],current->letter,seq[k],j);
	}
	while(j!=0){
		k = j-1;
		if(print_Message){
			if(j>0) printf("%li -> %i = %i + %i\n",current->ml - alMatrix[0],current->ml[j], alMatrix[0][j-1],SM1[codes[(int)current->letter]][codes[(int)seq[k]]]);
		}
		if(j>0 && current->ml[j] == alMatrix[0][j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]){
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
			printf("%c (j: %i)\n",seq[k],j);
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
//			printf("Number of ends: %i\n",end_num);
//			printf("Depth: %i\n", depth);
			printf("Matrixline: %li\n",alMatrix[0]-current->ml);
			printf("Nothing is true: This Case should not happen\n");
			poa_showAlignment(readseq,refseq,length);
			exit(1);

		}
	}
	if(print_Message && j!= 0) printf("j: %i\n",j);
}

static inline void poa_resetMatrix(int line, int len){
	// free the letters and set matrix to 0
	int i,j;
	for(i=1;i<line;i++){
		alMatrix_Letter[i]->ml = NULL;
		alMatrix_Letter[i]->score = 0;
		alMatrix_Letter[i]->junction = 0;
		for(j=1;j<=len;j++){
			alMatrix[i][j] = j * GAP_PENALTY;
		}
	}
}

void poa_handleErrors(){

}

/**
 * Aligns a read sequence to the part of the poa graph, provided by the overlap information of the string graph.
 * If read is proper it is part of the backbone the new part of the PO-graph is set to this area.
 * @param contig 	Is the PO-graph of this contig
 * @param read		Is the read to align with the part of the PO-graph
 * @param seq		Is the read sequence to align
 * @param backbone	Is a boolean value if the read was proper, than it is set as new reference point for the area in the PO-graph for the next read alignment
 */
char poa_heuristic_align2(struct Sequence* contig, struct reads* read, char* seq, char backbone, int insNum, int overhang){
	static char print_align = 0;
	static char print_Message = 0;
	static struct Letter_T** new_letters = NULL;
	int new_num = 0;
	static struct Letter_T** end_letters = NULL;
	int end_num = 0;

	if(!new_letters){
		new_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*10000); // Max breadth of graph = 100
		end_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*10000);
	}

	if(print_Message){
		if(backbone) printf("CHECKPOINT: Start PO_Alignment for PROPER Read\n");
		else printf("CHECKPOINT: Start PO_Alignment for CONTAINED Read\n");
	}

//	int i;
//	int j;

	// 1. Alignemnt matrix size definition and line assignment to the poa Nodes
	int line = poa_align_prepro(contig,strlen(seq),overhang);
	if(print_Message) printf("Build SW-Matrix of read: %i\n",read->ID);


	// 2. Initialize matrix Rows the form the first line(s) of the initial matrix before filling
	struct Letter_T* current = &Letters[contig->readleft];
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	new_num = poa_initMatrix(current, new_letters, seq);

//	if(print_Message) printf("Number of alternative start positions given by an alignment ring: %i\n",new_num);
//	if(print_Message) printf("Number of nodes used in partial graph: %i\n",line);

	// 3. Fill the Alignment Matrix
	struct Letter_T* end_node = &Letters[contig->readright];
	end_num = poa_fillMatrix(new_num,new_letters,seq,end_node,end_letters,overhang);

	if(end_num < 0){
		return 0;
	}

	clock_gettime(CLOCK_MONOTONIC, &ts_finish);
	sumMatrix += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));

	if(new_num){
		memcpy(&end_letters[end_num],new_letters,sizeof(struct Letter_T)*new_num);
		end_num += new_num;
	}
	if(print_Message) printf("Number of alternative ends: %i\n",end_num);

	// 4. Search best Alignmet end Point
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	int best_Letter = poa_searchEndPoint(line,seq,insNum,backbone,print_align,overhang,contig);
	if(!best_Letter) return 0;
	else current = alMatrix_Letter[best_Letter];
	int ID = current-Letters;
	if(ID < 0) ID *= -1;

	// 5. Make backtrace
	int len = strlen(seq);
	struct pairAlign align = poa_backtrace(contig,seq,current,print_Message,backbone); // Parameter 4: read->ID,
	if(!align.current){
		poa_resetMatrix(line,len);
		return 0;
	}

	// 6. Connect to Matrix origin
	poa_updateGraph(seq,&align,print_Message);
	char* refseq = align.refSeq;
	char* readseq = align.readSeq;
	int length = align.len;
//	j = align.j;
	current = align.current;

	clock_gettime(CLOCK_MONOTONIC, &ts_finish);
	sumTrace += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));

	if(backbone){
		contig->readright = ID;
		if(print_Message) printf("Set new end for cat the new backbone: %c (%i)\n",Letters[best_Letter].letter,best_Letter);
	}

	// Show Alignment
	if(print_Message) printf("Alignment Done\n");
//	print_align = 1;
	if(print_align){
		poa_showAlignment(readseq,refseq,length);
//		char* out = (char*)malloc(100);
//		sprintf(out,"output/error/Aln_%i.dot",insNum);
//		poa_part_toDot(out,contig);
//		free(out);
	}

	poa_resetMatrix(line,len);

	free(readseq);
	free(refseq);
	return 1;
}
