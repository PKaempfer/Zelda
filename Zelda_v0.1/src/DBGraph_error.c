/*
 ============================================================================
 Name        : DBGraph_error.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Error correction of the de Bruijn Graph by zipping paths
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "FileReader.h"
#include "kmer.h"
#include "DBGraph.h"
#include "DBGraph_error.h"

#define IterError
#define counterCorrection

int counterTrash = 45;
int maxDir;

void rekDot(int a, int layer, FILE *part){
	struct AdjListNode *edge;
	setVFlag(a);
	if(graph->array[a].head || graph->array[a].tail){
		fprintf(part,"%i [label=\"%i\"];\n",a,a);
		edge=graph->array[a].head;
		while(edge){
			fprintf(part,"%i -> %i [label=%c];\n",edge->dest,a,toCharTrans(edge->trans));
			if(layer<maxDir && !isVFlagged(edge->dest)) rekDot(edge->dest,layer+1,part);
			edge = edge->next;
		}
		edge = graph->array[a].tail;
		while(edge){
			fprintf(part,"%i -> %i [label=%c];\n",a,edge->dest,toCharTrans(edge->trans)+32);
			if(layer>-maxDir && !isVFlagged(edge->dest)) rekDot(edge->dest,layer-1,part);
			edge = edge->next;
		}
	}
}

void writePartialDot(int ori, int width){
	FILE *part = fopen("Part.dot","w");
	int layer = 0;
	maxDir = width;
	fprintf(part,"digraph dbg {\n");

	rekDot(ori,layer,part);

	fprintf(part,"}");
	fclose(part);
}

void cleanAllFlaggs(){
	int i;
	for(i=1;i<graph->V;i++){
		delVFlag(i);
	}
}

void doulbeEdgeTest(int a, int phase,int aa, int b, int ori){
	struct AdjListNode *s,*sNext;
	s = graph->array[a].tail;
	while(s){
		sNext = s->next;
		while(sNext){
			if(s->dest == sNext->dest && toCharTrans(s->trans) == toCharTrans(sNext->trans)){
				s = graph->array[a].tail;
				while(s){
					printf("-> %i (%c) --> a: %i b: %i ori: %i\n",s->dest,toCharTrans(s->trans),aa,b,ori);
					s = s->next;
				}
				printf("Double Tail Edge after Phase: %i\n",phase);
				exit(0);
			}
			sNext = sNext->next;
		}
		s = s->next;
	}

}

long deletedKmer = 0;


/**
 * Collapse Sibling Nodes a and b. a-node is destination and inherits all non-duplicate information from b
 */
void collapseNodes(int a, int b){
	char verbose = 0;
	deletedKmer++;
	struct AdjListNode *anode,*bnode,*bpar,*bchi;
	int c, delete;
	if(verbose) printf("a: %i-count: %i / b: %i-count: %i\n",a,graph->array[a].counter,b,graph->array[b].counter);
	if(graph->array[a].counter+graph->array[b].counter<256) graph->array[a].counter += graph->array[b].counter;
	else graph->array[a].counter = 255;

	// 0. take on the read links
//	printf("Check 0 a: %i b:%i ori: %i\n",a,b,ori);
	struct LinkListNode *linksA;
	if(graph->array[a].Link){
		linksA = graph->array[a].Link;
		while(linksA->next){
			linksA = linksA->next;
		}
		linksA->next = graph->array[b].Link;
		graph->array[b].Link = NULL;
	}
	else{
		graph->array[a].Link = graph->array[b].Link;
		graph->array[b].Link = NULL;
	}

	// 1. childeren of B's parents
//	printf("Check 1\n");
	bpar = graph->array[b].head;
	while(bpar){
		c = bpar->dest;
		delete = 0;
		bchi = graph->array[c].tail;
		anode = bchi;
//		doulbeEdgeTest(c,0,a,b,ori);
		while(bchi){
			if(bchi->dest == b){
				bnode = graph->array[c].tail;
				while(bnode){
//					printf("Problem: %i %i %i\n",bnode->dest,(bnode->trans & TRANS_MASK),(bchi->trans & TRANS_MASK));
					if(bnode->dest == a && (bnode->trans & TRANS_MASK) == (bchi->trans & TRANS_MASK)){
						if(bchi == graph->array[c].tail){
							graph->array[c].tail = bchi->next;
							free(bchi);
							bchi = graph->array[c].tail;
							delete = 1;
							break;
						}
						else{
							anode->next=bchi->next;
							free(bchi);
							bchi = anode->next;
							delete = 1;
							break;
						}
					}
					bnode = bnode->next;
				}
				if(!delete){
					bchi->dest = a;
					anode = bchi;
					bchi = bchi->next;
				}
				else delete = 0;
			}
			else{
				anode = bchi;
				bchi = anode->next;
			}
		}
//		doulbeEdgeTest(c,1,a,b,ori);
		bpar = bpar->next;
	}

	// 2. B's parents itself
	//	printf("Check 2\n");
	bnode = graph->array[b].head;
	while(bnode){
		delete = 0;
		anode = graph->array[a].head;
		while(anode){
			if(anode->dest == bnode->dest && (anode->trans & TRANS_MASK) == (bnode->trans & TRANS_MASK)){
				delete = 1;
			}
			anode = anode->next;
		}
		if(delete){
			graph->array[b].head = bnode->next;
			free(bnode);
			bnode = graph->array[b].head;
		}
		else{
			graph->array[b].head = bnode->next;
			bnode->next = graph->array[a].head;
			graph->array[a].head = bnode;
		}
		bnode = graph->array[b].head;
	}


	// 3. parents of B's childeren
	//	printf("Check 3\n");
	bpar = graph->array[b].tail;
	while(bpar){
		c = bpar->dest;
		delete = 0;
		bchi = graph->array[c].head;
		anode = bchi;
		while(bchi){
			if(bchi->dest == b){
				bnode = graph->array[c].head;
				while(bnode){
					if(bnode->dest == a && (bnode->trans & TRANS_MASK) == (bchi->trans & TRANS_MASK)){
						if(bchi == graph->array[c].head){
							graph->array[c].head = bchi->next;
							free(bchi);
							bchi = graph->array[c].head;
							delete = 1;
							break;
						}
						else{
							anode->next=bchi->next;
							free(bchi);
							bchi = anode->next;
							delete = 1;
							break;
						}
					}
					bnode = bnode->next;
				}
				if(!delete){
					bchi->dest = a;
					anode = bchi;
					bchi = bchi->next;
				}
				else delete = 0;
			}
			else{
				anode = bchi;
				bchi = anode->next;
			}
		}
		bpar = bpar->next;
	}

	// 4. B's children itself
	//	printf("Check 4\n");
	bnode = graph->array[b].tail;
	while(bnode){
		delete = 0;
		anode = graph->array[a].tail;
		while(anode){
			if(anode->dest == bnode->dest && (anode->trans & TRANS_MASK) == (bnode->trans & TRANS_MASK)){
				delete = 1;
			}
			anode = anode->next;
		}
		if(delete){
			graph->array[b].tail = bnode->next;
			free(bnode);
			bnode = graph->array[b].tail;
		}
		else{
			graph->array[b].tail = bnode->next;
			bnode->next = graph->array[a].tail;
			graph->array[a].tail = bnode;
		}
		bnode = graph->array[b].tail;
	}
//	doulbeEdgeTest(a,4,a,b,ori);
//	printf("Correction end\n");
//	printGraph();
}


void rekCorrection(int a, int b, int up){
//	printf("Rek Collapse: a: %i / b:%i\n",a,b);
	static uint32_t corrCount = 0;
	corrCount++;
	collapseNodes(a,b);
	if(corrCount%100000==0) printf("Number of current corrections: %i\n",corrCount);
	static int oldA, temp;
	oldA = a;
	static struct AdjListNode *node, *nextnode;
	if(up){
		if(graph->array[oldA].head && graph->array[oldA].head->next){
			node=graph->array[oldA].head;
			while(node){
				nextnode=node->next;
				while(nextnode){
					if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
						a = _min(node->dest,nextnode->dest);
						b = _max(node->dest,nextnode->dest);
#ifdef counterCorrection
						if(graph->array[a].counter < graph->array[b].counter){
							temp = a;
							a = b;
							b = temp;
						}
#endif
						if(graph->array[a].counter > counterTrash && graph->array[b].counter > counterTrash){
							nextnode = nextnode->next;
							continue;
						}
						else if(oldA == a || oldA == b){
//							return;
							nextnode = nextnode->next;
							continue;
						}
						else if(isChild(oldA,a) || isChild(oldA,b)){
//							printf("A oldA,a: %i,(%i,%i)\n",oldA,a,b);
//							return;
							nextnode = nextnode->next;
							continue;
						}
						else rekCorrection(a,b,1);
						return;
					}
					nextnode = nextnode->next;
				}
				node = node->next;
			}
		}
	}
	else{
		if(graph->array[oldA].tail && graph->array[oldA].tail->next){
			node=graph->array[oldA].tail;
			while(node){
				nextnode=node->next;
				while(nextnode){
					if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
						a = _min(node->dest,nextnode->dest);
						b = _max(node->dest,nextnode->dest);
#ifdef counterCorrection
						if(graph->array[a].counter < graph->array[b].counter){
							temp = a;
							a = b;
							b = temp;
						}
#endif
						if(graph->array[a].counter > counterTrash && graph->array[b].counter > counterTrash){
							nextnode = nextnode->next;
							continue;
						}
						else if(oldA == a || oldA == b) {
//							return;
							nextnode = nextnode->next;
							continue;
						}
						else if(isParent(oldA,a) || isParent(oldA,b)){
//							printf("B oldA,a: %i,(%i,%i)\n",oldA,a,b);
//							return;
							nextnode = nextnode->next;
							continue;
						}
						else rekCorrection(a,b,0);
						return;
					}
					nextnode = nextnode->next;
				}
				node = node->next;
			}
		}
	}
}

/** Recursive Error Correction Function - Not time efficient */
//static inline int collapse(int ori, int dir){
////	static int number = 0;
//	static int a,b, temp;
//	struct AdjListNode *node, *nextnode;
//	if(dir){
////		printf("while0: %i\n",ori);
//		if(graph->array[ori].head && graph->array[ori].head->next){
//			node = graph->array[ori].head;
//// 			nextnode =  graph->array[ori].head->next; // ->
//			while(node){
////				printf("while 1: %i\n",ori);
//				nextnode = node->next;
//				while(nextnode){
//					if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
//						// Control the coverage of the nodes (Needs to be implemented)
//						// Previously it has to be implemented in hash to adjacency transformation functions
//						a = _min(node->dest,nextnode->dest);
//						b = _max(node->dest,nextnode->dest);
//#ifdef counterCorrection
//						if(graph->array[a].counter < graph->array[b].counter){
//							temp = a;
//							a = b;
//							b = temp;
//						}
//#endif
//						if(graph->array[a].counter > counterTrash && graph->array[b].counter > counterTrash){
//							nextnode = nextnode->next;
//							continue;
//						}
//						else if(ori == a || ori == b){
//							nextnode = nextnode->next;
//							continue;
//						}
//						else if(isChild(ori,a) || isChild(ori,b)){
//							nextnode = nextnode->next;
//							continue;
//						}
//						else rekCorrection(a,b,1);
//						return 1;
//					}
//					nextnode = nextnode->next;
//				}
//				node = node->next;
////				printf("Nextnode\n");
//			}
//		}
//	}
//	else{
//		if(graph->array[ori].tail && graph->array[ori].tail->next){
//			node = graph->array[ori].tail;
////			nextnode =  graph->array[ori].tail->next; // ->
//			while(node){
//				nextnode = node->next;
//				while(nextnode){
//					if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
//						// Control the coverage of the nodes (Needs to be implemented)
//						// Previously it has to be implemented in hash to adjacency transformation functions
//						a = _min(node->dest,nextnode->dest);
//						b = _max(node->dest,nextnode->dest);
//#ifdef counterCorrection
//						if(graph->array[a].counter < graph->array[b].counter){
//							temp = a;
//							a = b;
//							b = temp;
//						}
//#endif
//						if(graph->array[a].counter > counterTrash && graph->array[b].counter > counterTrash){
//							nextnode = nextnode->next;
//							continue;
//						}
////						printf("Collapse Paths starting: ori: %i (a:%i / b%i)\n",ori,a,b);
//						else if(ori == a || ori == b){
//							nextnode = nextnode->next;
//							continue;
////							return 0;
//						}
//						else if(isParent(ori,a) || isParent(ori,b)){
//							nextnode = nextnode->next;
//							continue;
////							return 0;
//						}
//						else rekCorrection(a,b,0);
//						return 1;
//					}
//					nextnode = nextnode->next;
//				}
//				node = node->next;
////				printf("Nextnode\n");
//			}
//		}
//	}
//	return 0;
//}

struct timespec collapseSt;
struct timespec collapseEnd;
struct timespec iterSt;
struct timespec iterEnd;
long sumcollapse = 0;
long sumIter = 0;
long collapseCountUp = 0;
long collapseCountDown = 0;

static inline int iter_collapse(int ori, int dir){
//	static int number = 0;
	clock_gettime(CLOCK_MONOTONIC, &collapseSt);
	static int a,b, temp;
	struct AdjListNode *node, *nextnode;
	char found = 1;
	char foundany = 0;
	static uint32_t numNeighbor;
	if(dir){
		while(found){
			found = 0;
			numNeighbor = 0;
			node = graph->array[ori].head;
			if(node && node->next){
				while(node){
					numNeighbor++;
					if(graph->array[node->dest].counter <= counterTrash){
						found = 1;
//						break;
					}
					node = node->next;
				}
			}
			if(found){
//				if(numNeighbor>10) printf("High Number of parents: %i\n",numNeighbor);
				found = 0;
				node = graph->array[ori].head;
				while(node){
					nextnode = node->next;
					while(nextnode){
						if(graph->array[node->dest].counter > counterTrash && graph->array[nextnode->dest].counter > counterTrash){
							nextnode = nextnode->next;
							continue;
						}
						if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
							// Control the coverage of the nodes (Needs to be implemented)
							// Previously it has to be implemented in hash to adjacency transformation functions
							a = _min(node->dest,nextnode->dest);
							b = _max(node->dest,nextnode->dest);
#ifdef counterCorrection
							if(graph->array[a].counter > graph->array[b].counter){
								temp = a;
								a = b;
								b = temp;
							}
#endif
							if(ori == a || ori == b){
								nextnode = nextnode->next;
								continue;
							}
							else if(isChild(ori,a) || isChild(ori,b)){
								nextnode = nextnode->next;
								continue;
							}
							else{
								found=1;
								foundany = 1;
								clock_gettime(CLOCK_MONOTONIC, &iterSt);
								collapseNodes(a,b);
								clock_gettime(CLOCK_MONOTONIC, &iterEnd);
								sumIter += (((iterEnd.tv_sec * 1000000000) + iterEnd.tv_nsec) - ((iterSt.tv_sec * 1000000000) + iterSt.tv_nsec));
								collapseCountUp++;
								ori = a;
								break;
							}
						}
						nextnode = nextnode->next;
					}
					if(found == 1) break;
					else node = node->next;
				}
			}
		}
	}
	else{
		while(found){
			found = 0;
			numNeighbor = 0;
			node = graph->array[ori].tail;
			if(node && node->next){
				while(node){
					numNeighbor++;
					if(graph->array[node->dest].counter <= counterTrash){
						found = 1;
//						break;
					}
					node = node->next;
				}
			}
			if(found){
//				if(numNeighbor>10) printf("High Number of children: %i\n",numNeighbor);
				found = 0;
				node = graph->array[ori].tail;
	//			nextnode =  graph->array[ori].tail->next; // ->
				while(node){
					nextnode = node->next;
					while(nextnode){
						if(graph->array[node->dest].counter > counterTrash && graph->array[nextnode->dest].counter > counterTrash){
							nextnode = nextnode->next;
							continue;
						}
						if((node->trans & TRANS_MASK) == (nextnode->trans & TRANS_MASK) && node->dest != nextnode->dest){
							// Control the coverage of the nodes (Needs to be implemented)
							// Previously it has to be implemented in hash to adjacency transformation functions
							a = _min(node->dest,nextnode->dest);
							b = _max(node->dest,nextnode->dest);
#ifdef counterCorrection
							if(graph->array[a].counter > graph->array[b].counter){
								temp = a;
								a = b;
								b = temp;
							}
#endif
	//						printf("Collapse Paths starting: ori: %i (a:%i / b%i)\n",ori,a,b);
							if(ori == a || ori == b){
								nextnode = nextnode->next;
								continue;
	//							return 0;
							}
							else if(isParent(ori,a) || isParent(ori,b)){
								nextnode = nextnode->next;
								continue;
	//							return 0;
							}
							else{
								found=1;
								foundany = 1;
								clock_gettime(CLOCK_MONOTONIC, &iterSt);
								collapseNodes(a,b);
								clock_gettime(CLOCK_MONOTONIC, &iterEnd);
								sumIter += (((iterEnd.tv_sec * 1000000000) + iterEnd.tv_nsec) - ((iterSt.tv_sec * 1000000000) + iterSt.tv_nsec));
								collapseCountDown++;
								ori = a;
								break;
							}
						}
						nextnode = nextnode->next;
					}
					if(found == 1) break;
					else node = node->next;
				}
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &collapseEnd);
	sumcollapse += (((collapseEnd.tv_sec * 1000000000) + collapseEnd.tv_nsec) - ((collapseSt.tv_sec * 1000000000) + collapseSt.tv_nsec));
	if(foundany) return 1;
	else return 0;
}

void indelHandle(int ori){
	struct AdjListNode *node, *nextnode;
	node = graph->array[ori].head;
	while(node && node->dest == ori){
		printf("Edge pointing to edge origin (InDel Artifact) ->  Delete Edge\n");
		nextnode = node->next;
		free(node);
		graph->array[ori].head = nextnode;
		node = graph->array[ori].head;

	}
	while(node && node->next){
		nextnode = node->next;
		if(nextnode->dest == ori){
			printf("Edge pointing to edge origin (InDel Artifact) ->  Delete Edge\n");
			node->next = nextnode->next;
			free(nextnode);
		}
		node = node->next;
	}
	node = graph->array[ori].tail;
	while(node && node->dest == ori){
		printf("Edge pointing to edge origin (InDel Artifact) ->  Delete Edge\n");
		nextnode = node->next;
		free(node);
		graph->array[ori].tail = nextnode;
		node = graph->array[ori].tail;
	}
	while(node && node->next){
		nextnode = node->next;
		if(nextnode->dest == ori){
			printf("Edge pointing to edge origin (InDel Artifact) ->  Delete Edge\n");
			node->next = nextnode->next;
			free(nextnode);
		}
		node = node->next;
	}
}

void perfectErrorCorrection(){
//	struct AdjListNode *nextnode;
//	int a, b;
	char verbose = 0;
	int i;
	int round=0;
	int change = 0;
	int oldchange;
	int pro;
	char *dotFile = (char*)malloc(sizeof(char)*100);
	do{
		pro = 0;
		printf("ErrorCorrection Round %i\n",round);
		countKmers();
		change = 0;
		for(i=1;i<graph->V;i++){
			if(i%(graph->V/10)==0){
				pro++;
				printf("%.1f %% (%i) nodes Corrected\n",(float)pro*10,i);
				printf("Number of collapsed Kmers: up: %li / down: %li\n",collapseCountUp,collapseCountDown);
	    		printf("Iter time:      %.3f s\n",(float)sumcollapse/1000000000);
	    		printf("Collapse time:  %.3f s\n",(float)sumIter/1000000000);
			}
			oldchange = change;
#ifdef IterError
			while(iter_collapse(i,1)){
#else
			while(collapse(i,1)){
#endif
//				printf("Colappse: %i!\n",i);
				change++;
//				if(round>19) printSorrounding(i);
				if(change-oldchange>10) break;
			}
#ifdef IterError
			while(iter_collapse(i,0)){
#else
			while(collapse(i,0)){
#endif
//				printf("Colappse!Rev: %i\n",i);
				change++;
				if(change-oldchange>10) break;
			}
		}
		round++;
		if(round>5) break;
		if(verbose){
			sprintf(dotFile,"test_%i.dot",round);
			writeDot(dotFile);
		}
		printf("Deleted k-mers: %li\n",deletedKmer);
	} while(change);
	// Delete self pointing circles ?!?
	indelHandle(i);
}

