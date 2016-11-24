/*
 ============================================================================
 Name        : DBGraph_stringer_2.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Heart of Zelda: Converts the dBG to a String Graph by calculating
 	 	 	   transitively irreducible overlaps from the dBG.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DBGraph.h"
#include "DBGraph_stringer.h"
//#include "DBGraph_reduced.h"

#define TESTNODE 1051
#define TESTNODE2 1957
#define TESTREAD 403606
#define MINOVL 31

//#define UPSTREAMTEST

//#define CONTAIN_MASSAGE

// Functions for Stringer implementation based on the ideas of Gene Myers.

static char setaRead2(int rNid, int dir, struct myovlList *ovlGraph){
	if(!ovlGraph->read[rNid]){
		ovlGraph->read[rNid] = (struct aread*)malloc(sizeof(struct aread));
		ovlGraph->read[rNid]->dir = dir;
		ovlGraph->read[rNid]->first = NULL;
		ovlGraph->read[rNid]->length = readLenList[rNid];
		ovlGraph->read[rNid]->flag = PROPER;
//		printf("Set A-Read: dir = %i\n",dir);
		return 1;
	}
	else if(ovlGraph->read[rNid]->flag == PROPER){
		ovlGraph->read[rNid]->dir = dir;
		return 1;
	}
	else return 0;
}

static void setOVL(int Nid, int overhang,int ID, int ovlflag,struct myovlList *ovlGraph,char side){
	struct bread *newovl;
	newovl = (struct bread*)malloc(sizeof(struct bread));
	newovl->next = ovlGraph->read[Nid]->first;
	ovlGraph->read[Nid]->first = newovl;
	newovl->ID = ID;
	newovl->overhang = overhang;
	newovl->flag = ovlflag;
	newovl->sideflag = side;
	newovl->dest = NULL;
}

static char setbRead2(int rNid, int rnNid, char dir, struct myovlList *ovlGraph, int abover, int baover, int ovlflag){
	static int verbose = 0;
	static int verbose2 = 0;

	// a is proper and b is already contained
	if(!setaRead2(rnNid,dir,ovlGraph)) return 1;

	struct bread *newovl;
	// Test if the b-read was already set
	newovl = ovlGraph->read[rNid]->first;
	while(newovl){
		if(newovl->ID == rnNid) return 2;
		newovl = newovl->next;
	}

	// a is proper and b was or is proper now
	if(abover > 0 && baover > 0){
		setOVL(rnNid,baover,rNid,ovlflag+4,ovlGraph,1);
		setOVL(rNid,abover,rnNid,ovlflag,ovlGraph,0);
		if(verbose){
			if(rNid == 6216){
				printf("\n\n\n\n\n\n");
				printf("Set Proper Overlaps : %i --> %i",rNid,rnNid);
				newovl = ovlGraph->read[rNid]->first;
				while(newovl){
					printf("Overlap to read: %i\n",newovl->ID);
					newovl = newovl->next;
				}
				printf("\n\n\n\n\n\n");
			}
		}
		return 2;
	}
	else if(abover <= 0){
		if(verbose2) printf("Contained OVL %i  to  %i (ab: %i / ba: %i) dir1: %i dir2: %i\n",rnNid,rNid,abover,baover,ovlGraph->read[rNid]->dir,ovlGraph->read[rnNid]->dir);
		setOVL(rNid,abover,rnNid,ovlflag,ovlGraph,0);
		setOVL(rnNid,baover,rNid,ovlflag+4,ovlGraph,1);
		ovlGraph->read[rnNid]->flag = CONTAINED;
		//__________
		if(ovlGraph->read[rNid]->dir == ovlGraph->read[rnNid]->dir){
			ovlGraph->read[rnNid]->first->dir = 0;
		}
		else{
			ovlGraph->read[rnNid]->first->dir = 1;
		}
		//__________
		newovl = ovlGraph->read[rnNid]->first->next;
		while(newovl){
			ovlGraph->read[rnNid]->first->next = newovl->next;
			free(newovl);
			newovl = ovlGraph->read[rnNid]->first->next;

		}
#ifdef CONTAIN_MASSAGE
		printf("CONTAINMENT: %i is in %i (ab: %i -> ba: %i)\n",rnNid,rNid,abover,baover);
#endif
		return 1;
	}
	if(baover == 0 && abover > 0){
		if(verbose2) printf("Contained OVL %i  to  %i (ab: %i / ba: %i) dir1: %i dir2: %i\n",rNid,rnNid,abover,baover,ovlGraph->read[rNid]->dir,ovlGraph->read[rnNid]->dir);
		setOVL(rNid,abover,rnNid,ovlflag,ovlGraph,0);
		setOVL(rnNid,baover,rNid,ovlflag+4,ovlGraph,1);
		ovlGraph->read[rNid]->flag = CONTAINED;
		//__________
		if(ovlGraph->read[rNid]->dir == ovlGraph->read[rnNid]->dir){
			ovlGraph->read[rNid]->first->dir = 0;
		}
		else{
			ovlGraph->read[rNid]->first->dir = 1;
		}
		//__________
		newovl = ovlGraph->read[rNid]->first->next;
		while(newovl){
			ovlGraph->read[rNid]->first->next = newovl->next;
			free(newovl);
			newovl = ovlGraph->read[rNid]->first->next;

		}
		return 0;
	}

	printf("Hier geht was krachtig schief!!!!!!!!! -> abover: %i baover: %i (case: %i)\n",abover,baover,ovlflag);
	return 2;
}

void downstreamTree(int i, struct partGraph* childReadIds){
	//	redGraph->vFlag[i] = 0;
		childReadIds->array[0].ID = i;
		childReadIds->array[0].flag = 0;
		childReadIds->array[0].depth = 0;
		childReadIds->array[0].layer = 0;
		int j=0,k=1,l=0;
		int end = 0;
		struct edge *child;
		struct edge *newEdge;
		int jNode;

		while(childReadIds->array[j].tail){
			newEdge = childReadIds->array[j].tail;
			childReadIds->array[j].tail = childReadIds->array[j].tail->next;
			free(newEdge);
		}

		do{
			jNode = childReadIds->array[j].ID;
			child = redGraph->array[jNode].tail;
	//		printf("jNode: %i\n",jNode);
			while(child && childReadIds->array[j].depth < MAX_DEPTH){
	//			printf("Children: -> %i\n",child->dest);
				if(redGraph->vFlag[child->dest] != 1){
	//				printf("If Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
					redGraph->vFlag[child->dest] = 1;
					childReadIds->array[k].ID = child->dest;
					childReadIds->array[k].layer = childReadIds->array[j].layer+1;
					childReadIds->array[k].depth = redGraph->array[child->dest].len + childReadIds->array[j].depth;

					if(childReadIds->array[k].head){
						while(childReadIds->array[k].head->next){
							newEdge = childReadIds->array[k].head->next->next;
							free(childReadIds->array[k].head->next);
							childReadIds->array[k].head->next = newEdge;
						}
						childReadIds->array[k].head->dest = j;
					}
					else{
						newEdge = (struct edge*)malloc(sizeof(struct edge));
						newEdge->next = childReadIds->array[k].head;
						newEdge->dest = j;
						childReadIds->array[k].head = newEdge;
					}
					while(childReadIds->array[k].tail){
						newEdge = childReadIds->array[k].tail;
						childReadIds->array[k].tail = childReadIds->array[k].tail->next;
						free(newEdge);
					}
					newEdge = (struct edge*)malloc(sizeof(struct edge));
					newEdge->next = childReadIds->array[j].tail;
					newEdge->dest = k;
					childReadIds->array[j].tail = newEdge;
	//				printf("Edge depth: %i -> %i\n",childReadIds->array[k].ID,childReadIds->array[k].depth);
					k++;
				}
				else{
	//				printf("Else Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
					for(l=0;l<k;l++){
						if(childReadIds->array[l].ID == child->dest){
							newEdge = (struct edge*)malloc(sizeof(struct edge));
							newEdge->dest = l;
							newEdge->next = childReadIds->array[j].tail;
							childReadIds->array[j].tail = newEdge;

							newEdge = (struct edge*)malloc(sizeof(struct edge));
							newEdge->dest = j;
							newEdge->next = childReadIds->array[l].head;
							childReadIds->array[l].head = newEdge;
						}
					}
				}
				child = child->next;
			}
			j++;
			if(k>40 || j==k) end = 1;
		} while(!end);

		childReadIds->V = k;
	//	printf("TempGraph for node: %i (V: %i)\n",i,k);
	//	int layer = 0;
		for(l=0;l<k;l++){
	//		if(childReadIds->array[l].layer != layer){
	//			layer++;
	////			 printf("   ");
	//		}
	//		 printf("%i ",childReadIds->array[l].ID);
			redGraph->vFlag[childReadIds->array[l].ID] = 0;
		}
	//	 printf("\n");

	//	struct edge *edge;
	//	for(j=0;j<childReadIds->V;j++){
	//		printf("Test for right edges for node: %i",childReadIds->array[j].ID);
	//		printf("\nHead: ");
	//		edge = childReadIds->array[j].head;
	//		while(edge){
	//			printf("->%i ",childReadIds->array[edge->dest].ID);
	//			edge = edge->next;
	//		}
	//		printf("\nTail: ");
	//		edge = childReadIds->array[j].tail;
	//		while(edge){
	//			printf("->%i ",childReadIds->array[edge->dest].ID);
	//			edge = edge->next;
	//		}
	//		printf("\n");
	//	}

}


void upstreamTree(int i, struct partGraph* parentReadIds){
	parentReadIds->array[0].ID = i;
	parentReadIds->array[0].flag = 0;
	parentReadIds->array[0].depth = 0;
	parentReadIds->array[0].layer = 0;
	int j=0,k=1,l=0;
	int end = 0;
	struct edge *child;
	struct edge *newEdge;
	int jNode;

	while(parentReadIds->array[j].tail){
		newEdge = parentReadIds->array[j].tail;
		parentReadIds->array[j].tail = parentReadIds->array[j].tail->next;
		free(newEdge);
	}

	do{
		jNode = parentReadIds->array[j].ID;
		child = redGraph->array[jNode].head;
//		printf("jNode: %i\n",jNode);
		while(child && parentReadIds->array[j].depth < MAX_DEPTH){
//			printf("1 Parent: -> %i (flag: %i)\n",child->dest,redGraph->vFlag[child->dest]);
			if(redGraph->vFlag[child->dest] != 1){
//				printf("If Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
				redGraph->vFlag[child->dest] = 1;
//				redGraph->vFlag[j] = -1;
				parentReadIds->array[k].ID = child->dest;
				parentReadIds->array[k].layer = parentReadIds->array[j].layer+1;
				parentReadIds->array[k].depth = redGraph->array[child->dest].len + parentReadIds->array[j].depth;

				if(parentReadIds->array[k].head){
					while(parentReadIds->array[k].head->next){
						newEdge = parentReadIds->array[k].head->next->next;
						free(parentReadIds->array[k].head->next);
						parentReadIds->array[k].head->next = newEdge;
					}
					parentReadIds->array[k].head->dest = j;
				}
				else{
					newEdge = (struct edge*)malloc(sizeof(struct edge));
					newEdge->next = parentReadIds->array[k].head;
					newEdge->dest = j;
					parentReadIds->array[k].head = newEdge;
				}
				while(parentReadIds->array[k].tail){
					newEdge = parentReadIds->array[k].tail;
					parentReadIds->array[k].tail = parentReadIds->array[k].tail->next;
					free(newEdge);
				}
				newEdge = (struct edge*)malloc(sizeof(struct edge));
				newEdge->next = parentReadIds->array[j].tail;
				newEdge->dest = k;
				parentReadIds->array[j].tail = newEdge;
//				printf("Edge depth: %i -> %i\n",parentReadIds->array[k].ID,parentReadIds->array[k].depth);
				k++;
			}
//			else{
////				printf("Else Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
////				printf("2 Parent: -> %i\n",child->dest);
//				for(l=0;l<k;l++){
//					if(parentReadIds->array[l].ID == child->dest){
//						newEdge = (struct edge*)malloc(sizeof(struct edge));
//						newEdge->dest = l;
//						newEdge->next = parentReadIds->array[j].tail;
//						parentReadIds->array[j].tail = newEdge;
//
//						newEdge = (struct edge*)malloc(sizeof(struct edge));
//						newEdge->dest = j;
//						newEdge->next = parentReadIds->array[l].head;
//						parentReadIds->array[l].head = newEdge;
//					}
//				}
//			}
			child = child->next;
		}
		j++;
		if(k>40 || j==k) end = 1;
	} while(!end);

	parentReadIds->V = k;
//	printf("TempGraph for node: %i (V: %i)\n",i,k);
//	int layer = 0;
	for(l=0;l<k;l++){
//		if(parentReadIds->array[l].layer != layer){
//			layer++;
//			printf("   ");
//		}
//		printf("%i(%i) ",parentReadIds->array[l].ID,parentReadIds->array[l].depth);
		// Very Important, do not commend out
		redGraph->vFlag[parentReadIds->array[l].ID] = 0;
	}
//	 printf("\n");

//	struct edge *edge;
//	for(j=0;j<parentReadIds->V;j++){
//		printf("Test for right edges for node: %i",parentReadIds->array[j].ID);
//		printf("\nHead: ");
//		edge = parentReadIds->array[j].head;
//		while(edge){
//			printf("->%i ",parentReadIds->array[edge->dest].ID);
//			edge = edge->next;
//		}
//		printf("\nTail: ");
//		edge = parentReadIds->array[j].tail;
//		while(edge){
//			printf("->%i ",parentReadIds->array[edge->dest].ID);
//			edge = edge->next;
//		}
//		printf("\n");
//	}
}

void free_UpstreamTree(struct partGraph* parentReadIds){
	int i;
	struct edge* edge;
	for(i=0;i<1000;i++){
		while(parentReadIds->array[i].head){
			edge = parentReadIds->array[i].head->next;
			free(parentReadIds->array[i].head);
			parentReadIds->array[i].head = edge;
		}
		while(parentReadIds->array[i].tail){
			edge = parentReadIds->array[i].tail->next;
			free(parentReadIds->array[i].tail);
			parentReadIds->array[i].tail = edge;
		}
	}
	free(parentReadIds->array);
	free(parentReadIds);
}

// End start is again mandatory

void printTag(struct partGraph *parentReadIds,struct ReadNode* rNode){
	static struct ReadNode* tagReadNode;
	int i;

	printf("rNode: %i\n",rNode->read->ID);
	printf("Tagged Path:\n");
	for(i = 0 ; i < parentReadIds->V; i++){
		if(parentReadIds->array[i].flag){
			printf("Taglist of Vertex: %i\n",parentReadIds->array[i].ID);
			tagReadNode = redGraph->array[parentReadIds->array[i].ID].headread;
			while(tagReadNode){
				if(tagReadNode->flag){
					printf(" -> %i (%i / %i)\n",tagReadNode->read->ID,tagReadNode->pos,tagReadNode->dir);
				}
				tagReadNode = tagReadNode->next;
			}
			printf("\n");
		}
	}

	printf("\n");
}

void copyTag(struct partGraph *parentReadIds, struct partGraph *parentReadIdsTmp){
	int i,j;

	for(i = 0; i < parentReadIdsTmp->V;i++){
		for(j = 0; j < parentReadIds->V;j++){
			if(parentReadIdsTmp->array[i].ID == parentReadIdsTmp->array[j].ID){

			}
		}
	}
}

void printUpstreamTree(struct partGraph *tree){
	int i;
	struct edge *edge;

	printf("Tree for node: %i\n",tree->array[0].ID);
	for(i=0;i<tree->V;i++){
		printf("Node: %i at %i (depth: %i) flag: %i\n",tree->array[i].ID,i,tree->array[i].depth,tree->array[i].flag);
		printf("Parents:");
		edge = tree->array[i].tail;
		while(edge){
			printf("  -> %i (%i)",tree->array[edge->dest].ID,edge->dest);
			edge = edge->next;
		}
		printf("\nChildren:");
		edge = tree->array[i].head;
		while(edge){
			printf("  -> %i (%i)",tree->array[edge->dest].ID,edge->dest);
			edge = edge->next;
		}
		printf("\n");
	}
	printf("\n");
	printf("\n");

//	for(i=0;i<tree->V;i++){
//		edge = tree->array[i].tail;
//		while(edge){
//			printf("-> %i",tree->array[edge->dest].ID);
//			edge = edge->next;
//		}
//		printf("   ");
//
//	}
//	printf("\n");


}

struct KannteNode* tagKannte;

/**
 * Function proofs the correctness of the read direction in the dBG. The found end have to be downstream in the right distance of the second end.
 *
 * @param i
 * @param stp			Start position of the read in the node (NULL initially)
 * @param endp			End position of the read in the node (NULL initially)
 * @param stk			Node of the start point of the read (NULL initially)
 * @param endk			Node of the end point of the read (NULL initially)
 * @param rnode			The read start or end point of a read at a certain position
 * @param parentReadIds	Partial tree upwards the actual position
 * @return				1 if the reads has the right direction at this position (next end is upstream), 0 otherwise
 */
int endStart(int tag,int *stp, int *endp, int *stk, int *endk, struct ReadNode *rnode, struct partGraph *childReadIds,int *reallen){
	static struct KannteNode *kannte;
	static struct edge *edge;
	static struct ReadNode* tagReadNode;
	kannte = rnode->read->headkannte;
	static int mypos;
	mypos = 0;
	static int j;
	static int in;
	static int Id;
	in = 0;
	// Find correct start pos
//	printf("rnode: %i on %i\n",rnode->read->ID,i);

	static int mylen = 0;
	mylen = readLenList[rnode->read->ID];
	int relen = 0;
//	int nodeID;


	for(j=0;j<childReadIds->V;j++){
		kannte = rnode->read->headkannte;
		while(kannte){
			if(childReadIds->array[j].ID == kannte->dest){
//				if(rnode->read->ID == 37404){
////					printUpstreamTree(childReadIds);
//					printf("Test (V: %i)\n",childReadIds->V);
//					printf("\tCompare: %i == %i and %i == %i (kannteDest: %i)\n",rnode->pos,kannte->pos, childReadIds->array[j].ID,(*stk),kannte->dest);
//				}
//				printf("\tCompare: %i == %i and %i == %i\n",rnode->pos,kannte->pos, childReadIds->array[j].ID,(*stk));
				if(!(rnode->pos == kannte->pos && childReadIds->array[j].ID==(*stk))){
					// kante ist andere seite
					(*endk) = childReadIds->array[j].ID;
					(*endp) = kannte->pos;
					relen = (childReadIds->array[j].depth - redGraph->array[childReadIds->array[j].ID].len) + kannte->pos + (*stp);
//					(*endp) = kannte->pos - childReadIds->array[j].depth;
					// +/- 1 in case of InDels
					if(relen == mylen){
						if(tag){

							// Tag all ReadNodes on the path explicitly
							// Depth of the tagged area could be used as value to set the minOLVlen

							if(rnode->read->ID == TESTREAD){
								printf("Before Tags\n");
//								printUpstreamTree(childReadIds);
								printTag(childReadIds,rnode);
							}
							tagKannte = kannte;
							tagReadNode = kannte->ReadNode;
							childReadIds->array[j].flag = 1;

							if(rnode->read->ID == TESTREAD) printf("j = 0 ? : %i (Startpos: %i)\n",j,tagReadNode->pos);

//							if(rnode->read->ID == TESTREAD){
//								printTag(childReadIds,rnode);
//							}

							while(tagReadNode){
								tagReadNode->flag = 1;
								if(j==0 && tagReadNode->pos <= rnode->pos){
//									printf("j == 0\n");
									break;
								}
								tagReadNode = tagReadNode->next;
							}

							if(j!=0){
								if(rnode->read->ID == TESTREAD) printf("Tag from: %i -> %i\n", childReadIds->array[j].ID,childReadIds->array[0].ID);
								for(in = j; in > 0; in--){
									if(childReadIds->array[in].flag){
//										edge = childReadIds->array[in].tail;

										edge = childReadIds->array[in].head;
										while(edge){
											if(edge->dest != 0 && edge->dest < in){
												if(rnode->read->ID == TESTREAD) printf("Tagged child node: %i\n",childReadIds->array[edge->dest].ID);
												if(rnode->read->ID == TESTREAD) printf("Tag in Node: %i (in (%i) -> edge->dest = %i)\n",childReadIds->array[edge->dest].ID,in,edge->dest);
												childReadIds->array[edge->dest].flag = 1;
												Id = childReadIds->array[edge->dest].ID;
												tagReadNode = redGraph->array[Id].headread;
												while(tagReadNode){
													tagReadNode->flag = 1;
													tagReadNode = tagReadNode->next;
												}
											}
											edge = edge->next;
										}

									}
								}

								Id = childReadIds->array[0].ID;
								tagReadNode = redGraph->array[Id].headread;
								childReadIds->array[0].flag = 1;
								while(tagReadNode && tagReadNode->pos >= rnode->pos){
									tagReadNode->flag = 1;
									if(rnode->read->ID == TESTREAD) printf("Tagged pos %i: %i (of %i) -> Read: %i\n",Id,tagReadNode->pos,rnode->pos,tagReadNode->read->ID);
									tagReadNode = tagReadNode->next;
								}
							}


						}
						else{
							if(kannte->ReadNode->flag){
								(*reallen) = relen;
								if((*stk) == TESTNODE || (*stk) == TESTNODE2) printf("(j: %i)ID(%i): %i -> stp: %i (k: %i) endp: %i (k: %i) (mylen: %i / relen: %i)\n",j,tag,rnode->read->ID,(*stp),(*stk),(*endp),(*endk), mylen, relen);
								return 1;
							}
							else{
								(*reallen) = relen;
								if((*stk) == TESTNODE || (*stk) == TESTNODE2) printf("(j: %i)ID(%i) NOT FLAGGED: %i -> stp: %i (k: %i) endp: %i (k: %i) (mylen: %i / relen: %i)\n",j,tag,rnode->read->ID,(*stp),(*stk),(*endp),(*endk), mylen, relen);
								return 2;
							}
							return 0;
						}
						if((*stk) == TESTNODE || (*stk) == TESTNODE2) printf("(j: %i)ID(%i): %i -> stp: %i (k: %i) endp: %i (k: %i) (mylen: %i / relen: %i)\n",j,tag,rnode->read->ID,(*stp),(*stk),(*endp),(*endk), mylen, relen);
						(*reallen) = relen;
//						printf("Return true (rNode)\n");
						return 1;
					}
				}
			}
			kannte = kannte->next;
		}
	}
//	printf("Nothing found!!!\n");
	return 0;
}

/**
 *
 * @param ovlGraph 		Graph represents all transitively irreducible overlaps (initially NULL)
 * @return				0 if some error occurs 1 otherwise
 */
int stringer3(struct myovlList *ovlGraph){
	char verbose = 0;
	char verbose2 = 0;
	char verboseJOvl = 0;

	int i,m,k;
	int rNid = 0, rnNid = 0;
	int kanId;
	int starts=0, ends=0;
	int found = 0;
	int anyfound = 0;
	int ovls = 0;
	int seekexpand = 50;

	struct seekOvl *seek = (struct seekOvl*)malloc(sizeof(struct seekOvl));
	seek->newlen = 0;
	seek->oldlen = 0;
	seek->cbefpos = (int*)malloc(sizeof(int)*seekexpand);
	seek->childs  = (int*)malloc(sizeof(int)*seekexpand);
	seek->parents = (int*)malloc(sizeof(int)*seekexpand);
	seek->pbefpos = (int*)malloc(sizeof(int)*seekexpand);

//	int *oldch = (int*)malloc(sizeof(int)*20);
//	int *newch = (int*)malloc(sizeof(int)*20);

	struct edge *edge;

	struct partGraph *parentReadIds = (struct partGraph*)malloc(sizeof(struct partGraph));
	struct partGraph *parentReadIdsTemp = (struct partGraph*)malloc(sizeof(struct partGraph));
	parentReadIds->V = 0;
	parentReadIdsTemp->V = 0;
	parentReadIds->array = (struct node*)malloc(sizeof(struct node)*1000);
	parentReadIdsTemp->array = (struct node*)malloc(sizeof(struct node)*1000);

	for(i=0;i<1000;i++){
		parentReadIdsTemp->array[i].head = NULL;
		parentReadIds->array[i].head = NULL;
		parentReadIdsTemp->array[i].tail = NULL;
		parentReadIds->array[i].tail = NULL;
		parentReadIdsTemp->array[i].flag = 0;
		parentReadIds->array[i].flag = 0;
	}

	struct ReadNode *rNode,*nextrNode;
	struct ReadNode *tagReadNode;
	int astp, bstp, aendp, bendp; // positions
	int astk, bstk, aendk, bendk; // kannte
	int areallen,breallen;
	int abover,baover;
	int stdif, befpos;
	int tmp;
	int eS;
	int ovlLen;
	struct readEdge* readedge;

	// Redundant
	for(i=1;i<=redGraph->V;i++){
		redGraph->vFlag[i] = 0;
	}

	printf("redGraph->V: %i\n",ovlGraph->V);
	for(i=1;i<=redGraph->V;i++){
		if(redGraph->array[i].headread){
			if(verbose2) printf("Node : %i -> V: %i\n",i,parentReadIds->V);
			for(m=0;m<parentReadIds->V;m++){
				parentReadIds->array[m].flag = 0;
			}
			upstreamTree(i,parentReadIds);
//			if(i==TESTNODE){
//				printf("Before While:\n");
			if(verbose2) printUpstreamTree(parentReadIds);
//			}
#ifndef UPSTREAMTEST
			rNode = redGraph->array[i].headread;
			while(rNode){
				rNid = rNode->read->ID;
//				printf("rNid: %i\n",rNid);
				if(rNid<0) printf("Integer Overflow\n");
//				if(!ovlGraph->read[rNid]){
//					rNode = rNode->next;
//					continue;
//				}
				if(ovlGraph->read[rNid]->flag == CONTAINED){
					rNode = rNode->next;
					continue;
				}
				if(verbose) printf("rNode: %i\n",rNode->read->ID);
//				if(i == TESTNODE) printf("rNode: %i\n",rNode->read->ID);
//				if(i==TESTNODE){
//					printf("In While:\n");
				if(verbose) printUpstreamTree(parentReadIds);
//				}
//				critk = -1;
				rNid = rNode->read->ID;
				astp = redGraph->array[i].len - rNode->pos;
//				reallen = rNode->pos;
				astk = i;
				if(endStart(1,&astp,&aendp,&astk,&aendk,rNode,parentReadIds,&areallen)){
					if(setaRead2(rNid,rNode->dir,ovlGraph)){
						ovlGraph->read[rNid]->stk = astk;
						ovlGraph->read[rNid]->endk = aendk;
						befpos = rNode->pos;
						astp = 0;
						tmp = aendp;
						aendp = areallen;
						stdif = 0;
						found = 0;
						starts++;

						nextrNode = rNode->next;
						while(nextrNode && (!found)){
							rnNid = nextrNode->read->ID;
							bstp = redGraph->array[i].len - nextrNode->pos;
							bstk = i;
							stdif += befpos - nextrNode->pos;
							if(stdif > (maxRlen - MINOVL)){
								found = 1;
								break;
							}
							befpos = nextrNode->pos;
							if((eS = endStart(0,&bstp,&bendp,&bstk,&bendk,nextrNode,parentReadIds,&breallen))){
								if(eS == 2){
									if(!(tmp == bendp && aendk == bendk)){
										nextrNode = nextrNode->next;
										continue;
									}
									else{
//										printf("Rescue Overlap by eS == 2\n");
									}
								}
								bstp = astp - stdif;
								bendp = bstp + breallen;

								ovls++;
								abover = aendp - bendp;
								baover = astp - bstp;
								// MinOvlLen
								ovlLen = nK + (_min(areallen,breallen) - _min(abover,baover));
//								if(astk == TESTNODE || astk == TESTNODE2) printf("Overlap-Length: %i\n",ovlLen);
//								printf("Overlap-Length: %i\n",ovlLen);
//								printf("MinOverlap-Length: %i\n",ovlLen > MINOVL);
								if(ovlLen > MINOVL){
									if(verbose) printf("Overlap-Length: %i\n",ovlLen);
									k = setbRead2(rNid,rnNid,nextrNode->dir,ovlGraph,abover,baover,1);
									ovlGraph->read[rnNid]->stk = bstk;
									ovlGraph->read[rnNid]->endk = bendk;
									if(k == 2){
//										if(astk == TESTNODE || astk == TESTNODE2)
										if(verbose) printf("Overlap (PROPER): %i (dir: %i) -> % i (dir: %i) --->>> k: %i (abover: %i baover: %i)\n",rNid,rNode->dir,rnNid,nextrNode->dir,k,abover,baover);
										found = 1;
										break;
									}
									if(k == 1){
//										if(astk == TESTNODE || astk == TESTNODE2)
										if(verbose) printf("Overlap (b is CONTAINED): %i -> %i --->>> k: %i (abover: %i baover: %i)\n",rNid,rnNid,k,abover,baover);
										nextrNode = nextrNode->next;
										continue;
									}
									if(k == 0){
										found = 1;
//										if(astk == TESTNODE || astk == TESTNODE2)
										if(verbose) printf("Overlap (a is CONTAINED): %i -> %i --->>> k: %i (abover: %i baover: %i)\n",rNid,rnNid,k,abover,baover);
										break;
									}
								}
							}
							// Q: When stop searching? (Max readlen?)
							nextrNode = nextrNode->next;
						}
						if(!found){
							stdif = rNode->pos;
							seek->childs[seek->newlen] = i;
							seek->cbefpos[seek->newlen] = stdif;
							seek->newlen++;

							while(seek->newlen){
								memcpy(seek->parents,seek->childs,sizeof(int)*seek->newlen);
								memcpy(seek->pbefpos,seek->cbefpos,sizeof(int)*seek->newlen);
								seek->oldlen = seek->newlen;
								seek->newlen = 0;

								while(seek->oldlen){
									seek->oldlen--;
									edge = redGraph->array[seek->parents[seek->oldlen]].tail;
									anyfound = 0;
									while(edge){
										readedge = redGraph->array[aendk].readedge;
										while(readedge){
											if(readedge->KannteID == edge->dest){
												anyfound++;
												break;
											}
											readedge = readedge->next;
										}
										edge = edge->next;
									}
									edge = redGraph->array[seek->parents[seek->oldlen]].tail;
									while(edge){
										//___
										readedge = redGraph->array[aendk].readedge;
										while(readedge){
											if(readedge->KannteID == edge->dest) break;
											readedge = readedge->next;
										}
//										if(!readedge){
//											found = 1;
//											edge = edge->next;
//											continue;
//										}
										//___
										found = 0;
										kanId = edge->dest;
										befpos = redGraph->array[kanId].len;
										stdif = seek->pbefpos[seek->oldlen];
										for(m=0;m<parentReadIdsTemp->V;m++){
											parentReadIdsTemp->array[m].flag = 1;
										}
										upstreamTree(kanId,parentReadIdsTemp);
										nextrNode = redGraph->array[kanId].headread;
										while(nextrNode){
											stdif += befpos - nextrNode->pos;
											if(stdif > (maxRlen - MINOVL)){
												found = 1;
												break;
											}
											rnNid = nextrNode->read->ID;
											bstp = redGraph->array[kanId].len - nextrNode->pos;
											bstk = kanId;


											if((eS = endStart(0,&bstp,&bendp,&bstk,&bendk,nextrNode,parentReadIdsTemp,&breallen))){
												if(eS == 2){
													if(!(tmp == bendp && aendk == bendk)){
														nextrNode = nextrNode->next;
														continue;
													}
													else{
				//										printf("Rescue Overlap by eS == 2\n");
													}
												}
												bstp = astp - stdif;
												bendp = bstp + breallen;
	//											printf("Correct orientation\n");
												ovls++;
												abover = aendp - bendp;
												baover = astp - bstp;
												// MinOvlLen
												ovlLen = nK + (_min(areallen,breallen) - _min(abover,baover));
//												if(astk == TESTNODE || astk == TESTNODE2) printf("Overlap-Length: %i\n",ovlLen);
												if(ovlLen > MINOVL && (readedge || !anyfound)){
													k = setbRead2(rNid,rnNid,nextrNode->dir,ovlGraph,abover,baover,1);
													if(k == 2){
														if(verboseJOvl){
															printf("-> Next Overlap (PROPER): %i -> %i --->>> k: %i (abover: %i (%i - %i) baover: %i (%i - %i))\n",rNid,rnNid,k,abover,aendp,bendp,baover,astp,bstp);
															printf("InfoTest: --->  a: %i -> %i // b: %i -> %i\n",aendk,astk,bendk,bstk);
														}
														found = 1;
														break;
													}
													if(k == 1){
														if(verboseJOvl)
															printf("-> Next Overlap (b is CONTAINED): %i -> %i --->>> k: %i (abover: %i baover: %i)\n",rNid,rnNid,k,abover,baover);
														nextrNode = nextrNode->next;
														continue;
													}
													if(k == 0){
														found = 1;
														if(verboseJOvl)
															printf("-> Next Overlap (a is CONTAINED): %i -> %i --->>> k: %i (abover: %i baover: %i)\n",rNid,rnNid,k,abover,baover);
														break;
													}
												}
											}


											nextrNode = nextrNode->next;
										}
										if(!found && stdif < (maxRlen - MINOVL)){
											if(seek->newlen == seekexpand-1){
												seekexpand *= 2;
												seek->childs = (int*)realloc(seek->childs,sizeof(int)*seekexpand);
												seek->cbefpos = (int*)realloc(seek->cbefpos,sizeof(int)*seekexpand);
												seek->parents = (int*)realloc(seek->parents,sizeof(int)*seekexpand);
												seek->pbefpos = (int*)realloc(seek->pbefpos,sizeof(int)*seekexpand);
											}
											seek->childs[seek->newlen] = edge->dest;
											seek->cbefpos[seek->newlen] = seek->pbefpos[seek->oldlen] + redGraph->array[edge->dest].len;
											seek->newlen++;
//											printf("Not Found, go deeper, depth: %i (newlen: %i) (from: %i to %i)\n",stdif,seek->newlen,seek->parents[seek->oldlen],edge->dest);
										}
										edge = edge->next;
									}
								}
								if(seek->newlen > 500){
									seek->newlen = 0;
								}
							}
						}
					}

//					printf("Delete Tags\n");
					tmp = 0;
					for(m=0; m<parentReadIds->V; m++){
						if(parentReadIds->array[m].flag) tmp = m;
					}

					// EndNode
					if(tmp){
						tagReadNode = tagKannte->ReadNode;
						while(tagReadNode){
							tagReadNode->flag = 0;
							tagReadNode = tagReadNode->next;
						}
						parentReadIds->array[tmp].flag = 0;

						for(m = tmp-1 ; m > 0 ; m--){
							if(parentReadIds->array[m].flag){
								tagReadNode = redGraph->array[parentReadIds->array[m].ID].headread;
								while(tagReadNode){
									tagReadNode->flag = 0;
									tagReadNode = tagReadNode->next;
								}
								parentReadIds->array[m].flag = 0;
							}
						}

						tagReadNode = redGraph->array[parentReadIds->array[0].ID].headread;
						while(tagReadNode && tagReadNode->flag){
							tagReadNode->flag = 0;
							tagReadNode = tagReadNode->next;
						}
						parentReadIds->array[0].flag = 0;
					}
					else{
//							tagReadNode = redGraph->array[0].headread;
						tagReadNode = tagKannte->ReadNode;
						if(tagReadNode->read->ID == TESTREAD) printf("Startpos for deletion: %i\n",tagReadNode->pos);
						while(tagReadNode && tagReadNode->pos >= rNode->pos){
//								if(!tagReadNode->flag){
//									printf("!!! FALSCH !!!");
//								}
							tagReadNode->flag = 0;
							tagReadNode = tagReadNode->next;
						}
						parentReadIds->array[0].flag = 0;
					}

//					printTag(parentReadIds,rNode);
				}
				else{
					ends++;
				}

				rNode = rNode->next;
			}
#endif
		}
	}

	free_UpstreamTree(parentReadIds);
	free_UpstreamTree(parentReadIdsTemp);
	free(seek->cbefpos);
	free(seek->childs);
	free(seek->parents);
	free(seek->pbefpos);
	free(seek);

	printf("STARTS: %i\n",starts);
	printf("ENDS  : %i\n",ends);

	return 1;

}

void downstreamTreeCon(int i, struct partGraph* childReadIds){
	//	redGraph->vFlag[i] = 0;
		childReadIds->array[0].ID = i;
		childReadIds->array[0].flag = 0;
		childReadIds->array[0].depth = 0;
		childReadIds->array[0].layer = 0;
		int j=0,k=1,l=0;
		int end = 0;
		struct edge *child;
		struct edge *newEdge;
		int jNode;

		while(childReadIds->array[j].tail){
			newEdge = childReadIds->array[j].tail;
			childReadIds->array[j].tail = childReadIds->array[j].tail->next;
			free(newEdge);
		}

		do{
			jNode = childReadIds->array[j].ID;
			child = redGraph->array[jNode].tail;
	//		printf("jNode: %i\n",jNode);
			while(child && childReadIds->array[j].depth < MAX_DEPTH){
	//			printf("Children: -> %i\n",child->dest);
				if(redGraph->vFlag[child->dest] != 1){
	//				printf("If Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
					redGraph->vFlag[child->dest] = 1;
					childReadIds->array[k].ID = child->dest;
					childReadIds->array[k].layer = childReadIds->array[j].layer+1;
					childReadIds->array[k].depth = redGraph->array[child->dest].len + childReadIds->array[j].depth;

					if(childReadIds->array[k].head){
						while(childReadIds->array[k].head->next){
							newEdge = childReadIds->array[k].head->next->next;
							free(childReadIds->array[k].head->next);
							childReadIds->array[k].head->next = newEdge;
						}
						childReadIds->array[k].head->dest = j;
					}
					else{
						newEdge = (struct edge*)malloc(sizeof(struct edge));
						newEdge->next = childReadIds->array[k].head;
						newEdge->dest = j;
						childReadIds->array[k].head = newEdge;
					}
					while(childReadIds->array[k].tail){
						newEdge = childReadIds->array[k].tail;
						childReadIds->array[k].tail = childReadIds->array[k].tail->next;
						free(newEdge);
					}
					newEdge = (struct edge*)malloc(sizeof(struct edge));
					newEdge->next = childReadIds->array[j].tail;
					newEdge->dest = k;
					childReadIds->array[j].tail = newEdge;
	//				printf("Edge depth: %i -> %i\n",childReadIds->array[k].ID,childReadIds->array[k].depth);
					k++;
				}
//				else{
//	//				printf("Else Casus(i,j,k,l: %i %i %i %i)\n",i,j,k,l);
//					for(l=0;l<k;l++){
//						if(childReadIds->array[l].ID == child->dest){
//							newEdge = (struct edge*)malloc(sizeof(struct edge));
//							newEdge->dest = l;
//							newEdge->next = childReadIds->array[j].tail;
//							childReadIds->array[j].tail = newEdge;
//
//							newEdge = (struct edge*)malloc(sizeof(struct edge));
//							newEdge->dest = j;
//							newEdge->next = childReadIds->array[l].head;
//							childReadIds->array[l].head = newEdge;
//						}
//					}
//				}
				child = child->next;
			}
			j++;
			if(k>40 || j==k) end = 1;
		} while(!end);

		childReadIds->V = k;
	//	printf("TempGraph for node: %i (V: %i)\n",i,k);
	//	int layer = 0;
		for(l=0;l<k;l++){
	//		if(childReadIds->array[l].layer != layer){
	//			layer++;
	////			 printf("   ");
	//		}
	//		 printf("%i ",childReadIds->array[l].ID);
			redGraph->vFlag[childReadIds->array[l].ID] = 0;
		}
	//	 printf("\n");

	//	struct edge *edge;
	//	for(j=0;j<childReadIds->V;j++){
	//		printf("Test for right edges for node: %i",childReadIds->array[j].ID);
	//		printf("\nHead: ");
	//		edge = childReadIds->array[j].head;
	//		while(edge){
	//			printf("->%i ",childReadIds->array[edge->dest].ID);
	//			edge = edge->next;
	//		}
	//		printf("\nTail: ");
	//		edge = childReadIds->array[j].tail;
	//		while(edge){
	//			printf("->%i ",childReadIds->array[edge->dest].ID);
	//			edge = edge->next;
	//		}
	//		printf("\n");
	//	}

}

int startEndCon(int tag,int *stp, int *endp, int *stk, int *endk, struct ReadNode *rnode, struct partGraph *childReadIds){
	static struct KannteNode *kannte;
	kannte = rnode->read->headkannte;
	static int relen;
	static int j;
	static int truedir=0;
	static int falsedir=0;
//	static nothing=0;

	static int mylen = 0;
	mylen = readLenList[rnode->read->ID];

	for(j=0;j<childReadIds->V;j++){
		kannte = rnode->read->headkannte;
		while(kannte){
			if(childReadIds->array[j].ID == kannte->dest){
				if(!(rnode->pos == kannte->pos && childReadIds->array[j].ID==(*stk))){
					(*endk) = childReadIds->array[j].ID;
					if(!tag){
						(*endp) = kannte->pos;
						relen = rnode->pos + (childReadIds->array[j].depth - kannte->pos);
					}
					else {
						(*endp) += childReadIds->array[j].depth - kannte->pos;
						relen = (*endp) - (*stp);
						if(rnode->read->ID == 10555) printf("mylen: %i relen: %i\n",mylen,relen);
					}
//					relen = (childReadIds->array[j].depth - redGraph->array[childReadIds->array[j].ID].len) + kannte->pos + (*stp);
					if(relen == mylen){
						childReadIds->array[j].flag = 1;
						truedir++;
						(*endk) = j;
//						printf("ID: %i -> relen: %i = mylen: %i (depth: %i) rnodepos: %i kanntepos: %i (found: %i, notfound: %i)\n",rnode->read->ID,relen,mylen,childReadIds->array[j].depth,rnode->pos,kannte->pos,truedir,falsedir);
						tagKannte = kannte;
						return 2;
					}
					else falsedir++;
//					printf("ID: %i -> relen: %i = mylen: %i (depth: %i) rnodepos: %i kanntepos: %i (found: %i, notfound: %i)\n",rnode->read->ID,relen,mylen,childReadIds->array[j].depth,rnode->pos,kannte->pos,truedir,falsedir);
//					return 1;
				}
			}
			kannte = kannte->next;
		}
	}
//	nothing++;
//	printf("Nothing for ID: %i num: %i\n",rnode->read->ID, nothing);
	return 0;
}

int startEndCon2(int *stp, int *endp, int *stk, int *endk, struct ReadNode *rnode, struct partGraph *childReadIds,int *nodepos){
	static struct KannteNode *kannte;
	kannte = rnode->read->headkannte;
	static int relen;
	static int j;
	static int truedir=0;
	static int falsedir=0;
//	static nothing=0;

	static int mylen = 0;
	mylen = readLenList[rnode->read->ID];


	for(j=0;j<childReadIds->V;j++){
		kannte = rnode->read->headkannte;
		while(kannte){
			if(childReadIds->array[j].ID == kannte->dest){
				if(!(rnode->pos == kannte->pos && childReadIds->array[j].ID==(*stk))){
					(*endk) = childReadIds->array[j].ID;
					(*endp) = childReadIds->array[j].depth - kannte->pos;
					(*nodepos) = kannte->pos;
					relen = (*endp) - (*stp);
//					printf("ID: %i -> relen: %i = mylen: %i (depth: %i)\n",rnode->read->ID,relen,mylen,childReadIds->array[j].depth);
//					relen = (childReadIds->array[j].depth - redGraph->array[childReadIds->array[j].ID].len) + kannte->pos + (*stp);
					if(relen <= mylen + 5 && relen >= mylen - 5){
						childReadIds->array[j].flag = 1;
						truedir++;
						(*endk) = j;
//						printf("ID: %i -> relen: %i = mylen: %i (depth: %i) rnodepos: %i kanntepos: %i (found: %i, notfound: %i)\n",rnode->read->ID,relen,mylen,childReadIds->array[j].depth,rnode->pos,kannte->pos,truedir,falsedir);
						tagKannte = kannte;
						return 2;
					}
					else falsedir++;
					if(rnode->read->ID == 1144414) printf("ID: %i -> relen: %i = mylen: %i (depth: %i) rnodepos: %i kanntepos: %i (found: %i, notfound: %i)\n",rnode->read->ID,relen,mylen,childReadIds->array[j].depth,rnode->pos,kannte->pos,truedir,falsedir);
//					return 1;
				}
			}
			kannte = kannte->next;
		}
	}
//	nothing++;
//	printf("Nothing for ID: %i num: %i\n",rnode->read->ID, nothing);
	return 0;
}

void tagPath(struct ReadNode *rnode, int endp, int endk, struct partGraph *childReadIds){
	int i;
	int ID;
	struct edge *edge;

	if(endk == 0){
		while(rnode && rnode->pos >= endp){
			rnode->flag = 1;
			rnode = rnode->next;
		}
	}
	else{
		while(rnode){
			rnode->flag = 1;
			rnode = rnode->next;
		}
		childReadIds->array[0].flag = 1;
		ID = childReadIds->array[endk].ID;
		rnode = redGraph->array[ID].headread;
		while(rnode && rnode->pos >= endp){
			rnode->flag = 1;
			rnode = rnode->next;
		}


		for(i = endk; i > 0; i--){
			if(childReadIds->array[i].flag){
				edge = childReadIds->array[i].head;
				while(edge && childReadIds->array[edge->dest].flag == 0){
					ID = childReadIds->array[edge->dest].ID;
					childReadIds->array[edge->dest].flag = 1;
					rnode = redGraph->array[ID].headread;
					while(rnode){
						rnode->flag = 1;
						rnode = rnode->next;
					}
					edge = edge->next;
				}
			}
		}
	}


}

void detagPath(struct ReadNode *rnode, int endp, int endk, struct partGraph *childReadIds){
	int i;
	int ID;
	struct edge *edge;

	if(endk == 0){
		while(rnode && rnode->pos >= endp){
			rnode->flag = 0;
			rnode = rnode->next;
		}
	}
	else{
		while(rnode){
			rnode->flag = 0;
			rnode = rnode->next;
		}
		childReadIds->array[0].flag = 0;

		ID = childReadIds->array[endk].ID;
		rnode = redGraph->array[ID].headread;
		while(rnode && rnode->pos >= endp){
			rnode->flag = 0;
			rnode = rnode->next;
		}

		for(i = endk; i > 0; i--){
			if(childReadIds->array[i].flag){
				childReadIds->array[i].flag = 0;
				edge = childReadIds->array[i].head;
				while(edge && childReadIds->array[edge->dest].flag == 1){
					ID = childReadIds->array[edge->dest].ID;
					rnode = redGraph->array[ID].headread;
					while(rnode){
						rnode->flag = 0;
						rnode = rnode->next;
					}
					edge = edge->next;
				}
			}
		}
	}
}


#define CONTAIN_MASSAGE

/**
 * Not set reads (reads not included in the myovlList, because no proper start-point was found) are set flagged as contained, without a cointained overlap
 */
inline static void catchNonSetReads(struct myovlList *ovlGraph){
	int i;
	for(i=1;i<=ovlGraph->V;i++){
		if(!ovlGraph->read[i]){
			ovlGraph->read[i] = (struct aread*)malloc(sizeof(struct aread));
			ovlGraph->read[i]->dir = 0;
			ovlGraph->read[i]->first = NULL;
			ovlGraph->read[i]->length = readLenList[i];
			ovlGraph->read[i]->flag = CONTAINED;
		}
	}
}

void tag_A_Contained(struct myovlList *ovlGraph){
	char verbose = 0;

	int i,j;
	int rNid,rnNid;
	int stk, endk, stp, endp;
	int bstk, bendk, bendp, bstp;
	int abover, baover;
	int nodepos;
	int detagnodepos;
	int flagged;
	struct partGraph *childReadIds = (struct partGraph*)malloc(sizeof(struct partGraph));
	childReadIds->V = 0;
	childReadIds->array = (struct node*)malloc(sizeof(struct node)*1000);
	for(i=0;i<1000;i++){
		childReadIds->array[i].head = NULL;
		childReadIds->array[i].tail = NULL;
		childReadIds->array[i].flag = 0;
	}

	struct ReadNode *rNode, *nextrNode;
	struct KannteNode* kannte;

	if(verbose) printf("V: %i\n",redGraph->V);

	for(i = 1; i <= redGraph->V; i++){
		redGraph->vFlag[i] = 0;
	}

	for(i = 1; i <= redGraph->V; i++){
		if(redGraph->array[i].headread){
//			printf("TESTTTTT\n");
			downstreamTreeCon(i,childReadIds);
			for(j=0;j<childReadIds->V;j++){
				childReadIds->array[j].flag = 0;
			}
//			printUpstreamTree(childReadIds);
			stk = i;
			rNode = redGraph->array[i].headread;
			while(rNode){
				rNid = rNode->read->ID;
//				stp = rNode->pos;
				stp = 0 - rNode->pos;
//				// Test for length
				if(startEndCon2(&stp,&endp,&stk,&endk,rNode,childReadIds,&nodepos) == 2){
					if(setaRead2(rNid,rNode->dir,ovlGraph)){
//						tagPath(rNode,endp,endk,childReadIds);
						tagPath(rNode,nodepos,endk,childReadIds);
						detagnodepos = nodepos;
//						printf("Tag: ID: %i nodepos: %i , endk: %i\n",rNode->read->ID,nodepos,endk);
//						if(rNid == 6216){
//							printf("Print Tag\n");
//							printTag(childReadIds,rNode);
//						}
//						printf("TAG PATH\n");
						nextrNode = rNode->next;
						while(nextrNode && nextrNode->flag){
							if(setaRead2(nextrNode->read->ID,nextrNode->dir,ovlGraph)){
								kannte = nextrNode->read->headkannte;
								flagged = 1;
								while(kannte){
									if(!kannte->ReadNode->flag) flagged = 0;
									kannte = kannte->next;
								}
								if(flagged && ovlGraph->read[nextrNode->read->ID]->flag != CONTAINED){
//									bstp = nextrNode->pos;
									bstp = 0 - nextrNode->pos;
									bstk = i;
									if(startEndCon2(&bstp,&bendp,&bstk,&bendk,nextrNode,childReadIds, &nodepos) == 2){
										abover = bendp - endp;
										baover = bstp - stp;
										setbRead2(rNid,nextrNode->read->ID,nextrNode->dir,ovlGraph,abover,baover,1);
//										printf("CONTAINMENT: %i is in %i (ab: %i (%i - %i) -> ba: %i)\n",nextrNode->read->ID,rNode->read->ID,abover,bendp,endp,baover);
//										ovlGraph->read[nextrNode->read->ID]->dir = !nextrNode->dir;
//										ovlGraph->read[nextrNode->read->ID]->flag = CONTAINED;
//										ovlGraph->read[nextrNode->read->ID]->endk = bstk;
//										ovlGraph->read[nextrNode->read->ID]->stk = bendk;
										// set as b read
									}
								}
							}
							nextrNode = nextrNode->next;
						}
						// Search in flagged child Node for further contained reads if the path is tagged by the path of the a-read
						if(!nextrNode && redGraph->array[i].tail){
#ifdef CONTAIN_MASSAGE
							if(verbose) printf("Search in ChildNodes for A-Read (%i) Containments in Node: %i\n",rNid,i);
#endif
//							printf("Print Tag\n");
//							printTag(childReadIds,rNode);
							for(j=1;j<childReadIds->V;j++){
								if(childReadIds->array[j].flag){
									rnNid = childReadIds->array[j].ID;
//									printf("Child %i is tag for Containment search\n",rnNid);
									nextrNode = redGraph->array[rnNid].headread;
									while(nextrNode && nextrNode->flag){
										if(setaRead2(nextrNode->read->ID,0,ovlGraph)){
											kannte = nextrNode->read->headkannte;
											flagged = 1;
//											printf("ID: %i\n",nextrNode->read->ID);
											while(kannte){
//												printf("-> %i (%i) ",kannte->pos,kannte->ReadNode->flag);
												if(!kannte->ReadNode->flag) flagged = 0;
												kannte = kannte->next;
											}
//											printf("\n");
											if(flagged && ovlGraph->read[nextrNode->read->ID]->flag != CONTAINED){
												bstp = childReadIds->array[j].depth - nextrNode->pos;
												bstk = rnNid;
												if(startEndCon2(&bstp,&bendp,&bstk,&bendk,nextrNode,childReadIds, &nodepos) == 2){
													abover = bendp - endp;
													baover = bstp - stp;
													setbRead2(rNid,nextrNode->read->ID,nextrNode->dir,ovlGraph,abover,baover,1);
#ifdef CONTAIN_MASSAGE
													if(verbose) printf("\t -> Further contained read: %i (node: %i) in %i (node: %i)\n", nextrNode->read->ID,rnNid,rNode->read->ID,i);
#endif
//													ovlGraph->read[nextrNode->read->ID]->dir = !nextrNode->dir;
//													ovlGraph->read[nextrNode->read->ID]->flag = CONTAINED;
//													ovlGraph->read[nextrNode->read->ID]->endk = bstk;
//													ovlGraph->read[nextrNode->read->ID]->stk = bendk;
												}
//												if(startEndCon2(1,&bstp,&bendp,&bstk,&bendk,nextrNode,childReadIds, &nodepos) == 2){
////													printf("\t -> Further contained read: %i (node: %i) in %i (node: %i)\n", nextrNode->read->ID,rnNid,rNode->read->ID,i);
//													ovlGraph->read[nextrNode->read->ID]->dir = !nextrNode->dir;
//													ovlGraph->read[nextrNode->read->ID]->flag = CONTAINED;
//													ovlGraph->read[nextrNode->read->ID]->endk = bstk;
//													ovlGraph->read[nextrNode->read->ID]->stk = bendk;
//			//										printf("CONTAINMENT: %i is in %i\n",nextrNode->read->ID,rNode->read->ID);
//													// set as b read
//												}
//												printf("b-pos: %i -> %i\n",bstp,bendp);
											}
										}
										nextrNode = nextrNode->next;
									}
								}
							}
						}
						// What ever this is: but second version does not work
//						detagPath(rNode,endp,endk,childReadIds);
						detagPath(rNode,detagnodepos,endk,childReadIds);
//						printf("DeTag: ID: %i nodepos: %i , endk: %i\n",rNode->read->ID,detagnodepos,endk);
//						printTag(childReadIds,rNode);

					}
//					// Test for Containment flag
				}
				rNode = rNode->next;
			}
		}
	}

	free_UpstreamTree(childReadIds);
	catchNonSetReads(ovlGraph);
}

void freeMyOvlList(struct myovlList* G, struct string_graph* S){
	struct bread* bread;
	struct bread* breadN;
	for(int i=0;i<=G->V;i++){
		if(G->read[i]){
			bread = G->read[i]->first;
			if(bread){
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

}

