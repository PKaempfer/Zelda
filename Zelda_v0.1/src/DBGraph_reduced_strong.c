/*
 ============================================================================
 Name        : DBGraph_reduced_strong.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Collapse unique paths in dBG
 ============================================================================
 */


#include "DBGraph_reduced.h"
#include <stdlib.h>
#include <stdio.h>
#include "DBGraph.h"

// Max length of spur to collapse
#define MAXSPURLEN nK

void collapseEdges_strong(int i, int dest, int up,int* nodeList, int nodenum){
	// First step -> Copy read information of the spur edge to the sibling edge
	// a is true --> b is spur
	// up -> nodes have a common parent; !up -> nodes have a common child
	struct ReadNode *anode, *abefnode,*bnode,*temp,*temp2;
	int alen,blen;
//	int abdiff;
	struct KannteNode *kannte;
//	struct KannteNode *befkannte;

//	anode = redGraph->array[dest].headread;
	abefnode = NULL;
	bnode = redGraph->array[i].headread;
	alen = redGraph->array[dest].len;
	blen = redGraph->array[i].len;
//	int templenb = blen;
	int j;

	if(up) j = 0;
	else j = (nodenum-1) * 2;
	bnode = redGraph->array[i].headread;
	anode = redGraph->array[nodeList[j]].headread;

	while(bnode){
//		printf("bnode: -> %i (ID: %i) flag: %i \n",bnode->pos,bnode->read->ID,bnode->flag);
		if(bnode->pos >= nodeList[j+1]){
			if(!anode){
				while(kannte){
					if(kannte->dest == i && kannte->pos == bnode->pos){
						kannte->dest = nodeList[j];
						kannte->pos = bnode->pos - nodeList[j+1];
					}
					kannte = kannte->next;
				}

				temp = bnode->next;
				redGraph->array[nodeList[j]].headread = bnode;
				bnode->pos = bnode->pos - nodeList[j+1];
				bnode->next = NULL;
				anode = bnode;

				bnode = temp;
			}
			else{
				while(anode){
//					printf("\t Anode: %i -> Pos: %i in Node %i (ges: %i)\n",anode->read->ID,anode->pos,nodeList[j],nodeList[j+1]);
					if(!anode->next || anode->next->pos <= bnode->pos - nodeList[j+1]){
//						printf("Transfer bnode to a\n");
						kannte = bnode->read->headkannte;
						while(kannte){
							if(kannte->dest == i && kannte->pos == bnode->pos){
								kannte->dest = nodeList[j];
								kannte->pos = bnode->pos - nodeList[j+1];
							}
							kannte = kannte->next;
						}
						temp = anode->next;
						anode->next = bnode;

						temp2 = bnode->next;
						bnode->next = temp;
						bnode->pos = bnode->pos - nodeList[j+1];

						anode = bnode;
						bnode = temp2;
						if(!bnode) break;
//						printf("Bnode: -> %i (ID: %i) flag: %i \n",bnode->pos,bnode->read->ID,bnode->flag);
					}
					else{
						anode = anode->next;
					}
				}
			}

		}
		else{
			printf("Else, no bnode update\n");
			if(up) j += 2;
			else j -= 2;
			anode = redGraph->array[nodeList[j]].headread;
		}
	}
	redGraph->array[i].headread = NULL;
}

static inline int* setDestList(int i, int dest, int ilen, int destlen, struct edge* shared, int up, int* list, int* listlen){
	int templena = ilen-destlen;
	int tempdest;
	int num = 0;
	char verbose = 0;
	struct edge* tempedge;
	if(up){
		tempedge = redGraph->array[dest].head;
		if(verbose) printf("Collapse upSpur: %i (Parent: %i) over more than on node: %i",i,shared->dest,dest);
	}
	else{
		tempedge = redGraph->array[dest].tail;
		if(verbose) printf("Collapse downSpur: %i (Parent: %i) over more than on node: %i",i,shared->dest,dest);
	}

	list[num++] = dest;
	if(up) list[num++] = 0;
	else list[num++] = ilen - destlen;


	while(tempedge && templena>0 && !tempedge->next){
		tempdest = tempedge->dest;
		list[num++] = tempdest;
		if(up)	list[num] = list[num-2] + destlen;
		else list[num] = list[num-2] - redGraph->array[tempdest].len;
		destlen = redGraph->array[tempdest].len;
		num++;
		templena -= destlen;
		if(verbose) printf(" --> %i",tempdest);
		if(up)	tempedge = redGraph->array[tempdest].head;
		else tempedge = redGraph->array[tempdest].tail;
	}
	if(templena<=0){
		if(verbose) printf("\n");
		(*listlen) = num/2;
	}
	else{
		if(verbose) printf(" --> FAIL\n");
		(*listlen) = -1;
	}
	return list;
}

inline static void deleteSpur(int spur, int up){
	struct edge* edge;
	if(up)	edge = redGraph->array[spur].head;
	else edge = redGraph->array[spur].tail;
	int dest = edge->dest;

	if(up){
		free(redGraph->array[spur].head);
		redGraph->array[spur].head = NULL;
	}
	else{
		free(redGraph->array[spur].tail);
		redGraph->array[spur].tail = NULL;
	}

	struct edge* destedge;
	if(up) destedge = redGraph->array[dest].tail;
	else destedge = redGraph->array[dest].head;
	if(destedge->dest == spur){
		if(up) redGraph->array[dest].tail = destedge->next;
		else redGraph->array[dest].head = destedge->next;
		free(destedge);
	}
	else{
		while(destedge->next){
			if(destedge->next->dest == spur){
				edge = destedge->next;
				destedge->next = destedge->next->next;
				free(edge);
				break;
			}
			destedge = destedge->next;
		}
	}

}

void reduceRedGraph_strong(){
	int i;
	int num, dest;
	int lena=0, lenb=0;
//	int templena;
//	int tempdest;
	struct edge *downedge;
	struct edge *upedge;
	struct edge *tempedge;

	// Linear List but 2 ints per node: ID and Length
	static int* destNodeList = NULL;
	static int  destNodeNum = 0;
	static int  destNodeMaxNum = 100;
	if(!destNodeList) destNodeList = (int*)malloc(sizeof(int*)*destNodeMaxNum*2);

	int loopnum = 1;

	for(i=1;i<=redGraph->V;i++){
		if(loopnum % 1000){
			printf("Loop: %i -> i=%i\n",loopnum,i);
			loopnum++;
		}

		if(redGraph->array[i].len < 2*nK){	// Don't know if nK is a good limit -> think about (seems to be good)
			lena = redGraph->array[i].len;
			upedge = redGraph->array[i].head;
			downedge = redGraph->array[i].tail;

			if((upedge && !upedge->next) && !downedge){
				tempedge = redGraph->array[upedge->dest].tail;
				num = 0;
				dest = -1;
				// Find the true destination
				while(tempedge){
					num++;
					if(tempedge->dest != i){
						dest = tempedge->dest;
						lenb = redGraph->array[dest].len;
					}
					tempedge = tempedge->next;
				}
				if(num == 2 && dest != -1){
					if(lenb >= lena){
						// collapse dest and i in
//						printf("Collapse downSpur: %i --> %i\n",i,dest);
	//					collapseEdges(i,dest,upedge->dest,1);
					}
					else{
						destNodeList = setDestList(i,dest,lena,lenb,upedge,0,destNodeList,&destNodeNum);
						if(destNodeNum != -1) {
//							printf("Nodes: %i\n",destNodeNum);
							collapseEdges_strong(i,dest,1,destNodeList,destNodeNum);
							deleteSpur(i,1);
//							exit(1);
						}
					}
				}
				else if(num>2){
					// Search for the true destination in all kids of ori whose have more further children
					num = 0;
					tempedge = redGraph->array[upedge->dest].tail;
					while(tempedge){
						if(redGraph->array[tempedge->dest].tail){
							dest = tempedge->dest;
							lenb = redGraph->array[dest].len;
							num++;
						}
						tempedge = tempedge->next;
					}
					if(num==1){
						if(lenb >= lena){
							// collapse dest and i in
//							printf("Collapse downSpur (more than on sibling): %i --> %i\n",i,dest);
		//					collapseEdges(i,dest,upedge->dest,1);
						}
						else{
//							printf("Collapse downSpur over more than on node (more than on sibling): %i --> %i\n",i,dest);
//							collapseEdges_strong(i,dest,upedge->dest,1);
						}

					}
				}
			}
			else if((downedge && !downedge->next) && !upedge){
				 // Don't know if nK is a good limit -> think about
				tempedge = redGraph->array[downedge->dest].head;
				num = 0;
				dest = -1;
				lenb = 0;
				while(tempedge){
					num++;
					if(tempedge->dest != i){
						dest = tempedge->dest;
						lenb = redGraph->array[dest].len;
					}
					tempedge = tempedge->next;
				}
				if(num == 2 && dest != -1){
					if(lenb >= lena){
						// collapse dest and i in, than delete i
//						printf("Collapse upSpur: %i --> %i\n",i,dest);
						collapseEdges(i,dest,downedge->dest,0);
					}
					else{
						destNodeNum = 0;
						destNodeList = setDestList(i,dest,lena,lenb,downedge,1,destNodeList,&destNodeNum);
						if(destNodeNum != -1) {
//							printf("Nodes: %i\n",destNodeNum);
							collapseEdges_strong(i,dest,0,destNodeList,destNodeNum);
							deleteSpur(i,0);
//							exit(1);
						}
					}
				}
				else if(num>2){
					// Search for the true destination in all kids of ori whose have more further children
					num = 0;
					tempedge = redGraph->array[downedge->dest].head;
					while(tempedge){
						if(redGraph->array[tempedge->dest].head){
							dest = tempedge->dest;
							lenb = redGraph->array[dest].len;
							num++;
						}
						tempedge = tempedge->next;
					}
					if(num==1){
						if(lenb >= lena){
//							printf("Collapse downSpur (more than one sibling): %i --> %i\n",i,dest);
	//						collapseEdges(i,dest,downedge->dest,0);
						}
						else{
							// collapse dest and i in, than delete i
//							printf("Collapse upSpur over more than one node(more than on sibling): %i --> %i\n",i,dest);
	//						collapseEdges_strong(i,dest,downedge->dest,0);
						}

					}
				}
			}
		}
	}
}
