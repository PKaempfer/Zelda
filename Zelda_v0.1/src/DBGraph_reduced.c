/*
 * DBGraph_reduced.c
 *
 *  Created on: Nov 14, 2014
 *      Author: kaempfpp
 */

#include "DBGraph_reduced.h"
#include <stdlib.h>
#include <stdio.h>
#include "DBGraph.h"
#include "uthash.h"
#include "../PKmiscLib/include/misc_hash.h"

int cNodeId;
struct redGraph* redGraph;
struct tempHash *temp = NULL; // for real time-check if unreduced junction is already in reduced graph;
struct Reads *allReads = NULL;

void collapseEdges(int i, int dest,int ori, int up){
	// First step -> Copy read information of the spur edge to the sibling edge
	// a is true --> b is spur
	// up -> nodes have a common parent; !up -> nodes have a common child
	struct ReadNode *anode, *abefnode,*bnode,*temp;
	int alen,blen,abdiff;
	struct KannteNode *kannte,*befkannte;

	anode = redGraph->array[dest].headread;
	abefnode = NULL;
	bnode = redGraph->array[i].headread;
	alen = redGraph->array[dest].len;
	blen = redGraph->array[i].len;
	if(up) abdiff = alen - blen;
	else abdiff = 0;

//	printf("(%i to %i) -> up: %i (diff: %i)\n",i,dest,up,abdiff);

	// a-node is empty
	if(!anode){
		redGraph->array[dest].headread = redGraph->array[i].headread;
		redGraph->array[i].headread = NULL;
		while(bnode){
//			printf("Change bnode pos while anode was empty (%i -> %i)\n",i,dest);
			kannte = bnode->read->headkannte;
			while(kannte){
				if(kannte->dest == i && kannte->pos == bnode->pos){
					printf("Change ReadInfo at the edge\n");
					printf("ID: %i (dir: %i)---> %i -> %i (Dest: %i -> %i, pos: %i -> %i) Up: %i\n",bnode->read->ID,bnode->dir,i,dest,kannte->dest,dest,kannte->pos,bnode->pos+abdiff,up);
					kannte->dest = dest;
					kannte->pos = bnode->pos+abdiff;
				}
				kannte = kannte->next;
			}
			bnode->pos+=abdiff;
			bnode = bnode->next;
		}
	}
	else{
		while(anode){
			while(bnode){
				if(bnode->pos+abdiff >= anode->pos || !anode->next){
//					if(dest==2 && abefnode) printf("Vor Transfer anode befa %i a %i  b: %i (%i -> %i)\n",abefnode->read->ID,anode->read->ID, bnode->read->ID ,anode->pos,bnode->pos);
//					else if(dest==2) printf("Vor Transfer anode befa NULL a %i  b: %i (%i -> %i)\n",anode->read->ID, bnode->read->ID ,anode->pos,bnode->pos);
					redGraph->array[i].headread = bnode->next;
//					if (dest == 2) printf("abdiff: %i\n",abdiff);
					if(bnode->pos+abdiff >= anode->pos){
						bnode->next = anode;
						if(!abefnode){
							redGraph->array[dest].headread = bnode;
							abefnode = bnode;
						}
						else{
							abefnode->next = bnode;
							abefnode = abefnode->next;
						}
					}
					else{
						temp = anode->next;
						anode->next = bnode;
						bnode->next = temp;
						if(!abefnode){
							abefnode = anode;
							anode = bnode;
						}
						else{
//							if (dest == 2) printf("abdiff: %i\n",abdiff);
							abefnode = bnode;
							abefnode = anode;
						}
					}


					kannte = bnode->read->headkannte;
					while(kannte){
						if(kannte->dest == i && kannte->pos == bnode->pos){
//							printf("Change ReadInfo at the edge\n");
//							printf("ID: %i (dir: %i)---> %i -> %i (Dest: %i -> %i, pos: %i -> %i) Up: %i\n",bnode->read->ID,bnode->dir,i,dest,kannte->dest,dest,kannte->pos,bnode->pos+abdiff,up);
							kannte->dest = dest;
							kannte->pos = bnode->pos+abdiff;
						}
						kannte = kannte->next;
					}
					bnode->pos = bnode->pos+abdiff;
					bnode = redGraph->array[i].headread;
				}
				else{
					break;
				}
			}
			abefnode = anode;
			anode = anode->next;
		}
	}

	// Second step -> Deletion of the edge (if there are more than one spur, don't cut the other once)
	struct edge *tempedge,*temptemp;
	if(up) tempedge = redGraph->array[ori].tail;
	else tempedge = redGraph->array[ori].head;

	if(up){
		free(redGraph->array[i].head);
		redGraph->array[i].head = NULL;
	}
	else{
		free(redGraph->array[i].tail);
		redGraph->array[i].tail = NULL;
	}
	if(tempedge->dest == i){ // to delete edge (i) is on first position
		if (up){
			redGraph->array[ori].tail = tempedge->next;
//			free(tempedge);
		}
		else{
			redGraph->array[ori].head = tempedge->next;
//			free(tempedge);
		}
		free(tempedge);
	}
	else{ // to delete edge (i) is not on first position
		while(tempedge->next){
			if(tempedge->next->dest == i){
				temptemp = tempedge->next->next;
				free(tempedge->next);
				tempedge->next = temptemp;
				return;
			}
			else tempedge = tempedge->next;
		}
//		if(up){
//			free(redGraph->array[ori].tail->next);
//			redGraph->array[ori].tail->next = NULL;
//		}
//		else{
//			free(redGraph->array[ori].head->next);
//			redGraph->array[ori].head->next = NULL;
//		}
	}
}

void reduceRedGraph(){
	int i;
	int num, dest;
	int lena=0, lenb=0;
	struct edge *downedge;
	struct edge *upedge;
	struct edge *tempedge;

	for(i=1;i<=redGraph->V;i++){
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
				if(num == 2 && dest != -1 && lenb >= lena){
					// collapse dest and i in
//					printf("Collapse downSpur: %i --> %i\n",i,dest);
					collapseEdges(i,dest,upedge->dest,1);
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
					if(num==1 && lenb >= lena){
//						printf("Collapse downSpur (more than on sibling): %i --> %i\n",i,dest);
						collapseEdges(i,dest,upedge->dest,1);
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
//				printf("Potential Upspur: %i --> %i\n",i,dest);
				if(num == 2 && dest != -1 && lenb >= lena){
					// collapse dest and i in, than delete i
//					printf("Collapse upSpur: %i --> %i\n",i,dest);
					collapseEdges(i,dest,downedge->dest,0);
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
					if(num==1 && lenb >= lena){
//						printf("Collapse downSpur (more than one sibling): %i --> %i\n",i,dest);
						collapseEdges(i,dest,downedge->dest,0);
					}
				}
			}
		}
	}
}

// RUCURSIVE FUNCTION (goes upwards unique paths)
void internalVerticalReduction(int c, int p){
	// c -> child, p -> parent
	int degree = 0;
//	printf(" -> %i",p);

	// Collapse the edges

	// 1. Erhöhe pos aller kannten von p um len(c) und gehe bis zum letzten, gehe alle reads durch und verändere die pos
	//    aus Kannte p durch erhöhen um len(c) und zeiger von p auf c --> Done

	struct KannteNode *readkannte;
	int pos;
	degree = redGraph->array[c].len;
	struct ReadNode *readnode = redGraph->array[p].headread;
	struct ReadNode *lastreadnode = NULL;
	while(readnode){
		pos = readnode->pos;
		readkannte = readnode->read->headkannte;
//		if(p==686) printf("Step one -> Readnode: %i\n",readnode->pos);
		while(readkannte){
			if(readkannte->dest == p && readkannte->pos == pos){
				readkannte->dest = c;
				readkannte->pos += degree;
			}
			readkannte = readkannte->next;
		}
//		if(p==686) printf("Step two\n");

		readnode->pos += degree;
		if(!readnode->next) lastreadnode = readnode;
		readnode = readnode->next;
	}
	redGraph->array[c].len += redGraph->array[p].len;

	// 2. Hänge alle Kannten von c an das ende von p
	if(lastreadnode) lastreadnode->next = redGraph->array[c].headread;
	else redGraph->array[p].headread = redGraph->array[c].headread;

	// 3. Hänge die headkannte von p and die headkannte von c und NULLe p-Headkannte
	redGraph->array[c].headread = redGraph->array[p].headread;
	redGraph->array[p].headread = NULL;

	// 4. Search for further reducable parent node
	struct edge *kedge;
	int nextp = -1;
	degree = 0;
	kedge = redGraph->array[p].head;
	while(kedge){
		degree++;
		nextp = kedge->dest;
		kedge = kedge->next;
	}
	if(degree == 1){
		// Kann ich mir schenken, auf grunde der nächsten abfrage weiss ich die antwort ja bereits
		kedge = redGraph->array[p].tail;
		while(kedge){
			degree++;
			kedge = kedge->next;
		}
		if(degree == 2){
			kedge = redGraph->array[nextp].tail;
			while(kedge){
				degree++;
				kedge = kedge->next;
			}
			if(degree==3) degree=-1;
		}
	}

	// 5. Delete edge and update child and parent nodes (c gets all parents of p, all parents of p gets c as child)
	free(redGraph->array[c].head);
	free(redGraph->array[p].tail);
	redGraph->array[c].head = NULL;
	redGraph->array[p].tail = NULL;

	struct edge *paredge;
	if(nextp != -1){
		paredge = redGraph->array[p].head;
		while(paredge){
			kedge = redGraph->array[paredge->dest].tail;
			while(kedge){
				if(kedge->dest == p) kedge->dest = c;
				kedge = kedge->next;
			}
			paredge = paredge->next;
		}
		kedge = redGraph->array[p].head;
		while(kedge){
			redGraph->array[p].head = kedge->next;
			kedge->next = redGraph->array[c].head;
			redGraph->array[c].head = kedge;
			kedge = redGraph->array[p].head;
		}
//		redGraph->array[c].head->dest = nextp;
	}
	else{
//		printf("\nProblemchen");
		return;
	}


	// Stuff to print the list for the moment -> can be deleted
//	int i=c;
//	struct edge *edge = kedge;
//	struct KannteNode *kannte;
//	printf("\n Adjacency list of vertex %d (len: %i / depth: %i)\n head ", i,redGraph->array[i].len,redGraph->array[i].depth);
//	edge = redGraph->array[i].head;
//	while(edge){
//	       printf("-> %i ", edge->dest);
//		edge = edge->next;
//	}
//	printf("\n tail ");
//	edge = redGraph->array[i].tail;
//	while(edge){
//	       printf("-> %i ", edge->dest);
//		edge = edge->next;
//	}
//	printf("\n reads ");
//	readnode = redGraph->array[i].headread;
//	while(readnode){
//		printf("\n -> %i (%i / %i)",readnode->read->ID,readnode->pos,readnode->dir);
//		kannte = readnode->read->headkannte;
//		printf(" other nodes containing this read: ");
//		while(kannte){
//			printf(" %i (%i)  ",kannte->dest, kannte->pos);
//			kannte = kannte->next;
//		}
//
//		readnode = readnode->next;
//	}
//	printf("\n");


	if(degree == -1) internalVerticalReduction(c,nextp);

}

/**
 * Seeks for unique paths between the nodes after the horizontal reduction and merge them to single nodes in the reduced graph. Iterative simplification of the graph
 * @return 1 - If the graph was reduces and a new iteration to find a point of horizontal reduction makes sense
 * 		   0 - Otherwise: Nothing could be reduced, go ahead without a new reduction iteration.
 */
int verticalReduction(){
	int succ = 0, inter = 0;
	int i;
	int dest;
	int indegree, outdegree;
	struct edge *edge;
	for (i = 1; i <= redGraph->V; i++){
//		printf("Node: %i\n",i);
	 	indegree = 0;
	 	outdegree = 0;
	 	edge = redGraph->array[i].tail;
	 	while(edge){
	 		dest = edge->dest;
	 		outdegree++;
	 		edge = edge->next;
	 	}
	 	if(outdegree != 1){
	 		// Hypothetical start point of a reduction and not an upper end point
	 		if(redGraph->array[i].head){
	 			inter = 1;
//	 			printf("Hypothetical start point for vertical reduction at node: %i (outdegree: %i)\n",i,outdegree);
	 		}
	 	}
	 	else{
	 		//  Node has one child, but the child has many parents
	 		edge = redGraph->array[dest].head;
	 		while(edge){
	 			indegree++;
	 			edge = edge->next;
	 		}
	 		if(indegree > 1){
	 			inter = 1;
	 		}
	 	}
	 	if(inter){
//	 		printf("Test: %i\n",i);
		 	edge = redGraph->array[i].head;
		 	indegree = 0;
		 	while(edge){
		 		dest = edge->dest;
		 		indegree++;
		 		edge = edge->next;
		 	}
		 	if(indegree == 1){
//		 		printf("Test 2 : %i\n",i);
		 		// Hypothetical startpoint of a reducable edge
		 		outdegree = 0;
		 		edge = redGraph->array[dest].tail;
		 		while(edge){
	 				outdegree++;
	 				edge = edge->next;
//	 				printf("Outdegree: %i\n",outdegree);
	 			}
	 			if(outdegree == 1){
//	 				printf("verticalReduction start node: %i\n", i);
	 				internalVerticalReduction(i,dest);
//	 				printf("\n");
	 				succ = 1;
	 			}
	 		}
//		 	printf("Test 3!!!\n");
	 	}
	 	inter = 0;
	}
	if(succ) return 1;
	else return 0;
}

int scaffoldRedGraph(){
	printf("CHECKPOINT: Scaffold RedGraph\n");
	struct edge *edge, *edgehead;
	struct ReadNode *readNode, *readdestNode;
	struct KannteNode *kannte;
//	struct Reads *read;
	int i,j;
	int dest; // A free node in the graph, can be used for node duplication
	int ori;  // Repeat node, which should be duplicated

	for(i=1;i<=redGraph->V;i++){

		if(redGraph->array[i].head && redGraph->array[i].tail){
			if(!redGraph->array[i].head->next && !redGraph->array[i].tail->next){
				if(redGraph->array[i].head->dest == redGraph->array[i].tail->dest){
					ori = redGraph->array[i].head->dest;
					for(j=1;j<redGraph->V;j++){
						if(!redGraph->array[j].head && !redGraph->array[j].tail && !redGraph->array[j].headread){
							dest = j;
							break;
						}
					}
					printf("Node %i is a unique sequence between the repeat node %i\n",i,redGraph->array[i].tail->dest);
					printf("This could be easily resolved by copy information to node %i\n",dest);

					// 1. Set all reads to both nodes
					// Distribute the reads on ori and dest nodes, regarding to the direction they come from, if known
					// To both , otherwise. If PE information is already provided, Distribute in regard to the direction of the second pair
					readNode = redGraph->array[ori].headread;
					while(readNode){
						struct ReadNode *newreadNode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
						newreadNode->dir = readNode->dir;
						newreadNode->pos = readNode->pos;
						newreadNode->flag = 0;
						newreadNode->read = readNode->read;
						newreadNode->next = NULL;
						kannte = newreadNode->read->headkannte;
						while(kannte){
							if(kannte->dest == ori && readNode->pos == kannte->pos){
								struct KannteNode *newKannte = (struct KannteNode*)malloc(sizeof(struct KannteNode));
								newKannte->dest = dest;
								newKannte->pos = kannte->pos;
								newKannte->next = kannte->next;
								newKannte->ReadNode = newreadNode;
								kannte->next = newKannte;
							}
							kannte = kannte->next;
						}

						if(!redGraph->array[dest].headread){
							redGraph->array[dest].headread = newreadNode;
							readdestNode = newreadNode;
						}
						else{
							readdestNode->next = newreadNode;
							readdestNode = readdestNode->next;
						}
						readNode = readNode->next;
					}
					redGraph->array[dest].len = redGraph->array[ori].len;

//					redGraph->array[dest].headread = redGraph->array[ori].headread;
//					redGraph->array[dest].len = redGraph->array[ori].len;
//					readNode = redGraph->array[dest].headread;
//					while(readNode){
//						ari = 0;
//						kannte = readNode->read->headkannte;
//						while(kannte){
//							if(kannte->dest == dest) ari = 1;
//							kannte = kannte->next;
//						}
//						if(!ari){
//							kannte = readNode->read->headkannte;
//							while(kannte){
//								if(kannte->dest == ori){
//									struct KannteNode *newKannte = (struct KannteNode*)malloc(sizeof(struct KannteNode));
//									newKannte->dest = dest;
//									newKannte->pos = kannte->pos;
//									newKannte->next = kannte->next;
//									kannte->next = newKannte;
//								}
//								kannte = kannte->next;
//							}
//
//						}
//						readNode = readNode->next;
//					}

					// 2. Set i's tail edge from ori to dest
					redGraph->array[i].tail->dest = dest;

					// 3. Set all ori's tail heads from ori to dest
					edge = redGraph->array[ori].tail;
					while(edge){
						if(edge->dest != i){
							edgehead = redGraph->array[edge->dest].head;
							while(edgehead){
								if(edgehead->dest == ori) edgehead->dest = dest;
								edgehead = edgehead->next;
							}
						}
						edge = edge->next;
					}

					// 4. Delete ori's head edge to i, Set dest's head edge to i
					edge = redGraph->array[ori].head;
					if(edge->next)	printf("4. EdgeDest: %i\n",edge->dest);
					if(edge->dest == i){
						redGraph->array[dest].head = edge;
						redGraph->array[ori].head = edge->next;
						edge->next = NULL;
						edgehead = edge;
					}
					else{
						while(edge->next){
							printf("While\n");
							if(edge->next->dest == i){
								edgehead = edge->next;
								edgehead->next = NULL;
								edge->next = edge->next->next;
								redGraph->array[dest].head = edgehead;
								// break;
							}
							else edge = edge->next;
						}
					}

					// 5. All tail edges of ori, except to i now as new tail edges of dest
					edge = redGraph->array[ori].tail;
					while(edge && edge->dest!=i){
						edgehead = edge;
						redGraph->array[ori].tail = edge->next;
						edge = edge->next;
						edgehead->next = redGraph->array[dest].tail;
						redGraph->array[dest].tail = edgehead;
					}
					if(edge){
						while(edge->next){
							if(edge->next->dest!=i){
								edgehead = edge->next;
								edgehead->next = redGraph->array[dest].tail;
								redGraph->array[dest].tail = edgehead;
								edge->next = edge->next->next;
							}
						}
					}
				}
			}
		}
	}


	return 1;
}

int outDegree(int node){
	struct AdjListNode *s = graph->array[node].head;
	int num = 0;
	while(s){
		num++;
		s = s->next;
	}
	return num;
}


void redGraphConnector(){
	int i;
	struct ReadNode* readnode;
	struct readEdge* readedge;
	struct KannteNode* kannte;
	for(i=1; i<=redGraph->V;i++){
		redGraph->array[i].readedge = NULL;
		readedge = redGraph->array[i].readedge;
		readnode = redGraph->array[i].headread;
		while(readnode){
			if(readnode->pos <= maxRlen){
				kannte = readnode->read->headkannte;
				while(kannte){
					if(kannte->dest != i){
						if(!readedge){
							readedge = (struct readEdge*)malloc(sizeof(struct readEdge));
							readedge->KannteID = kannte->dest;
							readedge->next = redGraph->array[i].readedge;
							redGraph->array[i].readedge = readedge;
						}
						else{
							if(readedge->KannteID != kannte->dest){
								readedge = redGraph->array[i].readedge;
								while(readedge){
									if(readedge->KannteID == kannte->dest){
										break;
									}
									readedge = readedge->next;
								}
								if(!readedge){
									readedge = (struct readEdge*)malloc(sizeof(struct readEdge));
									readedge->KannteID = kannte->dest;
									readedge->next = redGraph->array[i].readedge;
									redGraph->array[i].readedge = readedge;
								}
							}
						}
						break;
					}
					kannte = kannte->next;
				}
			}
			readnode = readnode->next;
		}
	}
}

/*	Writes the reduced graph as a list to console for checking and debugging */

//void printRedGraph_MaskCon(struct myovlList *ovlGraph){
//	int i;
//	struct edge *edge;
//	struct ReadNode *readnode;
//	struct KannteNode *kannte;
//
//	for(i=0;i<redGraph->V;i++){
////		if((i > 50 && i < 58) || ( i>135 && i<139)){ 		//// For testing
////		if(i==2){
//			if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
//				printf("\n Adjacency list of vertex %d (len: %i)\n head ", i,redGraph->array[i].len);
//				edge = redGraph->array[i].head;
//				while(edge){
//				    	   printf("-> %i ", edge->dest);
//					edge = edge->next;
//				}
//				printf("\n tail ");
//				edge = redGraph->array[i].tail;
//				while(edge){
//				       printf("-> %i ", edge->dest);
//					edge = edge->next;
//				}
//				printf("\n reads ");
//				readnode = redGraph->array[i].headread;
//				while(readnode){
//					if(readnode->flag)
//					printf("\n -> %i (%i / %i)",readnode->read->ID,readnode->pos,readnode->dir);
//					kannte = readnode->read->headkannte;
//					printf(" other nodes containing this read: ");
//					while(kannte){
//						printf(" %i (%i)  ",kannte->dest, kannte->pos);
//						kannte = kannte->next;
//					}
//
//					readnode = readnode->next;
//				}
//				printf("\n");
//			}
////		}
//
//	}
//	printf("\n");
//}

void printRedGraphToFile(char* filePath){
	FILE* graph = fopen(filePath,"w");

	int i;
	struct edge* edge;
	struct ReadNode* readnode;
	struct readEdge* readedge;
	struct KannteNode* kannte;

	for(i=1;i<=redGraph->V;i++){
//		if((i > 50 && i < 58) || ( i>135 && i<139)){ 		//// For testing
//		if(i==2){
			if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
				fprintf(graph,"\n Adjacency list of vertex %d (len: %i)\n head ", i,redGraph->array[i].len);
				edge = redGraph->array[i].head;
				while(edge){
					fprintf(graph,"-> %i ", edge->dest);
					edge = edge->next;
				}
				fprintf(graph,"\n tail ");
				edge = redGraph->array[i].tail;
				while(edge){
					fprintf(graph,"-> %i ", edge->dest);
					edge = edge->next;
				}
				fprintf(graph,"\n ReadEdges: ");
				readedge = redGraph->array[i].readedge;
				while(readedge){
					fprintf(graph,"-> %i ", readedge->KannteID);
					readedge = readedge->next;
				}
				fprintf(graph,"\n reads ");
				readnode = redGraph->array[i].headread;
				while(readnode){
					fprintf(graph,"\n -> %i (%i / %i)",readnode->read->ID,readnode->pos,readnode->dir);
					kannte = readnode->read->headkannte;
					fprintf(graph," other nodes containing this read: ");
					while(kannte){
						if(kannte->pos <= maxRlen)
							fprintf(graph," %i (%i)  ",kannte->dest, kannte->pos);
						kannte = kannte->next;
					}

					readnode = readnode->next;
				}
				fprintf(graph,"\n");
			}
//		}

	}
	fprintf(graph,"\n");


	fclose(graph);
}

void printRedGraph(){
	int i;
	struct edge *edge;
	struct ReadNode *readnode;
	struct KannteNode *kannte;

	for(i=1;i<=redGraph->V;i++){
//		if((i > 50 && i < 58) || ( i>135 && i<139)){ 		//// For testing
//		if(i==2){
			if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
				printf("\n Adjacency list of vertex %d (len: %i)\n head ", i,redGraph->array[i].len);
				edge = redGraph->array[i].head;
				while(edge){
				       printf("-> %i ", edge->dest);
					edge = edge->next;
				}
				printf("\n tail ");
				edge = redGraph->array[i].tail;
				while(edge){
				       printf("-> %i ", edge->dest);
					edge = edge->next;
				}
				printf("\n reads ");
				readnode = redGraph->array[i].headread;
				while(readnode){
					printf("\n -> %i (%i / %i)",readnode->read->ID,readnode->pos,readnode->dir);
					kannte = readnode->read->headkannte;
					printf(" other nodes containing this read: ");
					while(kannte){
						printf(" %i (%i)  ",kannte->dest, kannte->pos);
						kannte = kannte->next;
					}

					readnode = readnode->next;
				}
				printf("\n");
			}
//		}

	}
	printf("\n");
}

void printRedGraphEdge(){
	int i;
	struct edge *edge;

	for(i=0;i<redGraph->V;i++){
//		if((i > 50 && i < 58) || ( i>135 && i<139)){ 		//// For testing
//		if(i==80){
			if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
				printf("\n Adjacency list of vertex %d (len: %i)\n head ", i,redGraph->array[i].len);
				edge = redGraph->array[i].head;
				while(edge){
				       printf("-> %i ", edge->dest);
					edge = edge->next;
				}
				printf("\n tail ");
				edge = redGraph->array[i].tail;
				while(edge){
				       printf("-> %i ", edge->dest);
					edge = edge->next;
				}
				printf("\n");
			}
//		}

	}
	printf("\n");
}

void printRedGraphList(char *filename){
	FILE *listFile = fopen(filename,"w");
	int i;
	struct edge *edge;
	struct ReadNode *readnode;
	struct KannteNode *kannte;

	for(i=0;i<redGraph->V;i++){
//		if((i > 50 && i < 58) || ( i>135 && i<139)){ 		//// For testing
//		if(i==80){
			if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
				fprintf(listFile,"\n Adjacency list of vertex %d (len: %i)\n head ", i,redGraph->array[i].len);
				edge = redGraph->array[i].head;
				while(edge){
				       fprintf(listFile,"-> %i ", edge->dest);
					edge = edge->next;
				}
				fprintf(listFile,"\n tail ");
				edge = redGraph->array[i].tail;
				while(edge){
				       fprintf(listFile,"-> %i ", edge->dest);
					edge = edge->next;
				}
				fprintf(listFile,"\n reads ");
				readnode = redGraph->array[i].headread;
				while(readnode){
					fprintf(listFile,"\n -> %i (%i / %i)",readnode->read->ID,readnode->pos,readnode->dir);
					kannte = readnode->read->headkannte;
					fprintf(listFile," other nodes containing this read: ");
					while(kannte){
						fprintf(listFile," %i (%i)  ",kannte->dest, kannte->pos);
						kannte = kannte->next;
					}

					readnode = readnode->next;
				}
				fprintf(listFile,"\n");
			}
//		}
	}
	fprintf(listFile,"\n");
	fclose(listFile);
}

/*	Writes the reduces Graph to dot-file
	Can be opend by xdot or GraphViz 	*/

void printRedDot(char *fileName){
	int i;
	FILE *redDot = fopen(fileName,"w");

	struct edge *edge;
	struct ReadNode *read;
	int readcount;

	fprintf(redDot,"digraph redDBG {\n");
	for(i=1;i<=redGraph->V;i++){

		readcount = 0;
		read = redGraph->array[i].headread;
		while(read){
			readcount++;
			read = read->next;
		}

		if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
			fprintf(redDot,"%i [label=\"%i (%i/%i)\"];\n",i,i,redGraph->array[i].len,readcount);
			edge = redGraph->array[i].tail;
			while(edge){
				fprintf(redDot,"%i -> %i;\n",i,edge->dest);
				edge = edge->next;
			}
		}
//		else printf("Empty Node\n");
	}
	fprintf(redDot,"}\n");
	fclose(redDot);
	printf("Wrote dot file: %s\n",fileName);
}

///* 	Called by travToRedOvl
//	Transfers the graph from the error corrected original deBruijn graph,
//	represented by adjacency list in reduced one.
//	Traversal follows the path till it reaches a vertex with an in- or outgree unequal to one
//	Than it stops and creates a new reduces edge (Kannte) and give it an new index if it doesn't already exists */
//
void collectInternalNodes(int i,int Knoten){
	int j,outdegree;
int goahead = 0;
	int depth = 1;
	struct AdjListNode *parent,*child;
	struct Kannte *s = &redGraph->array[Knoten];
	struct tempHash *p;

	struct LinkListNode *link;
	struct Reads *read;
	struct KannteNode *kanntennode;
	int key,dir;

	do{
		goahead = 0;
		parent = graph->array[i].head;
		child = graph->array[i].tail;

		if(((parent && !parent->next) && (child && !child->next)) || !parent){
			link = graph->array[i].Link;
			while(link){
				if(__builtin_clz(link->ID)){	// start
					key = link->ID;
					dir = 0;
				}
				else{							// end
					key = link->ID & DEL_READ_END;
					dir = 1;
				}
				HASH_FIND(hhb,allReads,&key,sizeof(int),read);
				if(!read){
					read = (struct Reads*)malloc(sizeof(struct Reads));
					read->ID = key;
					read->headkannte = NULL;
					HASH_ADD(hhb,allReads,ID,sizeof(int),read);
				}

				kanntennode = (struct KannteNode*)malloc(sizeof(struct KannteNode));
				kanntennode->dest = Knoten;
				kanntennode->pos = depth;
				kanntennode->next = read->headkannte;
				read->headkannte = kanntennode;

				struct ReadNode *readnode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
				kanntennode->ReadNode = readnode;
				readnode->dir = dir;
				readnode->pos = depth;
				readnode->flag = 0;
				readnode->read = read;
				readnode->next = redGraph->array[Knoten].headread;
				redGraph->array[Knoten].headread = readnode;

				link = link->next;
			}

			if((parent && !parent->next) && (child && !child->next)){
				depth++;
				i = parent->dest;
				goahead=1;
			}
		}
	} while(goahead);
	s->len = depth;

	// connect the new node to the parents and the parents to the new node
//		printf("P: %i\n",i);
	HASH_FIND( hhb ,temp, &i, sizeof(int), p);
	if(p){
		for(j=0;j<p->numIDs;j++){
			struct edge* newedge = (struct edge*)malloc(sizeof(struct edge));
			newedge->dest = Knoten;
			newedge->next = redGraph->array[p->newID[j]].tail;
			redGraph->array[p->newID[j]].tail = newedge;

			newedge = (struct edge*)malloc(sizeof(struct edge));
			newedge->dest = p->newID[j];
			newedge->next = s->head;
			s->head = newedge;
		}
	}
	else{
		outdegree = outDegree(i);
		if(outdegree){
			p = (struct tempHash*)malloc(sizeof(struct tempHash));
			p->oldID = i;
			p->numIDs = outdegree;
			p->newID = (int*)malloc(sizeof(int)*(outdegree+1));
			for(j=0;j<outdegree;j++){
//				redGraph->array[cNodeId] = (struct Knoten*)malloc(sizeof(struct Knoten));
				p->newID[j] = cNodeId++;
			}
			HASH_ADD(hhb,temp,oldID,sizeof(int),p);
			for(j=0;j<p->numIDs;j++){
				struct edge* newedge = (struct edge*)malloc(sizeof(struct edge));
				newedge->dest = Knoten;
				newedge->next = redGraph->array[p->newID[j]].tail;
				redGraph->array[p->newID[j]].tail = newedge;

				newedge = (struct edge*)malloc(sizeof(struct edge));
				newedge->dest = p->newID[j];
				newedge->next = s->head;
				s->head = newedge;
			}
		}
	}
		// whats next, create next (empty) edges, save them
}

/* 	Traverse once through the original, corrected graph and searches for junctions.
	If one is found, it creates a new reduces edge and calls collectInternalNodes
	till the path is not reducible anymore. */

void travToRedOVL(){
	int i,j,outdegree;
	cNodeId = 0;

	// create redGraph
	redGraph = (struct redGraph*)malloc(sizeof(struct redGraph));
	redGraph->array = (struct Kannte*)malloc(sizeof(struct Kannte)*(graph->V/10));
	redGraph->vFlag = (char*)malloc(sizeof(char)*(graph->V/10));
	struct Kannte *newArray;
	struct KannteNode *kanntennode;
	struct LinkListNode *link;
	struct Reads *read;
	int key,dir;
	char *newFlags;
	redGraph->V = (graph->V)/10;
	for(j=0;j<redGraph->V;j++){
		redGraph->array[j].len = 0;
		redGraph->array[j].head = NULL;
		redGraph->array[j].tail = NULL;
		redGraph->array[j].headread = NULL;
		redGraph->array[j].maxpos = -1;
		redGraph->array[j].minpos = -1;
		redGraph->vFlag[j]=0;
	}


	struct tempHash *s;
	struct AdjListNode *chnode,*panode;
	for(i=1; i<graph->V; i++){
		chnode = graph->array[i].tail;
		panode = graph->array[i].head;
		if((chnode && !chnode->next) && (panode && !panode->next)){
			// in and out == 1
			// internal node, catched by collectInternalNodes
			continue;
		}
		else if(!panode){											// out = 0
			// parent-dead-end (upstream)
			continue;
		}
		else{
			// Junction or child-dead-end
//			printf("Junction at: %i\n",i);
			HASH_FIND( hhb ,temp, &i, sizeof(int), s);
			if(s==NULL){	// New junction, not inndexed yet
				outdegree = outDegree(i);
				if(outdegree){
					s = (struct tempHash*)malloc(sizeof(struct tempHash));
					s->oldID = i;
					s->numIDs = outdegree;
					s->newID = (int*)malloc(sizeof(int)*(outdegree+1));
					for(j=0;j<outdegree;j++){
						s->newID[j] = cNodeId++;
					}
					HASH_ADD(hhb,temp,oldID,sizeof(int),s);
					for(j=0;j<s->numIDs;j++){
						link = graph->array[i].Link;
						while(link){
							if(__builtin_clz(link->ID)){	// start
								key = link->ID;
								dir = 0;
							}
							else{							// end
								key = link->ID & DEL_READ_END;
								dir = 1;
							}
							HASH_FIND(hhb,allReads,&key,sizeof(int),read);
							if(!read){
								read = (struct Reads*)malloc(sizeof(struct Reads));
								read->ID = key;
								read->headkannte = NULL;
								HASH_ADD(hhb,allReads,ID,sizeof(int),read);
							}
							kanntennode = (struct KannteNode*)malloc(sizeof(struct KannteNode));
							kanntennode->dest = s->newID[j];
							kanntennode->pos = 0;
							kanntennode->next = read->headkannte;
							read->headkannte = kanntennode;

							struct ReadNode *readnode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
							kanntennode->ReadNode = readnode;
							readnode->dir = dir;
							readnode->pos = 0;
							readnode->flag = 0;
							readnode->read = read;
							readnode->next = redGraph->array[s->newID[j]].headread;
							redGraph->array[s->newID[j]].headread = readnode;

							link = link->next;

						}

						// Collect information of the first node
//						printf("Start internal node collection at node: %i\n",i);
						collectInternalNodes(panode->dest,s->newID[j]);
						panode = panode->next;
					}
				}
			}
			else{ // Junction is already indexed
				for(j=0;j<s->numIDs;j++){
					link = graph->array[i].Link;
					while(link){
						if(__builtin_clz(link->ID)){	// start
							key = link->ID;
							dir = 0;
						}
						else{							// end
							key = link->ID & DEL_READ_END;
							dir = 1;
						}
						HASH_FIND(hhb,allReads,&key,sizeof(int),read);
						if(!read){
							read = (struct Reads*)malloc(sizeof(struct Reads));
							read->ID = key;
							read->headkannte = NULL;
							HASH_ADD(hhb,allReads,ID,sizeof(int),read);
						}
						kanntennode = (struct KannteNode*)malloc(sizeof(struct KannteNode));
						kanntennode->dest = s->newID[j];
						kanntennode->pos = 0;
						kanntennode->next = read->headkannte;
						read->headkannte = kanntennode;

						struct ReadNode *readnode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
						kanntennode->ReadNode = readnode;
						readnode->dir = dir;
						readnode->pos = 0;
						readnode->flag = 0;
						readnode->read = read;
						readnode->next = redGraph->array[s->newID[j]].headread;
						redGraph->array[s->newID[j]].headread = readnode;

						link = link->next;
					}
					// Collect information of the first node
//					printf("Start internal node collection at node: %i\n",i);
					collectInternalNodes(panode->dest,s->newID[j]);
					panode = panode->next;
				}
			}
		}

		if(cNodeId > redGraph->V-10){ // Resize reduced graph (Size depends hardly by k-size and the percentage of error)
			printf("Reallac redGraph Size (old: %i)\n",redGraph->V);
			newArray = (struct Kannte*)realloc(redGraph->array,(sizeof(struct Kannte))*(2*redGraph->V));
			newFlags = (char*)realloc(redGraph->vFlag,sizeof(char)*(2*redGraph->V));
			redGraph->array = newArray;
			redGraph->vFlag = newFlags;
			for(j=redGraph->V;j<(redGraph->V*2);j++){
				redGraph->array[j].len = 0;
				redGraph->array[j].head = NULL;
				redGraph->array[j].tail = NULL;
				redGraph->array[j].headread = NULL;
				redGraph->vFlag[j]=0;
			}
			redGraph->V*=2;
		}
	}
	redGraph->V = cNodeId;
	printf("End Graph-contraction. Remaining Nodes: %i\n",redGraph->V);

	// Delete and free the non-reduced Graph

}

///* 	Called by travToRedOvl
//	Transfers the graph from the error corrected original deBruijn graph,
//	represented by adjacency list in reduced one.
//	Traversal follows the path till it reaches a vertex with an in- or outgree unequal to one
//	Than it stops and creates a new reduces edge (Kannte) and give it an new index if it doesn't already exists */
//
int collectInternalNodes_v2(int i,int Knoten){
	int goahead = 1;
	int depth = 0;
	struct AdjListNode *parent;

	struct LinkListNode *link;
	struct Reads *read;
	struct KannteNode *kanntennode;
	int key,dir;

	do{
		parent = graph->array[i].head;
		if(!parent){
			goahead = 0;
		}
		else if(parent && parent->next){
			goahead = 0;
		}
		else{
			if(graph->array[parent->dest].tail && graph->array[parent->dest].tail->next) goahead = 0;
		}
		if(goahead){
			i = parent->dest;
			depth++;
		}
		else{
			redGraph->array[Knoten].len = depth+1;
			return i;
		}
		link = graph->array[i].Link;
		while(link){
			if(__builtin_clz(link->ID)){	// start
				key = link->ID;
				dir = 0;
			}
			else{							// end
				key = link->ID & DEL_READ_END;
				dir = 1;
			}
			HASH_FIND(hhb,allReads,&key,sizeof(int),read);
			if(!read){
				read = (struct Reads*)malloc(sizeof(struct Reads));
				read->ID = key;
				read->headkannte = NULL;
				HASH_ADD(hhb,allReads,ID,sizeof(int),read);
			}

			kanntennode = (struct KannteNode*)malloc(sizeof(struct KannteNode));
			kanntennode->dest = Knoten;
			kanntennode->pos = depth;
			kanntennode->next = read->headkannte;
			read->headkannte = kanntennode;

			struct ReadNode *readnode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
			kanntennode->ReadNode = readnode;
			readnode->dir = dir;
			readnode->pos = depth;
			readnode->flag = 0;
			readnode->read = read;
			readnode->next = redGraph->array[Knoten].headread;
			redGraph->array[Knoten].headread = readnode;

			link = link->next;
		}
	} while(1);
	return 0;
}

static inline int countNodes(){
	int i;
	int countNodes = 0;
	struct AdjListNode *chnode;
	for(i=1;i<=graph->V;i++){
		if(graph->array[i].head || graph->array[i].tail){
			chnode = graph->array[i].tail;
			if(!chnode){
//				printf("No child\n");
				countNodes++;
			}
			else if((chnode && chnode->next)){
//				printf("More than one child\n");
				countNodes++;
			}
			else if(chnode && !chnode->next){
				// One Child with more than one parent
				if(graph->array[chnode->dest].head && graph->array[chnode->dest].head->next){
//					printf("One Child with more than one parent\n");
					countNodes++;
				}
			}
		}
	}
	return countNodes;
}

inline static void initRedNode(int j){
			redGraph->array[j].len = 0;
			redGraph->array[j].head = NULL;
			redGraph->array[j].tail = NULL;
			redGraph->array[j].headread = NULL;
			redGraph->array[j].maxpos = -1;
			redGraph->array[j].minpos = -1;
			redGraph->vFlag[j]=0;
}

/* 	Traverse once through the original, corrected graph and searches for junctions.
	If one is found, it creates a new reduces edge and calls collectInternalNodes
	till the path is not reducible anymore. */

/** ToDo:  Rewrite function completely: Traverse graph, count future redGraph nodes and add to hash table to switch ID from AdjaGraph to redGraph
 * 		   From each start connect to parents and parent to itself, than connect the internal nodes above, till the end condition is satisfied.
 * 		   Try to get rid of UT-Hash -> Memory Overkill
 * 		   StartCondition: 	(1) No Child
 * 		   					(2) More than one Child
 * 		   					(3) One Child, with more than one Parent
 * 		   EndConidtion:	(1) No Parent
 * 		   					(2) More than one Parent
 * 		   					(3) One Parent with more than one Child
 */
void travToRedOVL_v2(){
	int i,j,outdegree;
	cNodeId = 1;

	struct AdjListNode *chnode,*panode;
	int nodeNum = countNodes();
	int go;

	printf("Count Node: %i\n",nodeNum);
	printf("Create temporary hashTable32\n");
	struct hashTable32* ht = hashTable32_create(BITS_TO_REPRESENT(nodeNum));
	printf("Number of new nodes representing a start point of a reduced graph: %i (in bits: %i)\n",nodeNum,BITS_TO_REPRESENT(nodeNum));

	// create redGraph
	redGraph = (struct redGraph*)malloc(sizeof(struct redGraph));
	redGraph->array = (struct Kannte*)malloc(sizeof(struct Kannte)*(nodeNum+1)); // Empty 0-Element
	redGraph->vFlag = (char*)malloc(sizeof(char)*(nodeNum+1));
	redGraph->V = nodeNum;
	struct Kannte *newArray;
	struct KannteNode *kanntennode;
	struct LinkListNode *link;
	struct Reads *read;

	for(i=1;i<=graph->V;i++){
		if(graph->array[i].head || graph->array[i].tail){
			chnode = graph->array[i].tail;
			if(!chnode){
				hashTale32_insert(ht,i,cNodeId);
				initRedNode(cNodeId++);
			}
			else if((chnode && chnode->next)){
				hashTale32_insert(ht,i,cNodeId);
				initRedNode(cNodeId++);
			}
			else if(chnode && !chnode->next){
				// One Child with more than one parent
				if(graph->array[chnode->dest].head && graph->array[chnode->dest].head->next){
					hashTale32_insert(ht,i,cNodeId);
					initRedNode(cNodeId++);
				}
			}
		}
	}

	hashTable32_stats(ht);

	int key,dir;
	int redNode,redNodeParent;
	int iParent;
	struct edge* tempedge;
	for(i=1;i<=graph->V;i++){
		chnode = graph->array[i].tail;
		if(!chnode) go = 1;
		else if((chnode && chnode->next)) go = 1;
		else if(chnode && !chnode->next){
			// One Child with more than one parent
			if(graph->array[chnode->dest].head && graph->array[chnode->dest].head->next) go = 1;
			else continue;
		}
		if(go){
			redNode = hashTable32_find(ht,i);
			link = graph->array[i].Link;
			while(link){
				if(__builtin_clz(link->ID)){	// start
					key = link->ID;
					dir = 0;
				}
				else{							// end
					key = link->ID & DEL_READ_END;
					dir = 1;
				}
				HASH_FIND(hhb,allReads,&key,sizeof(int),read);
				if(!read){
					read = (struct Reads*)malloc(sizeof(struct Reads));
					read->ID = key;
					read->headkannte = NULL;
					HASH_ADD(hhb,allReads,ID,sizeof(int),read);
				}
				kanntennode = (struct KannteNode*)malloc(sizeof(struct KannteNode));
				kanntennode->dest = redNode;
				kanntennode->pos = 0;
				kanntennode->next = read->headkannte;
				read->headkannte = kanntennode;

				struct ReadNode *readnode = (struct ReadNode*)malloc(sizeof(struct ReadNode));
				kanntennode->ReadNode = readnode;
				readnode->dir = dir;
				readnode->pos = 0;
				readnode->flag = 0;
				readnode->read = read;
				readnode->next = redGraph->array[redNode].headread;
				redGraph->array[redNode].headread = readnode;
				link = link->next;
			}
			redGraph->array[redNode].len = 1;
			iParent = collectInternalNodes_v2(i,redNode);
			panode = graph->array[iParent].head;
			while(panode){
				redNodeParent = hashTable32_find(ht,panode->dest);
				if(!redNodeParent){
					printf("ParentNode is not part of the HashTable (But have to be impossible)\n");
					exit(1);
				}
				tempedge = (struct edge*)malloc(sizeof(struct edge));
				tempedge->dest = redNodeParent;
				tempedge->next = redGraph->array[redNode].head;
				redGraph->array[redNode].head = tempedge;
				tempedge = (struct edge*)malloc(sizeof(struct edge));
				tempedge->dest = redNode;
				tempedge->next = redGraph->array[redNodeParent].tail;
				redGraph->array[redNodeParent].tail = tempedge;
				panode = panode->next;
			}
		}
	}
}

/* 	GraphTraversal to provide absolute node positions in the graph. Each node knows its position in the graph
	So, it's getting easier to say, a node is above or under an other node !!!
 	take signed integer to label the nodes starting with 0 */

//void reklabelGraph(int i){
//	redGraph->vFlag[i]=1;
//	struct edge *neighbor;
//	int depth = redGraph->array[i].depth;
//	neighbor = redGraph->array[i].head;
//	while(neighbor){
//		if(!redGraph->vFlag[neighbor->dest]){
//			redGraph->array[neighbor->dest].depth = depth + redGraph->array[i].len;
//			reklabelGraph(neighbor->dest);
//		}
//		neighbor = neighbor->next;
//	}
//	neighbor = redGraph->array[i].tail;
//	while(neighbor){
//		if(!redGraph->vFlag[neighbor->dest]){
//			redGraph->array[neighbor->dest].depth = depth - redGraph->array[neighbor->dest].len;
//			reklabelGraph(neighbor->dest);
//		}
//		neighbor = neighbor->next;
//	}
//	return;
//}

/* 	Calls reklabelGraph for one vertex in each component
	The componed will be label rekursively after first call
	Function is neccessary at the moment for fast search of the second read end,
	but is realy problematic because it can only handle trees without loops.
	The loop ends are labeld wrong because the end has a deeper depth than the following sequence*/

void labelGraph(){
	int i;
	for(i=1;i<redGraph->V;i++){
//		if(!redGraph->vFlag[i]){
//			redGraph->array[i].depth = 0;
//			reklabelGraph(i);
//		}
		redGraph->vFlag[i]=1;
	}
//	unsigned key_count = HASH_CNT(hhb,allReads);
//	fprintf(stderr,"number of remaining Reads: %u\n", key_count);
}

/*	Prepare the non-reduced graph for the reduction by deleting multiple edges between same nodes */

void prepareGraph(){
	int i;
	struct AdjListNode *s,*p,*q;

	for(i = 1; i < graph->V; i++){
		s = graph->array[i].head;
		p = graph->array[i].head;
		if(s && s->next){
			while(p){
				s = p;
				while(s->next){
					if(p->dest == s->next->dest){
						q = s->next->next;
						free(s->next);
						s->next = q;
					}
					else s = s->next;
				}
				p = p->next;
			}
		}
		s = graph->array[i].tail;
		p = graph->array[i].tail;
		if(s && s->next){
			while(p){
				s = p;
				while(s->next){
					if(p->dest == s->next->dest){
						q = s->next->next;
						free(s->next);
						s->next = q;
					}
					else s = s->next;
				}
				p = p->next;
			}
		}
	}
}
