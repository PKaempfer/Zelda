/*
 * DBgraph.c
 *
 *  Created on: Oct 13, 2014
 *      Author: kaempfpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FileReader.h"
#include "kmer.h"
#include "DBGraph.h"
#include "kmerHash.h"

#define LINKS

//struct adTab *graph;
//struct kmerList *upOld;
//struct kmerList *downOld;
//struct kmerList *upNew;
//struct kmerList *downNew;

struct Graph *graph;
uint32_t tabindex;
int minovlLen;

char* revRead(char *read){
	int i,j=0;
	int len = strlen(read);
	char *revread = NULL;
	revread = (char*)malloc(sizeof(char)*len+1);
	for(i=len-1; i>=0; i--){
		revread[j++] = compcodes[read[i]-65];
	}
	revread[len]='\0';
	return revread;
}

void revReadSt(char *read, char* revread){
	int i,j=0;
	int len = strlen(read);
//	printf("RevRead -> ReadLen: %i\n",len);
	for(i=len-1; i>=0; i--){
		revread[j++] = compcodes[read[i]-65];
	}
	revread[len]='\0';
//	printf("RevRead -> ReadLen: %i\n",strlen(revread));
}

void writeDot(char *fileName){
	FILE* dot;
	if(fileName){
		dot = fopen(fileName,"w");
	}
	else dot = fopen("output/test.dot","w");
	int i;
	struct AdjListNode *edge;
	struct LinkListNode *link;
	fprintf(dot,"digraph dbg {\n");
//	printf("digraph dbg {\n");
	for(i=1;i<graph->V;i++){
		if(graph->array[i].head || graph->array[i].tail){
#ifdef LINKS
			if(graph->array[i].Link){
				link = graph->array[i].Link;
				fprintf(dot,"%i [label=\"%i:",i,i);
				while(link){
					fprintf(dot," ->%i",link->ID & DEL_READ_END);
					link = link->next;
				}
				fprintf(dot,"\"];\n");
			}
			else fprintf(dot,"%i [label=\"%i\"];\n",i,i);
#else
			fprintf(dot,"%i [label=\"%i\"];\n",i,i);
#endif

//			printf("%i [label=\"%i\"];\n",i,i);
			edge=graph->array[i].head;
			while(edge){
				fprintf(dot,"%i -> %i [label=%c];\n",edge->dest,i,toCharTrans(edge->trans));
//				printf("%i -> %i [label=%c];\n",edge->dest,i,edge->trans);
				edge = edge->next;
			}
			edge = graph->array[i].tail;
			while(edge){
				fprintf(dot,"%i -> %i [label=%c];\n",i,edge->dest,toCharTrans(edge->trans)+32);
				edge = edge->next;
			}
		}
	}
	fprintf(dot,"}");
	fclose(dot);
}

void printGraph()
{
    int v;
    int end;
    for (v = 0; v < graph->V; ++v)
    {
        struct AdjListNode* pCrawl = graph->array[v].head;
        if(graph->array[v].head || graph->array[v].tail){
            printf("\n Adjacency list of vertex %d\n head ", v);
            while (pCrawl)
            {
                printf("-> %d (%c)", pCrawl->dest, toCharTrans(pCrawl->trans));
                pCrawl = pCrawl->next;
            }
            printf("\n tail ");
            pCrawl = graph->array[v].tail;
            while (pCrawl){
            	printf("-> %d (%c)", pCrawl->dest, toCharTrans(pCrawl->trans));
                pCrawl = pCrawl->next;
            }
            printf("\n");
            struct LinkListNode* pCrawlLink = graph->array[v].Link;
            if(pCrawlLink){
            	printf(" Read: ");
            	while(pCrawlLink){
            		if(__builtin_clz(pCrawlLink->ID)) end = 0;
            		else end = 1;
            		printf("-> %i (%i)",(pCrawlLink->ID & DEL_READ_END),end);
            		pCrawlLink = pCrawlLink->next;
            	}
            }
            printf("\n");
        }
    }
}


void createGraph(int V)
{
    graph = (struct Graph*) malloc(sizeof(struct Graph));
    graph->V = V;

    graph->array = (struct AdjList*) malloc((V+1) * sizeof(struct AdjList));
    graph->vFlag = (char*) malloc((V+1) * sizeof(char));

    int i;
    for (i = 0; i < V; ++i){
        graph->array[i].head = NULL;
    	graph->array[i].tail = NULL;
    	graph->array[i].Link = NULL;
    	graph->vFlag[i] = 0;
    }
}

void freeGraph(){
	printf("CHECKPOINT: Free AdjList\n");
	uint32_t i;
	struct AdjListNode* listNode;
	struct LinkListNode* linkNode;
	for(i = 0; i < graph->V; i++){
		listNode = graph->array[i].head;
		while(listNode){
			graph->array[i].head = listNode->next;
			free(listNode);
			listNode = graph->array[i].head;
		}
		listNode = graph->array[i].tail;
		while(listNode){
			graph->array[i].tail = listNode->next;
			free(listNode);
			listNode = graph->array[i].tail;
		}
		linkNode = graph->array[i].Link;
		while(linkNode){
			graph->array[i].Link = linkNode->next;
			free(linkNode);
			linkNode = graph->array[i].Link;
		}
	}
	free(graph->array);
	free(graph->vFlag);
	free(graph);
}

struct AdjListNode* newAdjListNode(int dest){
    struct AdjListNode* newNode = (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

struct LinkListNode* newLinkListNode(readID ID){
    struct LinkListNode* newNode = (struct LinkListNode*) malloc(sizeof(struct LinkListNode));
    newNode->ID = ID;
    newNode->next = NULL;
    return newNode;
}

void addEdge(int src, int dest, char base, char base2){
	if(src == dest) return;
    struct AdjListNode* newNodeParent = newAdjListNode(dest);
    newNodeParent->next = graph->array[src].head;
    newNodeParent->trans = toBinTrans(base);
    graph->array[src].head = newNodeParent;

    struct AdjListNode* newNodeChild = newAdjListNode(src);
    newNodeChild->next = graph->array[dest].tail;
    newNodeChild->trans = toBinTrans(base2);
    graph->array[dest].tail = newNodeChild;
}

void printSorrounding(int ori){
    struct AdjListNode* pCrawl = graph->array[ori].head;
    int end;
    printf("TransitionTable of Node %i: \n", ori);
    printf(" head:");
    while(pCrawl){
    	printf(" -> %i (%c)",pCrawl->dest,toCharTrans(pCrawl->trans));
    	pCrawl = pCrawl->next;
    }
    pCrawl = graph->array[ori].tail;
    printf("\n tail:");
    while(pCrawl){
    	printf(" -> %i (%c)",pCrawl->dest,toCharTrans(pCrawl->trans));
    	pCrawl = pCrawl->next;
    }
    printf("\n");
    struct LinkListNode* pCrawlLink = graph->array[ori].Link;
    if(pCrawlLink){
    	printf(" Read: ");
    	while(pCrawlLink){
    		if(__builtin_clz(pCrawlLink->ID)) end = 0;
    		else end = 1;
    		printf("-> %i (%i)",(pCrawlLink->ID & DEL_READ_END),end);
    		pCrawlLink = pCrawlLink->next;
    	}
    }
    printf("\n");

}

void countKmers(){
	int i;
	int count = 0;
	for(i=1;i<graph->V;i++){
		if(graph->array[i].head || graph->array[i].tail) count++;
	}
	printf("Number of left Kmers: %i\n",count);
}

void deleteComp(int start){
	struct AdjListNode *Node;
	int temp = start;
	for(;start<tabindex;start++){
		Node=graph->array[start].head;
		while(Node){
			graph->array[start].head=Node->next;
			free(Node);
			Node=graph->array[start].head;
		}
		Node=graph->array[start].tail;
		while(Node){
			graph->array[start].tail=Node->next;
			free(Node);
			Node=graph->array[start].tail;
		}
	}
	tabindex=temp;
}


//void hashToTabDFS(){
//    tabindex=1;
//    int comp = 0;
//    int oldIndex = 1;
//    int delcomp = 0, delNode = 0;
//
//#ifdef UT_HASH
//    struct hashTable *s, *t;
//    s = (struct hashTable*)malloc(sizeof(struct hashTable));
//    HASH_ITER(hhb, kmere, s, t){
//    	if(!s->trans){
//    		if(s->index==0){
//    			s->index=tabindex;
//    			tabindex++;
//    		}
//    		travDFS(s);
//
//    		if(tabindex-oldIndex<MINCOMP){
//    			graph->V-=tabindex-oldIndex;
//    			delNode += tabindex-oldIndex;
//    			deleteComp(oldIndex);
//    			delcomp++;
//    		}
//    		else{
//    			comp ++;
//    			printf("%i. Component, size: %i\n",comp,tabindex-oldIndex);
//    		}
//    		oldIndex =  tabindex;
//
//    	}
//    }
//    travDFS(NULL);
//#else
//
//    if(sizeof(KmerBitBuffer)*8==64){
////    	int j = INITHASHSIZE(bitnum);
////    	int i;
////    	for(i=0;i<j;i++){
////    		if(dbHash_oa[i].kmer != empty && !dbHash_oa[i].trans){
////        		if(dbHash_oa[i].index==0){
////        			dbHash_oa[i].index = tabindex;
////        			tabindex++;
////
////        		}
////        		travDFS(&dbHash_oa[i]);
////        		if(tabindex-oldIndex < MINCOMP){
////        			graph->V -= tabindex - oldIndex;
////        			delNode += tabindex - oldIndex;
////        			deleteComp(oldIndex);
////        			delcomp++;
////        		}
////        		else{
////        			comp++;
////        			printf("%i. Component, size: %i\n",comp,tabindex-oldIndex);
////        		}
////        		oldIndex = tabindex;
////    		}
////
////    	}
////    	travDFS(NULL);
//    }
//    else{
//        // My Hash Lib
//        struct hashkmer *s;
//        setIter();
//        while((s = iterKmer())){
//        	if(!s->trans){
//        		if(s->index==0){
//        			s->index=tabindex;
//        			tabindex++;
//        		}
//        		travDFS((struct hashTable*)s);
//
//        		if(tabindex-oldIndex<MINCOMP){
//        			graph->V-=tabindex-oldIndex;
//        			delNode += tabindex-oldIndex;
//        			deleteComp(oldIndex);
//        			delcomp++;
//        		}
//        		else{
//        			comp ++;
//        			printf("%i. Component, size: %i\n",comp,tabindex-oldIndex);
//        		}
//        		oldIndex =  tabindex;
//
//        	}
//        }
//        travDFS(NULL);
//    }
//
//#endif
//
//    printf("%i Nodes in %i Components\n",tabindex-1,comp);
//    printf("%i Nodes in %i Components deleted\n",delNode,delcomp);
//}

void addLinks(readID *ID, int index){
	int i=0;
	do{
		struct LinkListNode* newNodeParent = newLinkListNode(ID[i]);
		newNodeParent->next = graph->array[index].Link;
		graph->array[index].Link = newNodeParent;
		i++;
	} while(ID[i]);
}

//void goUpDFS(struct hashTable ***upStackold,struct hashTable ***downStackold, int *upPtr, int *downPtr, int *upSPtr, int *downSPtr, int upbool){
//	int up = (*upPtr);
//	int down = (*downPtr);
//	int upSize = (*upSPtr);
//	int downSize = (*downSPtr);
//	struct hashTable *s, *p;
//	struct hashTable **tmp,**upStack,**downStack;
//	downStack = (*downStackold);
//	upStack = (*upStackold);
//
//	static int edgeNum=0;
//	char b;
//	char transbase;
//	char transbaseDown;
//	char dir,pdir;
//	int updown;
//
//	struct readLink *reads;
//
//	do{
//		//if(!upbool)
//		if(upbool)	s = upStack[--up];
//		if(!upbool) s = downStack[--down];
//
//		if(!s->trans){
//
//			// Put read links to adjacency list
//			HASH_FIND(hhb, links, &s->kmer, sizeof(KmerBitBuffer),reads);
//			if(reads) addLinks(reads->read,ABS(s->index));
//
//			transbase = getTransBase(&s->kmer,s->index);
//
//			if(s->index>0) dir = 1;
//			else dir=-1;
//
//			for(b=0;b<4;b++){
//				pdir=dir;
//				p = hasParent(&s->kmer,b,&pdir);
//				if(p){
//					if(p->index == 0){
//						p->index = pdir * tabindex;
//						tabindex++;
//					}
//					if(dir == 1){
//						transbaseDown = getTransBaseDown(&p->kmer,p->index);
//						addEdge(ABS((s->index)),ABS((p->index)),transbase,transbaseDown);
//						edgeNum++;
//						if(!p->trans) upStack[up++] = p;
//						if(up==upSize-1){
//							tmp = (struct hashTable**)realloc(upStack,sizeof(struct hashTable*)*(upSize*2));
//							if(tmp){
//								upStack = tmp;
//								upSize *= 2;
//
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//
//						}
//					}
//					else{
//						if(!p->trans) downStack[down++] = p;
//						if(down==downSize-1){
//							tmp = (struct hashTable**)realloc(downStack,sizeof(struct hashTable*)*(downSize*2));
//							if(tmp){
//								downStack = tmp;
//								downSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//				}
//
//				pdir=dir;
//				p = hasChild(&s->kmer,b,&pdir);
//				if(p){
//					if(p->index == 0){
//						p->index = pdir * tabindex;
//						tabindex++;
//					}
//					if(dir == -1){
//						transbaseDown = getTransBaseDown(&p->kmer,p->index);
//						addEdge(ABS((s->index)),ABS((p->index)),transbase,transbaseDown);
//						edgeNum++;
//						if(!p->trans) upStack[up++] = p;
//						if(up==upSize-1){
//							tmp = (struct hashTable**)realloc(upStack,sizeof(struct hashTable*)*(upSize*2));
//							if(tmp){
//								printf("Realloc: New UpSize: %i (%p %p)\n",upSize,tmp,upStack);
//								upStack = tmp;
//								upSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//					else{
//						if(!p->trans) downStack[down++] = p;
//						if(down==downSize-1){
//							tmp = (struct hashTable**)realloc(downStack,sizeof(struct hashTable*)*(downSize*2));
//							if(tmp){
//								printf("Realloc: New DownSize: %i (%p %p)\n",downSize,tmp,downStack);
//								downStack = tmp;
//								downSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//				}
//
//				s->trans = 1;
//			}
//		}
//		if(upbool) updown = up;
//		else updown = down;
//	} while(updown);
//
////	printf("up: %i Upsize: %i Down: %i DownSize: %i\n", up, upSize, down, downSize);
//
//	(*upPtr) = up;
//	(*downPtr) = down;
//	(*upSPtr) = upSize;
//	(*downSPtr) = downSize;
//	(*downStackold) = downStack;
//	(*upStackold) = upStack;
//}
//
//void goUpDFS_2(struct hashTable ***upStackold,struct hashTable ***downStackold, int *upPtr, int *downPtr, int *upSPtr, int *downSPtr, int upbool){
//	int up = (*upPtr);
//	int down = (*downPtr);
//	int upSize = (*upSPtr);
//	int downSize = (*downSPtr);
//	struct hashkmer *s, *p;
//	struct hashkmer **tmp,**upStack,**downStack;
//	downStack = (struct hashkmer**)(*downStackold);
//	upStack = (struct hashkmer**)(*upStackold);
//
//	static int edgeNum=0;
//	char b;
//	char transbase;
//	char transbaseDown;
//	char dir,pdir;
//	int updown;
//
//	struct readLink *reads;
//
//	do{
//		//if(!upbool)
//		if(upbool)	s = upStack[--up];
//		if(!upbool) s = downStack[--down];
//
//		if(!s->trans){
//
//			// Put read links to adjacency list
//			HASH_FIND(hhb, links, &s->kmer, sizeof(KmerBitBuffer),reads);
//			if(reads) addLinks(reads->read,ABS(s->index));
//
//			transbase = getTransBase(&s->kmer,s->index);
//
//			if(s->index>0) dir = 1;
//			else dir=-1;
//
//			for(b=0;b<4;b++){
//				pdir=dir;
//				p = hasParent_2(&s->kmer,b,&pdir);
//				if(p){
//					if(p->index == 0){
//						p->index = pdir * tabindex;
//						tabindex++;
//					}
//					if(dir == 1){
//						transbaseDown = getTransBaseDown(&p->kmer,p->index);
//						addEdge(ABS((s->index)),ABS((p->index)),transbase,transbaseDown);
//						edgeNum++;
//						if(!p->trans) upStack[up++] = p;
//						if(up==upSize-1){
//							tmp = (struct hashkmer**)realloc(upStack,sizeof(struct hashkmer*)*(upSize*2));
//							if(tmp){
//								upStack = tmp;
//								upSize *= 2;
//
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//
//						}
//					}
//					else{
//						if(!p->trans) downStack[down++] = p;
//						if(down==downSize-1){
//							tmp = (struct hashkmer**)realloc(downStack,sizeof(struct hashkmer*)*(downSize*2));
//							if(tmp){
//								downStack = tmp;
//								downSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//				}
//
//				pdir=dir;
//				p = hasChild_2(&s->kmer,b,&pdir);
//				if(p){
//					if(p->index == 0){
//						p->index = pdir * tabindex;
//						tabindex++;
//					}
//					if(dir == -1){
//						transbaseDown = getTransBaseDown(&p->kmer,p->index);
//						addEdge(ABS((s->index)),ABS((p->index)),transbase,transbaseDown);
//						edgeNum++;
//						if(!p->trans) upStack[up++] = p;
//						if(up==upSize-1){
//							tmp = (struct hashkmer**)realloc(upStack,sizeof(struct hashkmer*)*(upSize*2));
//							if(tmp){
//								printf("Realloc: New UpSize: %i (%p %p)\n",upSize,tmp,upStack);
//								upStack = tmp;
//								upSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//					else{
//						if(!p->trans) downStack[down++] = p;
//						if(down==downSize-1){
//							tmp = (struct hashkmer**)realloc(downStack,sizeof(struct hashkmer*)*(downSize*2));
//							if(tmp){
//								printf("Realloc: New DownSize: %i (	 %p)\n",downSize,tmp);
//								downStack = tmp;
//								downSize *=2;
//							}
//							else printf("Error: could not reallocate Stack memory\n");
//						}
//					}
//				}
//
//				s->trans = 1;
//			}
//		}
//		if(upbool) updown = up;
//		else updown = down;
//	} while(updown);
//
////	printf("up: %i Upsize: %i Down: %i DownSize: %i\n", up, upSize, down, downSize);
//
//	(*upPtr) = up;
//	(*downPtr) = down;
//	(*upSPtr) = upSize;
//	(*downSPtr) = downSize;
//	(*downStackold) = (struct hashTable**) downStack;
//	(*upStackold) = (struct hashTable**) upStack;
//}

//void travDFS(struct hashTable* s){
//	static int upStacksize=1024;
//	static int downStacksize=1024;
//	int down = 0;
//	int up = 0;
//	static struct hashTable **downStack = NULL;
//	static struct hashTable **upStack = NULL;
//	if(!downStack) downStack = (struct hashTable**)malloc(sizeof(struct hashTable*)*1024);
//	if(!upStack) upStack = (struct hashTable**)malloc(sizeof(struct hashTable*)*1024);
//
//	if(s){
//		upStack[up++]=s;
//		do{
//#ifdef UT_HASH
//			goUpDFS(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,1);
//			if(down) goUpDFS(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,0);
//#else
//			goUpDFS_2(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,1);
//			if(down) goUpDFS_2(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,0);
//#endif
//		} while(up);
//	}
//	else{
//		free(upStack);
//		free(downStack);
//	}
//}

char isChild(int ori, int chi){
	struct AdjListNode *child = graph->array[ori].tail;
	while(child){
		if(child->dest == chi) return 1;
		child = child->next;
	}
	return 0;
}

char isParent(int ori, int par){
	struct AdjListNode *parent = graph->array[ori].head;
	while(parent){
		if(parent->dest == par) return 1;
		parent = parent->next;
	}
	return 0;
}