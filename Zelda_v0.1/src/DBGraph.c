/*
 ============================================================================
 Name        : DBGraph.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Meta-functions for dBG graph handling
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "FileReader.h"
#include "kmer.h"
#include "DBGraph.h"
#include "kmerHash.h"

#define LINKS

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
	for(i=len-1; i>=0; i--){
		revread[j++] = compcodes[read[i]-65];
	}
	revread[len]='\0';
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

			edge=graph->array[i].head;
			while(edge){
				fprintf(dot,"%i -> %i [label=%c];\n",edge->dest,i,toCharTrans(edge->trans));
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
    for (i = 0; i <= V; ++i){
        graph->array[i].head = NULL;
    	graph->array[i].tail = NULL;
    	graph->array[i].Link = NULL;
    	graph->vFlag[i] = 0;
    }
}

void freeGraph(){
//	printf("CHECKPOINT: Free AdjList\n");
	uint32_t i;
	struct AdjListNode* listNode;
	struct LinkListNode* linkNode;
	for(i = 0; i <= graph->V; i++){
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
	printf("\tNumber of left Kmers: %i\n",count);
}

void deleteComp(int start){
	struct AdjListNode *node;
	struct LinkListNode *link;
	int temp = start;
	for(;start<tabindex;start++){
		node=graph->array[start].head;
		while(node){
			graph->array[start].head=node->next;
			free(node);
			node=graph->array[start].head;
		}
		node=graph->array[start].tail;
		while(node){
			graph->array[start].tail=node->next;
			free(node);
			node=graph->array[start].tail;
		}
		link = graph->array[start].Link;
		while(link){
			graph->array[start].Link = link->next;
			free(link);
			link = graph->array[start].Link;
		}
	}
	tabindex=temp;
}


void addLinks(readID *ID, int index){
	int i=0;
	do{
		struct LinkListNode* newNodeParent = newLinkListNode(ID[i]);
		newNodeParent->next = graph->array[index].Link;
		graph->array[index].Link = newNodeParent;
		i++;
	} while(ID[i]);
}

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
