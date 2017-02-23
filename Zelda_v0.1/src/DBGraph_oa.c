/*
 ============================================================================
 Name        : DBGraph_oa.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Converts the dBG from (random access) hash table to cache
 	 	 	   coherent adjacency list and calculates edges
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include "DBGraph.h"
#include "DBGraph_oa.h"
#include "kmerHash.h"
#include "kmerHash_oa.h"
#include "kmer.h"

void addLinks_oa(struct readEnd* readend, int index){
	graph->array[index].Link = (struct LinkListNode*) readend;

//	while(readend){
//		struct LinkListNode* newNodeParent = newLinkListNode(readend->read);
//		newNodeParent->next = graph->array[index].Link;
//		graph->array[index].Link = newNodeParent;
//		readend = (struct readEnd*)readend->next;
//	}
}

uint32_t getKmer_oa(KmerBitBuffer kmer){
	static uint32_t bucket;
	bucket = my_hash(kmer);
	int i = 0;

	while(dbHash_oa[bucket].count && i < max_reprobes){
		if(dbHash_oa[bucket].kmer == kmer) return bucket;
		bucket++;
		i++;
	}
	return EMPTYBUCKET;
}

void hashToTabDFS_oa(){
    tabindex=1;
    int comp = 0;
    int oldIndex = 1;
    int delcomp = 0, delNode = 0, delEnd = 0;

	uint32_t j = INITHASHSIZE(bitnum);
	uint32_t i;
	uint32_t k;
	struct LinkListNode* readEnd;
	for(i=0;i<j;i++){
		if(dbHash_oa[i].count && !dbHash_oa[i].trans){
    		if(dbHash_oa[i].index==0){
    			dbHash_oa[i].index = tabindex;
    			tabindex++;

    		}
    		travDFS_oa(i);
    		if(tabindex-oldIndex < MINCOMP){
    			for(k=oldIndex;k<tabindex;k++){
    				readEnd = graph->array[k].Link;
    				while(readEnd){
    					delEnd++;
    					graph->array[k].Link = readEnd->next;
    					free((void*)readEnd);
    					readEnd = graph->array[k].Link;
    				}
    			}
//    			printf("%i Ends deleted\n",delEnd);
    			graph->V -= tabindex - oldIndex;
    			delNode += tabindex - oldIndex;
    			deleteComp(oldIndex);
    			delcomp++;
    		}
    		else{
    			comp++;
    		}
    		oldIndex = tabindex;
		}
	}
	travDFS_oa(EMPTYBUCKET);

    printf("%i Nodes in %i Components\n",tabindex-1,comp);
    printf("%i Ends deleted\n",delEnd);
    printf("%i Nodes in %i Components deleted\n",delNode,delcomp);
}

void travDFS_oa(uint32_t i){
	static uint32_t upStacksize=1024*256;
	static uint32_t downStacksize=1024*256;
	uint32_t down = 0;
	uint32_t up = 0;
	static uint32_t* downStack = NULL;
	static uint32_t* upStack = NULL;
	if(!downStack) downStack = (uint32_t*)malloc(sizeof(uint32_t)*downStacksize);
	if(!upStack) upStack = (uint32_t*)malloc(sizeof(uint32_t)*upStacksize);

	if(i == EMPTYBUCKET){
		printf("Free all intermediate stacks\n");
		free(upStack);
		free(downStack);
	}
	else if(dbHash_oa[i].count){
		upStack[up++]=i;
		do{
			goUpDFS_2_oa(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,1);
			if(down) goUpDFS_2_oa(&upStack, &downStack, &up, &down, &upStacksize, &downStacksize,0);
		} while(up);
	}

}


void goUpDFS_2_oa(uint32_t** upStackold,uint32_t** downStackold, uint32_t* upPtr, uint32_t* downPtr, uint32_t* upSPtr, uint32_t* downSPtr, int upbool){
	char verbose = 0;
	uint32_t up = (*upPtr);
	uint32_t down = (*downPtr);
	uint32_t upSize = (*upSPtr);
	uint32_t downSize = (*downSPtr);
	uint32_t s, p;
	uint32_t *tmp,*upStack,*downStack;
	downStack = (uint32_t*)(*downStackold);
	upStack = (uint32_t*)(*upStackold);

	static int edgeNum=0;
	char b;
	char transbase;
	char transbaseDown;
	char dir,pdir;
	int updown;
	KmerBitBuffer p_kmer, s_kmer;

	struct readEnd* reads;

	do{
		if(upbool)	s = upStack[--up];
		if(!upbool) s = downStack[--down];


		if(!dbHash_oa[s].trans){

			// Put read links to adjacency list
			reads = (struct readEnd*)dbHash_oa[s].ends;
			addLinks_oa(reads,ABS(dbHash_oa[s].index));
			graph->array[ABS(dbHash_oa[s].index)].counter = dbHash_oa[s].count;

			s_kmer = (KmerBitBuffer)dbHash_oa[s].kmer;
			transbase = getTransBase(&s_kmer,dbHash_oa[s].index);

			if(dbHash_oa[s].index>0) dir = 1;
			else dir=-1;

			for(b=0;b<4;b++){
				pdir=dir;
				s_kmer =  (KmerBitBuffer)dbHash_oa[s].kmer;
				p = hasParent_2_oa(&s_kmer,b,&pdir);
				if(p!=EMPTYBUCKET){
					p_kmer = (KmerBitBuffer)dbHash_oa[p].kmer;
					if(dbHash_oa[p].index == 0){
						dbHash_oa[p].index = pdir * tabindex;
						tabindex++;
					}
					if(dir == 1){
						p_kmer = (KmerBitBuffer)dbHash_oa[p].kmer;
						transbaseDown = getTransBaseDown(&p_kmer,dbHash_oa[p].index);
						addEdge(ABS((dbHash_oa[s].index)),ABS((dbHash_oa[p].index)),transbase,transbaseDown);
						edgeNum++;
						if(!dbHash_oa[p].trans){
							upStack[up++] = p;
						}
						// Resize Up-Stack
						if(up==upSize-1){
							tmp = (uint32_t*)realloc(upStack,sizeof(uint32_t)*(upSize*2));
							if(tmp){
								upStack = tmp;
								upSize *= 2;

							}
							else{
								printf("Error: could not reallocate Stack memory\n");
								exit(1);
							}
						}
					}
					else{
						if(!dbHash_oa[p].trans) downStack[down++] = p;
						// Resize Down-Stack
						if(down==downSize-1){
							tmp = (uint32_t*)realloc(downStack,sizeof(uint32_t)*(downSize*2));
							if(tmp){
								downStack = tmp;
								downSize *=2;
							}
							else printf("Error: could not reallocate Stack memory\n");
						}
					}
				}

				pdir=dir;
				s_kmer = (KmerBitBuffer)dbHash_oa[s].kmer;
				p = hasChild_2_oa(&s_kmer,b,&pdir);
				if(p!=EMPTYBUCKET){
//					printf("child found!\n");
					if(dbHash_oa[p].index == 0){
						dbHash_oa[p].index = pdir * tabindex;
						tabindex++;
					}
					if(dir == -1){
						p_kmer = (KmerBitBuffer)dbHash_oa[p].kmer;
						transbaseDown = getTransBaseDown(&p_kmer,dbHash_oa[p].index);
						addEdge(ABS((dbHash_oa[s].index)),ABS((dbHash_oa[p].index)),transbase,transbaseDown);
						edgeNum++;
						if(!dbHash_oa[p].trans) upStack[up++] = p;
						// Resize Up-Stack
						if(up==upSize-1){
							tmp = (uint32_t*)realloc(upStack,sizeof(uint32_t)*(upSize*2));
							if(tmp){
								if(verbose) printf("Realloc: New UpSize: %i (%p %p)\n",upSize,tmp,upStack);
								upStack = tmp;
								upSize *=2;
							}
							else printf("Error: could not reallocate Stack memory\n");
						}
					}
					else{
						if(!dbHash_oa[p].trans) downStack[down++] = p;
						// Resize Down-Stack
						if(down==downSize-1){
							tmp = (uint32_t*)realloc(downStack,sizeof(uint32_t)*(downSize*2));
							if(tmp){
								if(verbose) printf("Realloc: New DownSize: %i (	%p %p)\n",downSize,tmp,downStack);
								downStack = tmp;
								downSize *=2;
							}
							else printf("Error: could not reallocate Stack memory\n");
						}
					}
				}

				dbHash_oa[s].trans = 1;
			}
		}
		if(upbool) updown = up;
		else updown = down;
	} while(updown);

//	printf("up: %i Upsize: %i Down: %i DownSize: %i\n", up, upSize, down, downSize);

	(*upPtr) = up;
	(*downPtr) = down;
	(*upSPtr) = upSize;
	(*downSPtr) = downSize;
	(*downStackold) = (uint32_t*) downStack;
	(*upStackold) = (uint32_t*) upStack;
}
