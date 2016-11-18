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
	while(readend){
		struct LinkListNode* newNodeParent = newLinkListNode(readend->read);
		newNodeParent->next = graph->array[index].Link;
		graph->array[index].Link = newNodeParent;
		readend = (struct readEnd*)readend->next;
	}
}

uint32_t getKmer_oa(KmerBitBuffer kmer){
	static uint32_t bucket;
	bucket = my_hash(kmer);
//	printf("Search at Bucket. %i\n",bucket);
	int i = 0;

	while(dbHash_oa[bucket].count && i < max_reprobes){
//		printf("Search for The Seq\n");
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
    int delcomp = 0, delNode = 0;

	uint32_t j = INITHASHSIZE(bitnum);
	uint32_t i;
	for(i=0;i<j;i++){
		if(dbHash_oa[i].count && !dbHash_oa[i].trans){
    		if(dbHash_oa[i].index==0){
    			dbHash_oa[i].index = tabindex;
    			tabindex++;

    		}
//    		printf("Start Component building at position %i\n",i);
    		travDFS_oa(i);
//    		printf("Finished component\n");
    		if(tabindex-oldIndex < MINCOMP){
//    			printf("Component to small: %i -> Delete Component\n",tabindex-oldIndex);
//    			printf("size of Graph: %i\n",graph->V);
    			graph->V -= tabindex - oldIndex;
    			delNode += tabindex - oldIndex;
    			deleteComp(oldIndex);
    			delcomp++;
//    			printf("Deletion Complete\n");
    		}
    		else{
    			comp++;
//    			printf("%i. Component, size: %i\n",comp,tabindex-oldIndex);
    		}
    		oldIndex = tabindex;
		}
	}
	travDFS_oa(EMPTYBUCKET);

    printf("%i Nodes in %i Components\n",tabindex-1,comp);
    printf("%i Nodes in %i Components deleted\n",delNode,delcomp);
}

void travDFS_oa(uint32_t i){
	static uint32_t upStacksize=1024;
	static uint32_t downStacksize=1024;
	uint32_t down = 0;
	uint32_t up = 0;
	static uint32_t* downStack = NULL;
	static uint32_t* upStack = NULL;
	if(!downStack) downStack = (uint32_t*)malloc(sizeof(uint32_t)*1024);
	if(!upStack) upStack = (uint32_t*)malloc(sizeof(uint32_t)*1024);

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
//		printf("New Run (Up: %i / Down: %i)\n",up,down);
		if(upbool)	s = upStack[--up];
		if(!upbool) s = downStack[--down];

		if(!dbHash_oa[s].trans){

			// Put read links to adjacency list
			reads = (struct readEnd*)dbHash_oa[s].ends;
			addLinks_oa(reads,ABS(dbHash_oa[s].index));

			s_kmer = (KmerBitBuffer)dbHash_oa[s].kmer;
			transbase = getTransBase(&s_kmer,dbHash_oa[s].index);

			if(dbHash_oa[s].index>0) dir = 1;
			else dir=-1;

			for(b=0;b<4;b++){
				pdir=dir;
				s_kmer =  (KmerBitBuffer)dbHash_oa[s].kmer;
				p = hasParent_2_oa(&s_kmer,b,&pdir);
				if(p!=EMPTYBUCKET){
//					printf("current Kmer: %s\n",toSeq(s_kmer));
					p_kmer = (KmerBitBuffer)dbHash_oa[p].kmer;
//					printf("Parent found: %s (edge transition: %c)\n",toSeq(p_kmer),rev_codes[(int)b]);
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
//							printf("pop up parent\n");
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
								printf("Realloc: New UpSize: %i (%p %p)\n",upSize,tmp,upStack);
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
								printf("Realloc: New DownSize: %i (	%p %p)\n",downSize,tmp,downStack);
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
