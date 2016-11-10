/*
 * ConsensusCaller.c
 *
 *  Created on: Feb 26, 2016
 *      Author: kaempfpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "ConsensusCaller.h"
#include "DBGraph_scaffold.h"
//#include "DBGraph_stringer.h"
#include "readDB.h"
#include "DBGraph.h"


#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define MIN_CONTIG_LEN 100
#define MIN_SCAFF_LEN 100
#define MATRIX_MAX_BR 20

static char status_char[] = { 'W', 'C', 'S', 'P', 'J' };
struct timespec ts_start;
struct timespec ts_finish;
long sumMatrix = 0;
long sumTrace = 0;

void buildBackBone(struct myovlList* ovlgraph, struct string_graph* S, struct reads* reads){
	printf("MemoPos of OvlGraph: %p\n",ovlgraph);
    int       nends  = S->nends;
    unsigned char*     status = S->status;
    struct readend*  side   = S->side;
    struct overhang* edge   = S->edge;

    int       v, w, e, o;

    int i;
    int run = 1;
//    int breadID;
    struct bread* bread;
    char* adecomp;
    char* bdecomp;

    for (v = 1; v <= nends; v += 2)
    {
        if ( status[VERT(v)] == WIDOWED ||
             status[VERT(v)] == CONTAINED )
        {
            continue;
        }

        if(status[VERT(v)]==JUNCTION)
        printf("%5d(%i): %c\n", VERT(v),S->ID[VERT(v)], status_char[status[VERT(v)]]);

        for (o = 0; o < 2; o++){
        	w = v + o;

            for (e = side[w - 1].top + 1; e <= side[w].top; e++){
//            	if(status[VERT(v)]==JUNCTION) printf("     %c -> %i[%3d]\n", (o ? 'e' : 'b'), edge[e].target, edge[e].length);

            	if(edge[e].length >= MIN_CONTIG_LEN){
            		printf("Build new Contig\n");
                    printf("     %c -> %s [len: %3d]\n",
                           (o ? 'e' : 'b'), vtx_name(edge[e].target, 5), edge[e].length);

            		adecomp = decompressRead(reads[VERT(v)].seq,reads[VERT(v)].len);
            		printf("aread (%i): %s\n",VERT(v),adecomp);
            		free(adecomp);

            		bread = ovlgraph->read[VERT(v)]->first;
            		printf("bread: %i\n",bread->ID);
            		while(run && bread){
            			if(bread->flag == PROPER){
            				bdecomp = decompressRead(reads[bread->ID].seq,reads[bread->ID].len);
            				if(ovlgraph->read[VERT(v)]->dir != ovlgraph->read[bread->ID]->dir){
            					bdecomp = revRead(bdecomp);
            				}
            				printf("aread (%i): ",VERT(v));
            				for(i=0;i<bread->overhang;i++) printf(" ");
            				printf("%s\n",bdecomp);
            				free(bdecomp);

            				run = 0;
            			}
            			bread = bread->next;
            		}
//            		ovlgraph->read[VERT(v)]->first;
            	}
            	else{
            		printf("     <100: %c -> %s [len: %3d]\n", (o ? 'e' : 'b'), vtx_name(edge[e].target, 5), edge[e].length);
            	}

            }
        }
    }
}

// Build backbone on the basis of the overlap graph (uncollapsed string graph)
void buildBackBone3(struct myovlList* G, struct reads* reads){
    int i,j;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    char* adecomp;
    char* bdecomp;
    int dir;
    int bdir;
    int overhang;
    int nextoverhang = 0;
    int multidir;

    for(i=1; i < G->V; i++){
    	printf("i %i\n",i);
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		printf("JUNCTION found\n");
			dir=G->read[i]->dir;
    		bread = G->read[i]->first;
    		while(bread){
				breadID = bread->ID;
    			if(bread && G->read[breadID]->flag != CONTAINED && bread->dest && bread->dest->len >= MIN_CONTIG_LEN){
    				if(dir){
    					if(bread->sideflag) 	multidir = 3;
    					else					multidir = 0;
    				}
    				else{
    					if(bread->sideflag)		multidir = 1;
    					else					multidir = 2;
    				}
    				// IF (multidir % 2 == 1) THAN turn b-reads (dir==1) ELSE turn b-reads (dir==0)
    				// IF (multidir > 1) THAN turn a-read
    				bdir = bread->sideflag;
   					printf("Contig Path found: %i -> %i (len: %i) Dir: %i Side: %i\n",i,bread->dest->ID,bread->dest->len,dir,bread->sideflag);
       				adecomp = decompressRead(reads[i].seq,reads[i].len);
       				if(multidir>1) adecomp = revRead(adecomp);
       				printf("a %i (%i):\t%s\n",G->read[i]->dir,i,adecomp);
        			free(adecomp);
        			printf("b %i (%i):\t",G->read[breadID]->dir,breadID);
        			overhang = bread->overhang;
        			for(j=0;j<overhang;j++){
        				printf(" ");
        			}
        			bdecomp = decompressRead(reads[breadID].seq,reads[breadID].len);
        			if(multidir%2==1 && G->read[breadID]->dir) bdecomp = revRead(bdecomp);
        			if(multidir%2==0 && !G->read[breadID]->dir) bdecomp = revRead(bdecomp);
        			printf("%s\n",bdecomp);
        			free(bdecomp);
        			internb = bread;
    				while(G->read[breadID]->flag != JUNCTION){
    					internb = G->read[breadID]->first;
    					while(internb){
    						if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
    			    			printf("b %i (%i):\t",G->read[internb->ID]->dir,internb->ID);
    			    			nextoverhang = internb->overhang;
    			    			overhang += internb->overhang;
    			    			for(j=0;j<overhang;j++){
    			    				printf(" ");
    			    			}
    			    			bdecomp = decompressRead(reads[internb->ID].seq,reads[internb->ID].len);
    		        			if(multidir%2==1 && G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
    		        			if(multidir%2==0 && !G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
    			    			printf("%s\n",bdecomp);
    			    			free(bdecomp);
    							breadID = internb->ID;
//    							break;
    						}
    						if(G->read[internb->ID]->flag == CONTAINED){
    							printf("c %i (%i):\t",G->read[internb->ID]->dir,internb->ID);
    			    			for(j=0;j<overhang-nextoverhang;j++){
    			    				printf(" ");
    			    			}
    			    			bdecomp = decompressRead(reads[internb->ID].seq,reads[internb->ID].len);
    		        			if(multidir%2==0 && G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
    		        			if(multidir%2==1 && !G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
    			    			printf("%s\n",bdecomp);
    			    			free(bdecomp);
    						}
    						internb = internb->next;
    					}
    				}
    			}
    			bread = bread->next;
    		}
    	}
    }
}

struct contigList* realBackbone(struct myovlList* G, struct reads* reads){
	// Take junction read -> transform to POG nodes and edges
	// concatenate the overhangs of the proper flagged reads to the POG
    int i,j;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    char* backbone;
    char* adecomp;
    char* bdecomp;
//    char* readseq = (char*)malloc(maxReadLen+1);

    int dir;
    int bdir;
    int overhang;
    int nextoverhang;
    int multidir;

    struct contigList* cList = (struct contigList*)malloc(sizeof(struct contigList));
    cList->num = 0;
    cList->maxnum = 100;
    cList->contig = (struct contig*)malloc(sizeof(struct contig)*cList->maxnum);

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		printf("JUNCTION found\n");
			dir=G->read[i]->dir;
    		bread = G->read[i]->first;
    		while(bread){
				breadID = bread->ID;
    			if(bread && G->read[breadID]->flag != CONTAINED && bread->dest && bread->dest->len >= MIN_CONTIG_LEN){
//   					printf("Contig Path found: %i -> %i (len: %i) Dir: %i Side: %i\n",i,bread->dest->ID,bread->dest->len,dir,bread->sideflag);
//    				// ContigList is full, resize List
    				if(cList->num+1 == cList->maxnum){
    					printf("ContigList full, Realloc memory\n");
    					cList->maxnum *= 2;
    					cList->contig = (struct contig*)realloc(cList->contig,sizeof(struct contig)*cList->maxnum);
    				}
    				printf("Test\n");
    				cList->contig[cList->num].stReadID = i;
    				cList->contig[cList->num].len = bread->dest->len+reads[i].len;
    				printf("Len: %i\n",cList->contig[cList->num].len);
    				backbone = (char*)malloc(cList->contig[cList->num].len+1000);
    				if(!backbone){
    					printf("Could not allocate memory for backbone sequence\n");
    					exit(1);
    				}
//    				backbone[0] = '\0';
    				printf("Dirselection\n");
    				if(dir){
    					if(bread->sideflag) 	multidir = 3;
    					else					multidir = 0;
    				}
    				else{
    					if(bread->sideflag)		multidir = 1;
    					else					multidir = 2;
    				}
    				// IF (multidir % 2 == 1) THAN turn b-reads (dir==1) ELSE turn b-reads (dir==0)
    				// IF (multidir > 1) THAN turn a-read
    				bdir = bread->sideflag;
       				adecomp = decompressRead(reads[i].seq,reads[i].len);
//    				decompressRead(reads[i].seq,readseq,reads[i].len);
       				if(multidir>1) adecomp = revRead(adecomp);
     				printf("a %i (%i):\t%s\n",G->read[i]->dir,i,adecomp);
//       				strcpy(cList->contig[cList->num].seq,adecomp);
        			free(adecomp);
        			printf("b %i (%i):\t",G->read[breadID]->dir,breadID);
        			overhang = bread->overhang;
        			for(j=0;j<overhang;j++){
        				printf(" ");
        			}
        			bdecomp = decompressRead(reads[breadID].seq,reads[breadID].len);
        			if(multidir%2==1 && G->read[breadID]->dir) bdecomp = revRead(bdecomp);
        			if(multidir%2==0 && !G->read[breadID]->dir) bdecomp = revRead(bdecomp);
        			printf("%s\n",bdecomp);
        			printf("Cat: %i -> %s\n",reads[breadID].len - bread->overhang, &bdecomp[(reads[breadID].len - bread->overhang)]);
        			strcat(backbone,&bdecomp[(reads[breadID].len - bread->overhang)]);
//        			printf("Seq: %i -> %s\n",strlen(cList->contig[cList->num].seq),cList->contig[cList->num].seq);
        			free(bdecomp);
        			internb = bread;
    				while(G->read[breadID]->flag != JUNCTION){
    					internb = G->read[breadID]->first;
    					while(internb){
    						if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
//    			    			printf("b %i (%i):\t",G->read[internb->ID]->dir,internb->ID);
    			    			nextoverhang = internb->overhang;
    			    			overhang += internb->overhang;
//    			    			for(j=0;j<overhang;j++){
//    			    				printf(" ");
//    			    			}
    			    			bdecomp = decompressRead(reads[internb->ID].seq,reads[internb->ID].len);
    		        			if(multidir%2==1 && G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
    		        			if(multidir%2==0 && !G->read[internb->ID]->dir) bdecomp = revRead(bdecomp);
//    			    			printf("%s\n",bdecomp);
//    		        			printf("Cat: %i -> %s\n",reads[breadID].len - bread->overhang, &bdecomp[(reads[breadID].len - internb->overhang)]);
//    			    			strcat(cList->contig[cList->num].seq,&bdecomp[reads[internb->ID].len-internb->overhang]);
//    			    			printf("Seq: %i\n",strlen(cList->contig[cList->num].seq));
    			    			free(bdecomp);
    							breadID = internb->ID;
    							break;
    						}

    						internb = internb->next;
    					}
    				}
    				free(backbone);
//    				cList->contig[cList->num].endReadID =  breadID;
//    				cList->num++;
    			}
    			bread = bread->next;
    		}
    	}
    }
    return cList;
}

struct contigList* realBackbone2(struct myovlList* G, struct reads* reads){
	// Take junction read -> transform to POG nodes and edges
	// concatenate the overhangs of the proper flagged reads to the POG
    int i;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    char* backbone;
    printf("MaxReadLen: %i\n",maxReadLen);
    char* readseq = (char*)malloc(sizeof(char)*(maxReadLen+1));
    char* revreadseq = (char*)malloc(maxReadLen+1);

    int dir;
    int bdir;
    int overhang;
    int nextoverhang;
    int multidir;

    struct contigList* cList = (struct contigList*)malloc(sizeof(struct contigList));
    cList->num = 0;
    cList->maxnum = 100;
    cList->contig = (struct contig*)malloc(sizeof(struct contig)*cList->maxnum);

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		printf("JUNCTION found\n");
			dir=G->read[i]->dir;
    		bread = G->read[i]->first;
    		while(bread){
				breadID = bread->ID;
    			if(bread && G->read[breadID]->flag != CONTAINED && bread->dest && bread->dest->len >= MIN_CONTIG_LEN){
//   					printf("Contig Path found: %i -> %i (len: %i) Dir: %i Side: %i\n",i,bread->dest->ID,bread->dest->len,dir,bread->sideflag);
//    				// ContigList is full, resize List
    				if(cList->num+1 == cList->maxnum){
    					printf("ContigList full, Realloc memory\n");
    					cList->maxnum *= 2;
    					cList->contig = (struct contig*)realloc(cList->contig,sizeof(struct contig)*cList->maxnum);
    				}
//    				printf("Test\n");
    				cList->contig[cList->num].stReadID = i;
//    				cList->contig[cList->num].len = bread->dest->len+reads[i].len;
    				backbone = (char*)malloc(bread->dest->len+reads[i].len+1000);
    				if(!backbone){
    					printf("Could not allocate memory for backbone sequence\n");
    					exit(1);
    				}
    				backbone[0] = '\0';
//    				printf("Dirselection\n");
    				if(dir){
    					if(bread->sideflag) 	multidir = 3;
    					else					multidir = 0;
    				}
    				else{
    					if(bread->sideflag)		multidir = 1;
    					else					multidir = 2;
    				}
    				// IF (multidir % 2 == 1) THAN turn b-reads (dir==1) ELSE turn b-reads (dir==0)
    				// IF (multidir > 1) THAN turn a-read
    				bdir = bread->sideflag;
       				decompressReadSt(reads[i].seq,readseq,reads[i].len);
//    				decompressRead(reads[i].seq,readseq,reads[i].len);
       				if(multidir>1){
       					revReadSt(readseq,revreadseq);
//       					printf("   Read: %s\n",readseq);
//       					printf("RevRead: %s\n",revreadseq);
       					strcpy(readseq,revreadseq);
       				}
//     				printf("a %i (%i):\t%s\n",G->read[i]->dir,i,readseq);
       				strcpy(backbone,readseq);
//       				printf("a - Seqlen: %i\n",strlen(backbone));
//        			printf("b %i (%i):\t",G->read[breadID]->dir,breadID);
        			overhang = bread->overhang;
//        			for(j=0;j<overhang;j++){
//        				printf(" ");
//        			}
        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
        			if(multidir%2==1 && G->read[breadID]->dir){
        				revReadSt(readseq,revreadseq);
        				strcpy(readseq,revreadseq);
        			}
        			if(multidir%2==0 && !G->read[breadID]->dir){
            			revReadSt(readseq,revreadseq);
            			strcpy(readseq,revreadseq);
        			}
//        			printf("%s\n",readseq);
//        			printf("Cat: %i -> %s\n",reads[breadID].len - bread->overhang, &readseq[(reads[breadID].len - bread->overhang)]);
        			strcat(backbone,&readseq[(reads[breadID].len - bread->overhang)]);
//        			printf("Seqlen: %i\n",strlen(backbone));
//        			printf("Seq: %i -> %s\n",strlen(cList->contig[cList->num].seq),cList->contig[cList->num].seq);
        			internb = bread;
    				while(G->read[breadID]->flag != JUNCTION){
    					internb = G->read[breadID]->first;
    					while(internb){
    						if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
//    			    			printf("b %i (%i):\t",G->read[internb->ID]->dir,internb->ID);
    			    			nextoverhang = internb->overhang;
    			    			overhang += internb->overhang;
//    			    			for(j=0;j<overhang;j++){
//    			    				printf(" ");
//    			    			}
    			    			breadID = internb->ID;
    		        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
    		        			if(multidir%2==1 && G->read[breadID]->dir){
    		        				revReadSt(readseq,revreadseq);
//    		       					printf("   Read: %s\n",readseq);
//    		       					printf("RevRead: %s\n",revreadseq);
    		        				strcpy(readseq,revreadseq);
    		        			}
    		        			if(multidir%2==0 && !G->read[breadID]->dir){
    		            			revReadSt(readseq,revreadseq);
//    		       					printf("   Read: %s\n",readseq);
//    		       					printf("RevRead: %s\n",revreadseq);
    		            			strcpy(readseq,revreadseq);
    		        			}
//    			    			printf("%s\n",readseq);
//    		        			printf("Cat: %i -> %s\n",reads[breadID].len - internb->overhang, &readseq[(reads[breadID].len - internb->overhang)]);
    			    			strcat(backbone,&readseq[reads[internb->ID].len-internb->overhang]);
//    			    			printf("Seqlen: %i\n",strlen(backbone));
//    			    			printf("Seq: %i\t%s\n",strlen(backbone),backbone);
    							break;
    						}

    						internb = internb->next;
    					}
    				}
    				cList->contig[cList->num].len = strlen(backbone);
    				printf("Len: %i\n",cList->contig[cList->num].len);
    				cList->contig[cList->num].seq = (char*)malloc(strlen(backbone)+1);
    				strcpy(cList->contig[cList->num].seq,backbone);
    				free(backbone);
    				cList->contig[cList->num].endReadID =  breadID;
    				cList->num++;
    			}
    			bread = bread->next;
    		}
    	}
    }

    free(readseq);
    free(revreadseq);

    return cList;
}

//LPOLetter_T* LPOnodes  		=  NULL;   	// LPOLetter_T array -> First malloc in make_poa
//uint32_t     usedNodes 		=     0;	// Current number of elements in LPOLetter_T array
//uint32_t     maxNodes  		= 10000;	// Initial max size of LPOLetter_T array
//uint32_t 	 contigNum 		=     0;	// Number of contigs
//uint32_t 	 max_contigNum 	=   100;	// Max number of contigs before resize the struct (Sequence_T)
//
//void poa_initBackbone(Sequence_T* contig, struct reads* read, char* seq){
//	int i;
//	uint32_t currentID = usedNodes;
//	uint32_t leftID = currentID;
//
//	contig->length = strlen(seq);
//	contig->nsource_seq = 1;
//	contig->letter = &LPOnodes[currentID];
//
//	LPOLetter_T* current = contig->letter;
//	currentID++;
//
//	LPOLetter_T* left = current;
//	current = &LPOnodes[currentID];
//
//	for(i=0;i<contig->length;i++){
//		current->letter = seq[i];
//		current->source.ipos = i;
//		current->source.iseq = read->ID;
//		current->source.more = NULL;
//
//		current->left.ipos = leftID;
//		current->left.more = NULL;
//		left->right.ipos = currentID;
//		left->right.more = NULL;
//
//		left = current;
//		leftID = currentID;
//		currentID++;
//		current = &LPOnodes[currentID];
//	}
//}

struct Letter_T* Letters     =      NULL;
uint32_t          numNodes    =         0;	// Current number of elements in LPOLetter_T array
uint32_t          maxNumNodes =  10000000;	// Initial max size of LPOLetter_T array

int **alMatrix = NULL;
struct Letter_T** alMatrix_Letter = NULL;

/**
 * Check size of Letters-array and resize if full
 */
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
				Letters[j].source.next = NULL;
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


/**
 *	Initialize the linear backbone poa-graph, beginning with the junction read. Connect with the contig struct.
 * @param contig	The Contig initially beginning with this given junction read
 * @param read		The junction read representing the start of a contig
 * @param seq		Sequence of the junction read
 */
void poa_initBackbone2(struct Sequence* contig, struct reads* read, char* seq){
	int i;
	uint32_t leftID;

	contig->length = strlen(seq);
	contig->nsource_seq = 1;

	struct Letter_T* left = NULL;
	struct Letter_T* current;
	contig->startLetter.dest = numNodes;
	contig->startLetter.next = NULL;
	contig->readleft = numNodes;

	for(i=0;i<contig->length;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->source.ipos = i;
		current->source.iseq = read->ID;
		current->source.next = NULL;
		current->counter = 1;

		if(left){
			current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
			current->left->counter = 1;
			current->left->dest = leftID;
			current->left->next = NULL;
//			printf("LeftID: %i\n",leftID);
			left = &Letters[leftID];
			left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
			left->right->counter = 1;
			left->right->dest = numNodes;
			left->right->next = NULL;
		}

		left = current;
		leftID = numNodes;
		numNodes++;
		poa_LetterSizeCheck();
	}
	contig->readright = numNodes-1;
	printf("Number of set numNodes: %i\n",numNodes);
}

void poa_catBackbone(struct Sequence* contig, struct myovlList *G, struct reads* read, char* seq, int leftID, int rightID){
	struct bread* leftB;
	struct bread* rightB;

//	printf("LeftID: %i\n",leftID);

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

	int left_ovh, right_ovh;
	left_ovh = leftB->overhang;
	right_ovh = rightB->overhang;
	contig->length += right_ovh;

	int len = strlen(seq);
	int i =  len - right_ovh;
	leftID = contig->readright;
	struct Letter_T* left;
	struct Letter_T* current;

//	printf("LeftID: %i\n",leftID);

	for(;i<len;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->counter=1;
		// _________
		// This things not until alignment
		current->source.ipos = i;
		current->source.iseq = read->ID;
		current->source.next = NULL;
		// _________
//		printf("NewLetter: %i\n",numNodes);
		current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		current->left->counter = 1;
		current->left->dest = leftID;
		current->left->next = NULL;
//		printf("Set Right edge to leftID: %i\n",leftID);
		left = &Letters[leftID];
		left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		left->right->counter = 1;
		left->right->dest = numNodes;
		left->right->next = NULL;

		leftID = numNodes;
		numNodes++;
		poa_LetterSizeCheck();
	}
	contig->readright = numNodes-1;

	current = &Letters[contig->readleft];
	struct LetterEdge* edge = NULL;
	for(i=0;i<left_ovh;i++){
		edge = current->right;
		if(edge){
			current = &Letters[edge->dest];
		}
		else{
			printf("Letter has no edge!!!\n");
			exit(1);
		}
	}
	contig->readleft = edge->dest;
}


void poa_toDot(char* dotFile){
	FILE* dot = fopen(dotFile,"w");

	int i;
	struct Letter_T* current;
	struct LetterEdge* edge;
	fprintf(dot,"digraph POMSA {\n"); // Directedgraph - PartialOrderedMultipleSequenceAlignment

	for(i=0;i<numNodes;i++){
		Letters[i].junction = 1;
	}

	int cluster = 0;
	int id;

	for(i=0;i<numNodes;i++){
		current = &Letters[i];
		if(current->junction == 1){
			if(current->align_ring){
				fprintf(dot,"\tsubgraph cluster_%i{\n",cluster);
				fprintf(dot,"\t\tlabel=\"%i\";\n",cluster++);
//				fprintf(dot,"\t\tlabel=\" \";\n");
//				fprintf(dot,"\tsubgraph {\n");
//				fprintf(dot,"\t\trank=same;\n");
				while(current->junction){
					current->junction = 0;
					id = current - Letters;
					if(id<0) id *= -1;
					fprintf(dot,"\t\t%i [label=\"%c: %i\"];\n",id,current->letter,current->counter);
					current = current->align_ring;
				}
				fprintf(dot,"\t}\n");
			}
			else fprintf(dot,"\t%i [label=\"%c: %i\"];\n",i,current->letter,current->counter);
		}

		edge = current->right;
		while(edge){
			fprintf(dot,"\t%i -> %i [label=%i];\n",i,edge->dest,edge->counter);
			edge = edge->next;
		}
	}
	fprintf(dot,"}\n");

	fclose(dot);
}

void poa_part_toDot(char* dotFile, struct Sequence* contig){
	printf("CHECKPOINT: Write Suspect Part of the POA Graph\n");

	FILE* dot = fopen(dotFile,"w");

//	struct Letter_T* current = &Letters[contig->readleft];

	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	fprintf(dot,"digraph POMSA {\n"); // Directedgraph - PartialOrderedMultipleSequenceAlignment
	int cluster = 0;
	int id;

	struct Letter_T** newLetters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000);
	int newNum = 0;
	struct Letter_T** oldLetters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000);
	int oldNum = 0;

	newLetters[newNum++] = current;

	char* edgecharp = (char*)malloc(10000);
	char* edgechar = (char*)malloc(10000);

	while(newNum){
//		printf("While: newNum = %i\n",newNum);
		memcpy(oldLetters,newLetters,sizeof(struct Letters_T*)*newNum);
		oldNum = newNum;
		newNum = 0;

		while(oldNum){
			current = oldLetters[--oldNum];
			if(!current->vFlag){
				if(current->align_ring){
					fprintf(dot,"\tsubgraph cluster_%i{\n",cluster);
					fprintf(dot,"\t\tlabel=\"%i\";\n",cluster++);
					*edgechar = '\0';
					while(!current->vFlag){
						id = current - Letters;
						if(id<0) id *= -1;
						fprintf(dot,"\t\t%i [label=\"%c: %i\"];\n",id,current->letter,current->counter);
						current->vFlag = 1;

						edge = current->right;
						while(edge){
							sprintf(edgecharp,"\t%i -> %i [label=%i];\n",id,edge->dest,edge->counter);
							strcat(edgechar,edgecharp);
//							fprintf(dot,"\t%i -> %i [label=%i];\n",id,edge->dest,edge->counter);
							if(!Letters[edge->dest].vFlag){
								newLetters[newNum++] = &Letters[edge->dest];
							}
							edge = edge->next;
						}
						current = current->align_ring;
					}
					fprintf(dot,"\t}\n");
					fprintf(dot,"%s",edgechar);
				}
				else{
					id = current - Letters;
					if(id<0) id *= -1;
					if(current == &Letters[contig->readright]){
						fprintf(dot,"\t%i [label=\"%c: %i\",color=green];\n",id,current->letter,current->counter);
					}
					else{
						fprintf(dot,"\t%i [label=\"%c: %i\"];\n",id,current->letter,current->counter);
					}
					current->vFlag = 1;
					edge = current->right;
					while(edge){
						fprintf(dot,"\t%i -> %i [label=%i];\n",id,edge->dest,edge->counter);
						if(!Letters[edge->dest].vFlag){
							newLetters[newNum++] = &Letters[edge->dest];
						}
						edge = edge->next;
					}
				}

			}
		}
//		printf("While: newNum = %i\n",newNum);
	}

	fprintf(dot,"}\n");

	fclose(dot);

	current = &Letters[contig->startLetter.dest];
	newNum = 0;
	oldNum = 0;
	newLetters[newNum++] = current;
	while(newNum){
//		printf("While: newNum = %i\n",newNum);
		memcpy(oldLetters,newLetters,sizeof(struct Letters_T*)*newNum);
		oldNum = newNum;
		newNum = 0;
		while(oldNum){
			current = oldLetters[--oldNum];
			if(current->vFlag){
				if(current->align_ring){
					while(current->vFlag){
						current->vFlag = 0;
						edge = current->right;
						while(edge){
							if(Letters[edge->dest].vFlag){
								newLetters[newNum++] = &Letters[edge->dest];
							}
							edge = edge->next;
						}
						current = current->align_ring;
					}
				}
				else{
					current->vFlag = 0;
					edge = current->right;
					while(edge){
						if(Letters[edge->dest].vFlag){
							newLetters[newNum++] = &Letters[edge->dest];
						}
						edge = edge->next;
					}
				}
			}
		}
	}
	printf("Part writing successful\n");
//	exit(1);
}

void poa_consensus(struct Sequence* contig){
	printf("CHECKPOINT: PO-MSA to Contig\n");

	char* seq = (char*)malloc(contig->length + 100);
	int i=0;

	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	struct LetterEdge* bestedge;


	while(1){
		if(current->counter < 5) seq[i++] = current->letter+32;
		else seq[i++] = current->letter;
		edge = current->right;
		if(edge){
			bestedge = edge;
			while(edge){
				if(bestedge->counter < edge->counter){
					bestedge = edge;
				}
				edge = edge->next;
			}
			current = &Letters[bestedge->dest];
		}
		else{
			break;
		}
	}

	seq[i]='\0';
//	printf(">Correct_%s\n",contig->name);
//	int k;
//	for(k=0;k<i;k+=80){
//		printf("%.80s\n",&seq[k]);
//	}
//	printf("\n");
	contig->sequence = seq;
}

void poa_printContigs(struct POG* pog, char* contigFile){
	printf("CHECKPOINT: Write CorrectContigs in fasta\n");
	FILE* correctContigs = fopen(contigFile,"w");

	int i,j;
	int len;

	for(i=0;i<pog->contigNum;i++){
		len = strlen(pog->contig[i].sequence);
		fprintf(correctContigs,">%s_%i\n",pog->contig[i].name,len);
		for(j=0;j<len;j+=80){
			fprintf(correctContigs,"%.80s\n",&pog->contig[i].sequence[j]);
		}
	}
	fclose(correctContigs);
}

int findLeftMostJunction(int i){
	int leftMost = i;

    struct pathEdge* pathEdge;
    struct pathEdge* tempEdge;
    int tempdepth;

    // Path touring, run through the junctions over spanned by a certain (which threshold?) number of read pairs
    for(i=1;i<pathsNum;){
    	if(!paths[i].flag){
    		// walk through the graph till a left most unvisted path was found. Intelligent touring, not so easy.  Think intensively
        	pathEdge = paths[i].leftPath;
        	if(pathEdge){
        		tempdepth = 0;
        		while(pathEdge){
        			if(tempdepth < pathEdge->depth){
        				tempEdge = pathEdge;
        			}
        		}
        	}
    	}
    	else{
    		i++;
    	}

    }

    return leftMost;
}


struct POG* make_poa(struct myovlList* G, struct reads* reads){
	// Take junction read -> transform to POG nodes and edges
	// concatenate the overhangs of the proper flagged reads to the POG
	char verbose = 0;
    int i,j;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    printf("MaxReadLen: %i\n",maxReadLen);
    char* readseq = (char*)malloc(sizeof(char)*(maxReadLen+1));
    char* revreadseq = (char*)malloc(maxReadLen+1);
    char* name = (char*)malloc(1000);

//    Sequence_T *contigs = (Sequence_T*)malloc(sizeof(Sequence_T) * max_contigNum);
//    LPOnodes = (LPOLetter_T*)malloc(sizeof(LPOLetter_T)*maxNodes);

    // Init Graph structure
    struct POG* pog = (struct POG*)malloc(sizeof(struct POG));
    pog->contigNum = 0;
    pog->maxNum = 100;
    pog->contig = (struct Sequence*)malloc(sizeof(struct Sequence)*pog->maxNum);
    if(!Letters) Letters = (struct Letter_T*)malloc(sizeof(struct Letter_T)*maxNumNodes);
    for(i=0;i<maxNumNodes;i++){
    	Letters[i].left = NULL;
    	Letters[i].right = NULL;
    	Letters[i].junction = 0;
    	Letters[i].source.next = NULL;
    	Letters[i].score = 0;
    	Letters[i].vFlag = 0;
    }

    // Init Matrix (5 x maxreaden * maxreadlen)
    // E.g. for maxReadLen = 100 -> 500 x 100 matrix (50,000 array)
    alMatrix = (int**)malloc(sizeof(int*)*(maxReadLen*MATRIX_MAX_BR+1)); // Convention that the aligning part of the graph do not contain more than 5*maxReadLen nodes
    alMatrix_Letter = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*(maxReadLen*MATRIX_MAX_BR+1));
    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
    	alMatrix[i]=(int*)malloc(sizeof(int)*(maxReadLen+1));
    	alMatrix_Letter[i] = NULL;
    	for(j=0;j<=maxReadLen;j++){
//    		alMatrix[i][j] = (i+j) * GAP_PENALTY;
    		alMatrix[i][j] = j * GAP_PENALTY;
    	}
    }

    int dir;
    int bdir;
    int overhang;
    int nextoverhang;
    int multidir;
    int oldbreadID;
    int inserts = 0;
    struct bread* counterbread;

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
//    		printf("JUNCTION found\n");
			dir=G->read[i]->dir;
    		bread = G->read[i]->first;
    		while(bread){
				breadID = bread->ID;
    			if(bread && G->read[breadID]->flag != CONTAINED && bread->dest && !bread->dest->flag && bread->dest->len >= MIN_CONTIG_LEN){

    				// Flag the other side of the unique path as already used, to not report the reverse path
    				counterbread =  G->read[bread->dest->ID]->first;
    				while(counterbread){
    					if(counterbread->dest && counterbread->dest->ID == i){
    						printf("Flag reverse path\n");
    						counterbread->dest->flag = 1;
    						break;
    					}
    					counterbread = counterbread->next;
    				}

   					printf("Contig Path found: %i -> %i (len: %i) Dir: %i Side: %i\n",i,bread->dest->ID,bread->dest->len,dir,bread->sideflag);
//    				// ContigList is full, resize List
    				if(pog->contigNum == pog->maxNum){
    					printf("ContigList full, Realloc memory\n");
    					pog->contig = (struct Sequence*)realloc(pog->contig,sizeof(struct Sequence)*(pog->maxNum*2));
    					pog->maxNum *= 2;
    				}
    				pog->contig[pog->contigNum].nsource_seq = 0;
    				pog->contig[pog->contigNum].length = 0;
    				sprintf(name,"%i_%i_%i_%i",bread->dest->pathID,pog->contigNum,i,bread->dest->ID);
    				pog->contig[pog->contigNum].name = (char*)malloc(strlen(name)+1);
    				strcpy(pog->contig[pog->contigNum].name,name);
    				if(dir){
    					if(bread->sideflag) 	multidir = 3;
    					else					multidir = 0;
    				}
    				else{
    					if(bread->sideflag)		multidir = 1;
    					else					multidir = 2;
    				}
    				// IF (multidir % 2 == 1) THAN turn b-reads (dir==1) ELSE turn b-reads (dir==0)
    				// IF (multidir > 1) THAN turn a-read
    				bdir = bread->sideflag;
       				decompressReadSt(reads[i].seq,readseq,reads[i].len);
       				if(multidir>1){
       					revReadSt(readseq,revreadseq);
       					strcpy(readseq,revreadseq);
       				}
//       				poa_initBackbone(&contigs[contigNum],&reads[i],readseq);
       				poa_initBackbone2(&pog->contig[pog->contigNum],&reads[i],readseq);
        			overhang = bread->overhang;
        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
        			if(multidir%2==1 && G->read[breadID]->dir){
        				revReadSt(readseq,revreadseq);
        				strcpy(readseq,revreadseq);
        			}
        			if(multidir%2==0 && !G->read[breadID]->dir){
            			revReadSt(readseq,revreadseq);
            			strcpy(readseq,revreadseq);
        			}
        			poa_catBackbone(&pog->contig[pog->contigNum],G,&reads[breadID],readseq,i,breadID);
        			internb = bread;
//        			if(pog->contigNum == 4){
        				while(G->read[breadID]->flag != JUNCTION){
        					internb = G->read[breadID]->first;
        					while(internb){
        						if(G->read[internb->ID]->flag == CONTAINED){
        		        			decompressReadSt(reads[internb->ID].seq,readseq,reads[internb->ID].len);
        		        			if(multidir%2==1 && !G->read[internb->ID]->dir){
        		        				revReadSt(readseq,revreadseq);
        		        				strcpy(readseq,revreadseq);
        		        			}
        		        			if(multidir%2==0 && G->read[internb->ID]->dir){
        		            			revReadSt(readseq,revreadseq);
        		            			strcpy(readseq,revreadseq);
        		        			}
        		        			if(verbose) printf("ALINGING CONTAINED READ\n");
        		        			if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
        		        			if(verbose) printf("Read: %s\n",readseq);
//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
        		        			poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[internb->ID],readseq,0,inserts,0);
//        		        			inserts++;
        						}
        						internb = internb->next;
        					}
        					internb = G->read[breadID]->first;
        					while(internb){
        						if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
        							nextoverhang = internb->overhang;
        			    			overhang += internb->overhang;
        			    			oldbreadID = breadID;
        			    			breadID = internb->ID;
        		        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
        		        			if(multidir%2==1 && G->read[breadID]->dir){
        		        				revReadSt(readseq,revreadseq);
        		        				strcpy(readseq,revreadseq);
        		        			}
        		        			if(multidir%2==0 && !G->read[breadID]->dir){
        		            			revReadSt(readseq,revreadseq);
        		            			strcpy(readseq,revreadseq);
        		        			}
            						poa_catBackbone2(&pog->contig[pog->contigNum],G,&reads[breadID],readseq,oldbreadID,breadID);
    //        						poa_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts);
//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
            						if(verbose) printf("ALINGING PROPER READ\n");
            						if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
        		        			if(verbose) printf("Read: %s\n",readseq);
            						poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts,nextoverhang);
    //        						if(inserts == 290585)	poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,4077);
    //        						else poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,1);
            						inserts++;
    //        						if(inserts == 2) return pog;
        							break;
        						}
        						internb = internb->next;
        					}
        				}
//        			}
//        			if(pog->contigNum == 4)	poa_consensus(&pog->contig[pog->contigNum]);
        			poa_consensus(&pog->contig[pog->contigNum]);
    				pog->contigNum++;
    				printf("Matrix time:    %.3f s\n",(float)sumMatrix/1000000000);
    				printf("Backtrace time: %.3f s\n",(float)sumTrace/1000000000);
    				poa_toDot("output/poa.dot");
//    				exit(1);
    			}
    			bread = bread->next;
    		}
    	}
    }

    free(readseq);
    free(revreadseq);
    free(name);

    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
    	free(alMatrix[i]);
    }
    free(alMatrix);

    return pog;
}

//void testFunct(){
//    int i,j;
//    printf("MaxReadLen: %i\n",maxReadLen);
//    maxReadLen = 100;
//
//    Sequence_T *contigs = (Sequence_T*)malloc(sizeof(Sequence_T) * max_contigNum);
//    LPOnodes = (LPOLetter_T*)malloc(sizeof(LPOLetter_T)*maxNodes);
//
//    // Init Graph structure
//    struct POG* pog = (struct POG*)malloc(sizeof(struct POG));
//    pog->contigNum = 0;
//    pog->maxNum = 100;
//    pog->contig = (struct Sequence*)malloc(sizeof(struct Sequence)*pog->maxNum);
//    if(!Letters) Letters = (struct Letter_T*)malloc(sizeof(struct Letter_T)*maxNumNodes);
//    for(i=0;i<maxNumNodes;i++){
//    	Letters[i].left = NULL;
//    	Letters[i].right = NULL;
//    	Letters[i].junction = 0;
//    	Letters[i].source.next = NULL;
//    	Letters[i].score = 0;
//    }
//
//    // Init Matrix (5 x maxreaden * maxreadlen)
//    // E.g. for maxReadLen = 100 -> 500 x 100 matrix (50,000 array)
//    alMatrix = (int**)malloc(sizeof(int*)*(maxReadLen*MATRIX_MAX_BR+1)); // Convention that the aligning part of the graph do not contain more than 5*maxReadLen nodes
//    alMatrix_Letter = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*(maxReadLen*MATRIX_MAX_BR+1));
//    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
//    	alMatrix[i]=(int*)malloc(sizeof(int)*(maxReadLen+1));
//    	alMatrix_Letter[i] = NULL;
//    	for(j=0;j<=maxReadLen;j++){
//    		alMatrix[i][j] = j * GAP_PENALTY;
//    	}
//    }
//	pog->contig[pog->contigNum].length = 0;
//
//	struct reads* reads = (struct reads*)malloc(sizeof(struct reads)*3);
//	reads[0].ID = 0;
//	reads[1].ID = 1;
//	reads[2].ID = 2;
//	reads[0].len = 50;
//	reads[1].len = 50;
//	reads[2].len = 40;
//	reads[0].seq = (char*)malloc(51);
//	reads[1].seq = (char*)malloc(51);
//	reads[2].seq = (char*)malloc(41);
//	strcpy(reads[0].seq,"GAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGCCAC");
//	                   strcpy(reads[1].seq,"AGCAAATTAAAATTTTATTGACATAGGTCACTAAATACTTGATCCAATAT");
//	                         strcpy(reads[2].seq,"TTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCA");
//
//	struct myovlList *G = (struct myovlList*)malloc(sizeof(struct myovlList));
//	G->read = (struct aread**)malloc(sizeof(struct aread*)*3);
//	G->V = 3;
//	G->read[0] = (struct aread*)malloc(sizeof(struct aread));
//	G->read[0]->flag = PROPER;
//	G->read[0]->dir = 1;
//	G->read[0]->first = (struct bread*)malloc(sizeof(struct bread));
//	G->read[0]->first->ID = 1;
//	G->read[0]->first->overhang = 19;
//	G->read[0]->first->next = NULL;
//	G->read[1] = (struct aread*)malloc(sizeof(struct aread));
//	G->read[1]->flag = PROPER;
//	G->read[1]->dir = 1;
//	G->read[1]->first = (struct bread*)malloc(sizeof(struct bread));
//	G->read[1]->first->ID = 0;
//	G->read[1]->first->overhang = 19;
//	G->read[1]->first->next = (struct bread*)malloc(sizeof(struct bread));
//	G->read[1]->first->next->ID = 2;
//	G->read[1]->first->next->overhang = -4;
//	G->read[1]->first->next->next = NULL;
//	G->read[2] = (struct aread*)malloc(sizeof(struct aread));
//	G->read[2]->flag = PROPER;
//	G->read[2]->dir = 1;
//	G->read[2]->first = (struct bread*)malloc(sizeof(struct bread));
//	G->read[2]->first->ID = 1;
//	G->read[2]->first->overhang = 6;
//	G->read[2]->first->next = NULL;
//
//	poa_initBackbone2(&pog->contig[pog->contigNum],&reads[0],reads[0].seq);
//	poa_catBackbone2(&pog->contig[pog->contigNum],G,&reads[1],reads[1].seq,0,1);
//	printf("Insert First read\n");
//	poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[1],reads[1].seq,1,0);
//	poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[2],reads[2].seq,0,0);
//
//	pog->contigNum++;
//	char* contigPath = "./output/correct_contigs.fasta";
//	poa_printContigs(pog,contigPath);
//	poa_toDot("output/poa.dot");
//
//	exit(1);
//
//
//}

void scaffold_stats(struct scaffold_set* aS){
    int gesLen = 0;													// Sum over all scaffold length
	int *nStat = (int*)malloc(sizeof(int)*aS->num);				// List of Scaffold length

	int i;

    int v;
    int k=0;
    for (v = 0; v < aS->num; v++){
    	nStat[k++] = aS->scaff[v].len;
    }

	int temp;
	k = aS->num;
	int newk;
	// Sort Scaffolds, simple bubble-sort
	do{
		newk = 1;
		for(i = 0; i < k-1 ; ++i){
			if(nStat[i] < nStat[i+1]){
				temp = nStat[i+1];
				nStat[i+1] = nStat[i];
				nStat[i] = temp;
				newk = i+1;
			}
		}
		k = newk;
	} while(k>1);

	for(i=0;i<aS->num;i++) gesLen += nStat[i];

	int sum=0;
	int ns=0;

	struct scaffEdge* scaffEdge;
	int startJunction;
	printf("Scaffold Paths -> num: %i\n",aS->num);

	for(i=0;i<aS->num;i++){
//		if(aS->scaff[i].first){
			if(aS->scaff[i].first->targetJunction == paths[aS->scaff[i].first->ID].rightJunction){
				startJunction = paths[aS->scaff[i].first->ID].leftJunction;
			}
			else{
				startJunction = paths[aS->scaff[i].first->ID].rightJunction;
			}
			printf("Scaffold: %i (len: %i bp)\n",i,aS->scaff[i].len);
			scaffEdge = aS->scaff[i].first;
			printf(KRED"%i"KNRM,startJunction);
			while(scaffEdge){
				printf(" -> "KGRN"%i"KNRM,scaffEdge->ID);
				printf(" -> "KRED"%i"KNRM,scaffEdge->targetJunction);
				scaffEdge = scaffEdge->next;
			}
			printf("\n");
//		}
//		else{
//			printf("Scaffold: %i (len: %i bp)\n",i,aS->scaff[i].len);
//			printf(KRED"%i"KNRM,aS->scaff[i].startJunction);
//			printf(" -> "KGRN"%i"KNRM,aS->scaff[i].ID);
//			printf(" -> "KRED"%i"KNRM,aS->);
//
//		}

	}

	printf("\n");
	printf("Largest Scaffold: %i bp\n",nStat[0]);
	printf("Number of Scaffolds: %i\n",aS->num);
	printf("Total Length over all Scaffolds: %i\n",gesLen);
	for(i=0;i<aS->num;i++){
		sum += nStat[i];
		if(sum > (gesLen/10) && ns == 0){
			printf("N10: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (gesLen/4) && ns == 1){
			printf("N25: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (gesLen/2) && ns == 2){
			printf("N50: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (gesLen/4)*3 && ns == 3){
			printf("N75: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (gesLen/10)*9 && ns == 4){
			printf("N90: %i\n\n",nStat[i]);
			ns++;
		}
	}

}


struct pathEdge** uniqPath;
int uniqlenmax = 100;

void scaffold_uniqueStringTouring(int i, int rightdir, struct scaffold* scaff){
	paths[i].flag = 1;

	struct pathEdge* pathEdge;
	// Search for paths without a scaffold connection or with ambiguous connection to one side and tour to the other one
	printf("Start point for Path Scaffolding Found: %i\n",i);
	printf("\t -> %i Length: %i (Startpath)\n",i,paths[i].len);


	int uniqlen = 0;

	struct scaffEdge* scaffEdge;
	struct scaffEdge* ancestorEdge;
	int depth = 0;

	scaffEdge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
	scaffEdge->ID = i;
	scaffEdge->depth = depth++;
	scaffEdge->next = NULL;
	ancestorEdge = scaffEdge;
	scaff->first = scaffEdge;
	scaff->len = paths[i].len;
	scaff->type = 0;

	if(!rightdir){
		pathEdge = paths[i].leftPath;
		scaff->startJunction = paths[i].rightJunction;
		scaff->endJunction = paths[i].leftJunction;
		scaffEdge->targetJunction = paths[i].leftJunction;
	}
	else{
		pathEdge = paths[i].rightPath;
		scaff->startJunction = paths[i].leftJunction;
		scaff->endJunction = paths[i].rightJunction;
		scaffEdge->targetJunction = paths[i].rightJunction;
	}
	while(pathEdge){
		if(pathEdge->sibl){
			printf("Path becomes ambiguous!\n");
			break;
		}
		printf("\t -> %i (c: %i) Target: %i Length: %i\n",pathEdge->ID,pathEdge->counter,pathEdge->targetJunction,paths[pathEdge->ID].len);
		scaff->len += paths[pathEdge->ID].len;
		scaff->type = 1;
		scaff->endJunction = pathEdge->targetJunction;
		paths[pathEdge->ID].flag = 1;
		scaffEdge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
		scaffEdge->ID = pathEdge->ID;
		scaffEdge->depth = depth++;
		scaffEdge->targetJunction = pathEdge->targetJunction;
		scaffEdge->next = NULL;
		ancestorEdge->next = scaffEdge;
		ancestorEdge = scaffEdge;

		uniqPath[uniqlen++] = pathEdge;
		if(uniqlen == uniqlenmax){
			uniqlenmax *= 2;
			uniqPath = (struct pathEdge**)realloc(uniqPath,sizeof(struct pathEdge*)*uniqlenmax);
			if(!uniqPath){
				printf("Could not reallocate uniqPath\n Abort\n");
				exit(1);
			}
		}
		pathEdge = pathEdge->next;
	}

	int j=0;
	int k;

	while(j < uniqlen){
		paths[uniqPath[j]->ID].flag = 1;
		k=j+1;

		printf("Try to elongate from path: %i (Target: %i) right: %i left: %i\n",uniqPath[j]->ID,uniqPath[j]->targetJunction,paths[uniqPath[j]->ID].rightJunction,paths[uniqPath[j]->ID].leftJunction);
		if(paths[uniqPath[j]->ID].rightJunction == uniqPath[j]->targetJunction){
			pathEdge = paths[uniqPath[j]->ID].rightPath;
		}
		else if(paths[uniqPath[j]->ID].leftJunction == uniqPath[j]->targetJunction){
			pathEdge = paths[uniqPath[j]->ID].leftPath;
		}
		else{
			printf("No side found! Abort!\n");
			exit(1);
		}
		// search end of the recorded path for probable elongation
		while(k < uniqlen && pathEdge){
			if(uniqPath[k] == pathEdge){
				pathEdge = pathEdge->next;
			}
			else{
				pathEdge = pathEdge->sibl;
				while(pathEdge){
					if(uniqPath[k] == pathEdge){
						pathEdge = pathEdge->next;
						break;
					}
					pathEdge = pathEdge->sibl;
				}
			}
			k++;
		}
		j++;
		// now try to elongate
//		printf("Elongate\n");
		while(pathEdge){
			if(pathEdge->sibl){
				printf("Path becomes ambiguous!\n");
				break;
			}
			printf("\t -> %i (c: %i) Length: %i (Elong)\n",pathEdge->ID,pathEdge->counter,paths[pathEdge->ID].len);
			scaff->len += paths[pathEdge->ID].len;
			scaff->endJunction = pathEdge->targetJunction;
			paths[pathEdge->ID].flag = 1;
			scaffEdge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
			scaffEdge->ID = pathEdge->ID;
			scaffEdge->depth = depth++;
			scaffEdge->targetJunction = pathEdge->targetJunction;
			scaffEdge->next = NULL;
			ancestorEdge->next = scaffEdge;
			ancestorEdge = scaffEdge;
			uniqPath[uniqlen++] = pathEdge;

			if(uniqlen == uniqlenmax){
				uniqlenmax *= 2;
				uniqPath = (struct pathEdge**)realloc(uniqPath,sizeof(struct pathEdge*)*uniqlenmax);
				if(!uniqPath){
					printf("Could not reallocate uniqPath\n Abort\n");
					exit(1);
				}
			}
			pathEdge = pathEdge->next;
		}
	}
}


struct scaffold_set* scaffold_init(){
    int i;
    struct scaffold_set* aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
    aS->num = 0;
    aS->nummax = 1000;
    aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*aS->nummax);

    // 2 loops for graph touring:
    uniqPath = (struct pathEdge**)malloc(sizeof(struct pathEdge*)*uniqlenmax);

    for(i=1;i<pathsNum;i++){
    	if(!paths[i].flag){
    		// Left end of a path
        	if(!paths[i].leftPath && paths[i].rightPath){
        		scaffold_uniqueStringTouring(i,1,&aS->scaff[aS->num]);
        		aS->num++;
        	}
        	// right end of a path
        	else if(!paths[i].rightPath && paths[i].leftPath){
        		scaffold_uniqueStringTouring(i,0,&aS->scaff[aS->num]);
        		aS->num++;
        	}
        	// left side is ambiguous
        	else if(paths[i].leftPath && paths[i].leftPath->sibl){
//        		scaffold_uniqueStringTouring(i,1,&aS->scaff[aS->num]);
//        		aS->num++;
        	}
        	// right side is ambiguous
        	else if(paths[i].rightPath && paths[i].rightPath->sibl){
//        		scaffold_uniqueStringTouring(i,0,&aS->scaff[aS->num]);
//        		aS->num++;
        	}
        	else{
//        		scaffold_uniqueStringTouring(i,1,&aS->scaff[aS->num]);
//        		aS->num++;
        	}
        	if(aS->num == aS->nummax){
        		aS->nummax *= 2;
        		aS->scaff = (struct scaffold*)realloc(aS->scaff,sizeof(struct scaffold)*aS->nummax);
        		if(!aS->scaff){
        			printf("No realloc of struct scaffold possible\n Abort\n");
        			exit(1);
        		}
        	}
    	}

    }

    // 2. Search for paths not flagged after the first touring and call consensus by treating them as singletons

	return aS;
}

/**
 * Same implementation of the pao algorithms but graph touring over the paths instead of the junctions. Spanning junctions of scaffolds span over
 * @return
 */
struct POG* make_poaScaff(struct myovlList* G, struct reads* reads, char scaffolding){

	struct scaffold_set* aS;
	if(scaffolding){
		printf("Checkpoint: Init Scaffold Correction\n");
		aS = scaffold_init();
	}
	else{
		printf("Checkpoint: Init Contig Correction\n");
		aS = contigs_init(G); // ,reads
	}
    scaffold_stats(aS);

	char verbose = 0;
    int i,j;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    printf("MaxReadLen: %i\n",maxReadLen);
    char* readseq = (char*)malloc(sizeof(char)*(maxReadLen+1));
    char* revreadseq = (char*)malloc(maxReadLen+1);
    char* name = (char*)malloc(1000);

//    Sequence_T *contigs = (Sequence_T*)malloc(sizeof(Sequence_T) * max_contigNum);
//    LPOnodes = (LPOLetter_T*)malloc(sizeof(LPOLetter_T)*maxNodes);

    // Init Graph structure
    struct POG* pog = (struct POG*)malloc(sizeof(struct POG));
    pog->contigNum = 0;
    pog->maxNum = 100;
    pog->contig = (struct Sequence*)malloc(sizeof(struct Sequence)*pog->maxNum);
    if(!Letters) Letters = (struct Letter_T*)malloc(sizeof(struct Letter_T)*maxNumNodes);
    for(i=0;i<maxNumNodes;i++){
    	Letters[i].left = NULL;
    	Letters[i].right = NULL;
    	Letters[i].junction = 0;
    	Letters[i].source.next = NULL;
    	Letters[i].score = 0;
    	Letters[i].vFlag = 0;
    }

    // Init Matrix (5 x maxreaden * maxreadlen)
    // E.g. for maxReadLen = 100 -> 500 x 100 matrix (50,000 array)
    alMatrix = (int**)malloc(sizeof(int*)*(maxReadLen*MATRIX_MAX_BR+1)); // Convention that the aligning part of the graph do not contain more than 5*maxReadLen nodes
    alMatrix_Letter = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*(maxReadLen*MATRIX_MAX_BR+1));
    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
    	alMatrix[i]=(int*)malloc(sizeof(int)*(maxReadLen+1));
    	alMatrix_Letter[i] = NULL;
    	for(j=0;j<=maxReadLen;j++){
//    		alMatrix[i][j] = (i+j) * GAP_PENALTY;
    		alMatrix[i][j] = j * GAP_PENALTY;
    	}
    }

    int dir;
    int bdir;
    int overhang;
    int nextoverhang;
    int multidir;
    int oldbreadID;
    int inserts = 0;
//    struct bread* counterbread;


    //__________ Scaffold Consensus Calling
    struct scaffEdge* scaffEdge;
    int startJunction;
    for(i=0;i<aS->num;i++){
    	if(aS->scaff[i].len > MIN_SCAFF_LEN){
    		scaffEdge = aS->scaff[i].first;
    		printf("FirstEdge: %i\n",scaffEdge->ID);
    		startJunction = aS->scaff[i].startJunction;
    		dir = G->read[startJunction]->dir;
    		bread = G->read[startJunction]->first;
    		// Scaffold is a Singleton; Just on Contig (no scaffEdge)
    		if(aS->scaff[i].type == 1){
    			while(bread){
    				if(bread->dest) printf("destID: %i == endJunction: %i\n",bread->dest->ID, aS->scaff[i].endJunction);
    				if(bread->dest && bread->dest->ID == aS->scaff[i].endJunction) break;
    				else bread = bread->next;
    			}
//    			if(bread->dest->ID != 284801){
//    				continue;
//    			}
    		}
    		else{
        		while(bread){
        			printf("destPathID: %i == scaffEdgeID: %i\n",bread->dest->pathID, scaffEdge->ID);
        			if(bread->dest && bread->dest->pathID == scaffEdge->ID) break;
        			bread = bread->next;
        		}
    		}

    		if(bread){
    			breadID = bread->ID;
    			if(verbose) printf("breadID  : %i \n",breadID);
//    			if(verbose) printf("breadID  : %i \n",breadID);
				if(pog->contigNum == pog->maxNum){
					printf("ContigList full, Realloc memory\n");
					pog->contig = (struct Sequence*)realloc(pog->contig,sizeof(struct Sequence)*(pog->maxNum*2));
					pog->maxNum *= 2;
				}
				pog->contig[pog->contigNum].nsource_seq = 0;
				pog->contig[pog->contigNum].length = 0;
				sprintf(name,"Scaffold_%i_%i_%i_len:",i+1,aS->scaff[i].startJunction,aS->scaff[i].endJunction);
				pog->contig[pog->contigNum].name = (char*)malloc(strlen(name)+1);
				strcpy(pog->contig[pog->contigNum].name,name);
				if(dir){
					if(bread->sideflag) 	multidir = 3;
					else					multidir = 0;
				}
				else{
					if(bread->sideflag)		multidir = 1;
					else					multidir = 2;
				}
   				bdir = bread->sideflag;
   				decompressReadSt(reads[aS->scaff[i].startJunction].seq,readseq,reads[aS->scaff[i].startJunction].len);
   				if(multidir>1){
   					revReadSt(readseq,revreadseq);
   					strcpy(readseq,revreadseq);
   				}
//       			poa_initBackbone(&contigs[contigNum],&reads[i],readseq);
   				poa_initBackbone2(&pog->contig[pog->contigNum],&reads[aS->scaff[i].startJunction],readseq);
    			overhang = bread->overhang;
    			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
    			if(multidir%2==1 && G->read[breadID]->dir){
    				revReadSt(readseq,revreadseq);
    				strcpy(readseq,revreadseq);
    			}
    			if(multidir%2==0 && !G->read[breadID]->dir){
        			revReadSt(readseq,revreadseq);
        			strcpy(readseq,revreadseq);
    			}
    			poa_catBackbone(&pog->contig[pog->contigNum],G,&reads[breadID],readseq,startJunction,breadID);
    			poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts,overhang);
    			inserts++;
    			internb = bread;
				while(breadID != aS->scaff[i].endJunction){
					if(G->read[breadID]->flag == JUNCTION){
						printf("bread is a JUNCTION\n");
						scaffEdge = scaffEdge->next;
						if(!scaffEdge){
							printf("Catch: no more Edge, but not endJunction reached\nAbort\n");
							exit(1);
						}
						internb = G->read[breadID]->first;
						while(internb){
							if(G->read[internb->ID]->flag == CONTAINED){
			        			decompressReadSt(reads[internb->ID].seq,readseq,reads[internb->ID].len);
			        			if(multidir%2==1 && !G->read[internb->ID]->dir){
			        				revReadSt(readseq,revreadseq);
			        				strcpy(readseq,revreadseq);
			        			}
			        			if(multidir%2==0 && G->read[internb->ID]->dir){
			            			revReadSt(readseq,revreadseq);
			            			strcpy(readseq,revreadseq);
			        			}
			        			if(verbose) printf("ALINGING CONTAINED READ IN JUNCTION\n");
			        			if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose) printf("Read: %s\n",readseq);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
			        			poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[internb->ID],readseq,0,inserts,0);
	    						inserts++;
							}
							internb = internb->next;
						}
						internb = G->read[breadID]->first;
						while(internb){
							if(internb->dest && internb->dest->pathID == scaffEdge->ID){
								nextoverhang = internb->overhang;
				    			overhang += internb->overhang;
				    			oldbreadID = breadID;
				    			breadID = internb->ID;
			        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
			        			if(multidir%2==1 && G->read[breadID]->dir){
			        				revReadSt(readseq,revreadseq);
			        				strcpy(readseq,revreadseq);
			        			}
			        			if(multidir%2==0 && !G->read[breadID]->dir){
			            			revReadSt(readseq,revreadseq);
			            			strcpy(readseq,revreadseq);
			        			}
	    						poa_catBackbone2(&pog->contig[pog->contigNum],G,&reads[breadID],readseq,oldbreadID,breadID);
	//        						poa_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
	    						if(verbose) printf("ALINGING PROPER READ\n");
	    						if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose) printf("Read: %s\n",readseq);
	    						poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts,nextoverhang);
	//        						if(inserts == 290585)	poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,4077);
	//        						else poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,1);
	    						inserts++;
								break;
							}
							internb = internb->next;
						}
					}
					else{
						internb = G->read[breadID]->first;
						while(internb){
							if(G->read[internb->ID]->flag == CONTAINED){
			        			decompressReadSt(reads[internb->ID].seq,readseq,reads[internb->ID].len);
			        			if(multidir%2==1 && !G->read[internb->ID]->dir){
			        				revReadSt(readseq,revreadseq);
			        				strcpy(readseq,revreadseq);
			        			}
			        			if(multidir%2==0 && G->read[internb->ID]->dir){
			            			revReadSt(readseq,revreadseq);
			            			strcpy(readseq,revreadseq);
			        			}
			        			if(verbose) printf("ALINGING CONTAINED READ\n");
			        			if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose) printf("Read: %s\n",readseq);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
			        			poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[internb->ID],readseq,0,inserts,0);
	    						inserts++;
							}
							internb = internb->next;
						}
						internb = G->read[breadID]->first;
						while(internb){
							if(internb->sideflag == bdir && G->read[internb->ID]->flag != CONTAINED){
								nextoverhang = internb->overhang;
				    			overhang += internb->overhang;
				    			oldbreadID = breadID;
				    			breadID = internb->ID;
			        			decompressReadSt(reads[breadID].seq,readseq,reads[breadID].len);
			        			if(multidir%2==1 && G->read[breadID]->dir){
			        				revReadSt(readseq,revreadseq);
			        				strcpy(readseq,revreadseq);
			        			}
			        			if(multidir%2==0 && !G->read[breadID]->dir){
			            			revReadSt(readseq,revreadseq);
			            			strcpy(readseq,revreadseq);
			        			}
	    						poa_catBackbone2(&pog->contig[pog->contigNum],G,&reads[breadID],readseq,oldbreadID,breadID);
	//        						poa_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
	    						if(verbose) printf("ALINGING PROPER READ\n");
	    						if(verbose) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose) printf("Read: %s\n",readseq);
	    						poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts,nextoverhang);
	//        						if(inserts == 290585)	poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,4077);
	//        						else poa_heuristic_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,1);
	    						inserts++;
	//        						if(inserts == 2) return pog;
								break;
							}
							internb = internb->next;
						}
					}
//					printf("%i -> %i\n",breadID, aS->scaff[i].endJunction);
				}
    		}
    		else{
    			printf("Path for this Junction not found\n Abort!\n");
    		}
    		poa_consensus(&pog->contig[pog->contigNum]);
    		pog->contigNum++;
    		printf("Matrix time:    %.3f s\n",(float)sumMatrix/1000000000);
    		printf("Backtrace time: %.3f s\n",(float)sumTrace/1000000000);
    	}

    }

    // _________ Scaffold Consensus Calling

    free(readseq);
    free(revreadseq);
    free(name);

    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
    	free(alMatrix[i]);
    }
    free(alMatrix);

    return pog;
}


////#define SHOW_MATRIX
//
///**
// * Aligns a read sequence to the part of the poa graph, provided by the overlap information of the string graph.
// * If read is proper it is part of the backbone the new part of the PO-graph is set to this area.
// * @param contig 	Is the PO-graph of this contig
// * @param read		Is the read to align with the part of the PO-graph
// * @param seq		Is the read sequence to align
// * @param backbone	Is a boolean value if the read was proper, than it is set as new reference point for the area in the PO-graph for the next read alignment
// */
//void poa_align(struct Sequence* contig, struct read* read, char* seq, char backbone, int insNum){
//	char print_Message = 0;
//	if(print_Message) printf("CHECKPOINT: Start PO_Alignment\n");
//
//	int i,j,k;
//
//	struct Letter_T* current = &Letters[contig->readleft];
//	struct Letter_T* left;
//
////	int line = 1;
//	int line = poa_align_prepro(contig,strlen(seq),insNum);
//	int len = strlen(seq);
//
//	// build matrix along a BFS touring to the right-end node
//
//	static struct Letter_T** new_letters = NULL;
//	int new_num = 0;
//	static struct Letter_T** old_letters = NULL;
//	int old_num = 0;
//	static struct Letter_T** end_letters = NULL;
//	int end_num = 0;
//
//	if(!new_letters){
//		new_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000000); // Max breadth of graph = 100
//		old_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000000);
//		end_letters = (struct Letter_T**)malloc(sizeof(struct Letter_T*)*1000000);
//	}
//
//	struct LetterEdge* edge;
//
//	struct Letter_T* start_node = current;
//	struct Letter_T* end_node = &Letters[contig->readright];
//
//	if(print_Message) printf("Build SW-Matrix of read: %i\n",read->ID);
//#ifdef SHOW_MATRIX
//	if(insNum == 4077){
//		printf("\t");
//		for(i=0;i<len;i++){
//			printf("  %c",seq[i]);
//		}
//		printf("\n");
//	}
//
//#endif
//
//	clock_gettime(CLOCK_MONOTONIC, &ts_start);
//
//	// Set first non-gap matrix line for all nodes in the alignment-ring with the contig->readleft
//	do{
//		// TODO Think about possibility of non-having a right edge of the contig->readleft
////		if(!current->ml){
////			alMatrix_Letter[line] = current;
////			current->ml = alMatrix[line++];
////		}
//		new_letters[new_num++] = current;
//
//#ifdef SHOW_MATRIX
//		if(insNum == 4077){
//			printf("%c\t",current->letter);
//		}
//#endif
//
//		for(j=1;j<=len;j++){
//			// k is pos in seq, j-1, because j==0 is first gap position;
//			k = j-1;
//			current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(alMatrix[0][j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]),(alMatrix[0][j]+GAP_PENALTY));
//#ifdef SHOW_MATRIX
//			if(insNum == 4077){
//				if(current->ml[j]>9) printf(" %i",current->ml[j]);
//				else printf("  %i",current->ml[j]);
//			}
//#endif
//		}
//#ifdef SHOW_MATRIX
//		if(insNum == 4077){
//			printf("\n");
//		}
//#endif
//		if(current->align_ring && current->align_ring != start_node){
//			if(print_Message) printf("Found alternative start point in alignment ring:\n");
//			current = current->align_ring;
//		}
//		else{
//			break;
//		}
//	} while(1);
//
//	if(print_Message) printf("Number of alternative start positions given by an alignment ring: %i\n",new_num);
//	if(print_Message) printf("Number of nodes used in partial graph: %i\n",line);
//
//
//	// Look if the current visited node is part of a ring? and if it is, are the members of the ring part already part of the alignment? If not, they should become part of it!
//	// Answer: No don't do this, merge alignment rings: Doing during touring ;-)
//	int depth = 0;
//	int id;
//	char rightbool;
//	while(new_num && depth <= maxReadLen + 20){
////		printf("Go deeper: %i\n",depth);
//		memcpy(old_letters,new_letters,sizeof(struct Letter_T)*new_num);
//		old_num = new_num;
//		new_num = 0;
//
//		while(old_num){
//			old_num--;
//			if(old_letters[old_num] != end_node){
//				edge = old_letters[old_num]->right;
//				left = old_letters[old_num];
//				id = old_letters[old_num] - Letters;
//				if(id < 0) id *=-1;
////				printf("Letter: %c (id: %i) num: %i\n",old_letters[old_num]->letter,id,(int)old_letters[old_num]->junction);
//				if(old_letters[old_num]->junction == 1){
//					old_letters[old_num]->junction--;
//					rightbool = 1;
//				}
//				else if(old_letters[old_num]->junction > 1){
//					old_letters[old_num]->junction--;
//					rightbool = 0;
//				}
//				else{
//					printf("This case should not happen\n");
//					exit(1);
//				}
//
//				while(edge){
////					new_letters[new_num++] = &Letters[edge->dest];
//					if(rightbool){
//						new_letters[new_num++] = &Letters[edge->dest];
//					}
////					if(!Letters[edge->dest].ml){
////						alMatrix_Letter[line] = &Letters[edge->dest];
////						Letters[edge->dest].ml = alMatrix[line++];
////					}
//					current = &Letters[edge->dest];
//#ifdef SHOW_MATRIX
//					if(insNum == 4077){
//						printf("%c\t",current->letter);
//					}
//#endif
////					printf("Test 0: %c (id: %i)\n",current->letter,edge->dest);
////					if(current->ml) printf("ML is available: dpeth: %i\n",depth);
////					printf("TEst 1: %i\n",current->ml[0]);
////					printf("TEst 3: %i\n",left->ml[1]);
////					printf("TEst 4: %i\n",SM1[codes[(int)current->letter]][codes[(int)seq[1]]]);
////					printf("Test END\n");
//					for(j=1;j<=len;j++){
//						// k is pos in seq, j-1, because j==0 is first gap position;
//						k = j-1;
//						// Smith-Waterman Scoring function: Best of itself, left, diagonal, top -> Dynamic programming
//						current->ml[j] = max_func(current->ml[j],(current->ml[j-1]+GAP_PENALTY),(left->ml[j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]),(left->ml[j]+GAP_PENALTY));
//#ifdef SHOW_MATRIX
//						if(insNum == 4077){
//							if(current->ml[j]>9) printf(" %i",current->ml[j]);
//							else printf("  %i",current->ml[j]);
//						}
//#endif
//					}
//#ifdef SHOW_MATRIX
//					if(insNum == 4077){
//						printf("\n");
//					}
//#endif
//					edge = edge->next;
//				}
//			}
//			else{
//				old_letters[old_num]->junction = 0;
//				if(end_node->junction) end_node->junction = 0;
//				end_letters[end_num++] = end_node;
////				printf("EndNode found\n");
//			}
//		}
//		// Limit number of fields in matrix to compute. Give a maximum distance from the diagonal (e.g. 5bp)
//		depth++;
//	}
//
//	// Are all junctions ZERO again?
////	printf("Junctions are Zero?\n");
////	for(i=0;i<numNodes;i++){
////		if(Letters[i].junction) printf("ID: %i in num: %i\n",i,Letters[i].junction);
////	}
//
//	clock_gettime(CLOCK_MONOTONIC, &ts_finish);
//	sumMatrix += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));
//
//	if(new_num){
//		memcpy(&end_letters[end_num],new_letters,sizeof(struct Letter_T)*new_num);
//		end_num += new_num;
//	}
//
//
//	if(print_Message) printf("Number of alternative ends: %i\n",end_num);
//
//	clock_gettime(CLOCK_MONOTONIC, &ts_start);
//
//	// back tracing -> Core of this function
//	// search for highest end: could by danger if there is no penalty for gap insertions, because the score maximizes if the reference is long enough independent of the sequence similarity
//	int best_Letter = 0;
//	int best_Score = 0;
//	for(i=0;i<end_num;i++){
//		if(end_letters[i]->ml[len] > best_Score){
//			best_Letter = i;
//			best_Score = end_letters[i]->ml[len];
//		}
//	}
//
//	if(print_Message) printf("Start back tracing\n");
////	printf("Start back tracing\n");
//
//	current = end_letters[best_Letter];
//	j=len;
//	char* readseq = (char*)malloc(len*2);
//	int readlen=0;
//	char* refseq = (char*)malloc(len*2);
//	int reflen=0;
//	char leftbool = 0;
//
////	struct LetterEdge* counteredge;
//	struct Letter_T* newLetter;
//	struct Letter_T* newLetterRight = NULL;
//	int32_t newLetterRightID;
//	struct LetterEdge* newEdge;
//	struct Letter_T* current_Right;
//
//	while(1){
//		edge = current->left;
//		while(edge){
//			if(Letters[edge->dest].ml){
////				printf("edge->dest: %i (%i)\n",edge->dest,Letters[edge->dest].ml-alMatrix[0]);
//				// not sure
//				leftbool = 1;
//				left = &Letters[edge->dest];
//				k = j-1;
//
//				if(j>0 && current->ml[j] == left->ml[j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]){
//					// Entry from Diagonal
//					readseq[readlen++] = seq[k];
//					refseq[reflen++] = current->letter;
//
//					// Update POG
//					if(current->letter == seq[k]){
////						printf("Match\n");
//						// Same Letter -> Match: increase counter of existing letter
//						if(current->counter<255) current->counter++;
//						// Position Calculations depends on endianess.
//						// Check if the right letter was already a match or if it comes from a new path
//						if(newLetterRight){
//							if(current_Right == newLetterRight){
//								// Update edge counts
//								newEdge = current->right;
//								while(newEdge){
//									if(newLetterRightID == newEdge->dest){
//										if(newEdge->counter<255) newEdge->counter++;
//										break;
//									}
//									newEdge = newEdge->next;
//								}
//								newEdge = current_Right->left;
//								newLetterRightID = current-Letters;
//								if(newLetterRightID < 0) newLetterRightID *= -1;
//								while(newEdge){
//									if(newLetterRightID == newEdge->dest){
//										if(newEdge->counter<255) newEdge->counter++;
//										break;
//									}
//									newEdge = newEdge->next;
//								}
//								// last letter matched already
////								printf("\t->Last letter was a ref-match\n");
//							}
//							else{
//								// Last Letter was a mismatch and newly created
////								printf("\t->Last letter was NO ref-match\n");
//								newEdge = current->right;
//								current->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//								current->right->dest = newLetterRightID;
//								current->right->counter = 1;
//								current->right->next = newEdge;
////								printf("Letters: %p current: %p\n",Letters,current);
//								newLetterRightID = current - Letters;
//								if(newLetterRightID < 0) newLetterRightID *= -1;
////								printf("Old ID for next connection = %i\n",newLetterRightID);
//								newEdge = newLetterRight->left;
//								newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//								newLetterRight->left->dest = newLetterRightID;
//								newLetterRight->left->counter = 1;
//								newLetterRight->left->next = newEdge;
//							}
//						}
////						else{
////							printf("Diagonal and no precurser\n");
////						}
//						current_Right = current;
//						newLetterRight = current;
//						newLetterRightID = current-Letters;
//						if(newLetterRightID < 0) newLetterRightID *= -1;
////						printf("Letters: %p current: %p structSize: %i -> dif: %i\n",Letters,current,sizeof(struct Letter_T),current-Letters);
////						printf("Old ID for next connection = %i\n",newLetterRightID);
//
//					}
//					else{
//						// Differnt letters -> Mismatch: Create new Letter and new edges
////						printf("Mismatch\n");
//						newLetter = &Letters[numNodes];
//						newLetter->counter = 1;
//						newLetter->letter = seq[k];
//						newLetter->source.ipos = k;
//						newLetter->source.iseq = read->ID;
//						newLetter->source.next = NULL;
//						// Make, Close alignment ring
//						if(current->align_ring){
//							newLetter->align_ring = current->align_ring;
//							current->align_ring = newLetter;
//						}
//						else{
//							current->align_ring = newLetter;
//							newLetter->align_ring = current;
//						}
//						// Connect new letter and take the new latter as new current point -> Close path if next step is match
//						if(newLetterRight){
//							newEdge = newLetterRight->left;
//							newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//							newLetterRight->left->dest = numNodes;
//							newLetterRight->left->counter = 1;
//							newLetterRight->left->next = newEdge;
//
//							newLetter->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//							newLetter->right->dest = newLetterRightID;
////							printf("dest: %i\n",newLetter->right->dest);
//							newLetter->right->counter = 1;
//							newLetter->right->next = NULL;
//						}
//						newLetterRightID = numNodes;
//						newLetterRight = newLetter;
//						numNodes++;
//						poa_LetterSizeCheck();
//					}
//
//					current = left;
//					j--;
//					break;
//				}
//				else if(j>0 && current->ml[j] == current->ml[j-1] + GAP_PENALTY){
//					// Entry from left -> Gap in ref -> Stay in current matrix line;
//					if(print_Message) printf("Gap in REF\n");
//					readseq[readlen++] = seq[k];
//					refseq[reflen++] = '-';
//
//					// Update POG
//					newLetter = &Letters[numNodes];
//					newLetter->counter = 1;
//					newLetter->letter = seq[k];
//					newLetter->source.ipos = k;
//					newLetter->source.iseq = read->ID;
//					newLetter->source.next = NULL;
//					if(newLetterRight){
//						newEdge = newLetterRight->left;
//						newLetterRight->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//						newLetterRight->left->dest = numNodes;
//						newLetterRight->left->counter = 1;
//						newLetterRight->left->next = newEdge;
//
//						newLetter->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
//						newLetter->right->dest = newLetterRightID;
//						newLetter->right->counter = 1;
//						newLetter->right->next = NULL;
//					}
//					newLetterRightID = numNodes;
//					newLetterRight = newLetter;
//					numNodes++;
//					poa_LetterSizeCheck();
//
//					j--;
//					continue;
//				}
//				else if(current->ml[j] == left->ml[j] + GAP_PENALTY){
//					// Entry from above -> Gap in seq
//					if(print_Message) printf("Gap in Seq\n");
//					readseq[readlen++] = '-';
//					refseq[reflen++] = current->letter;
//					current_Right = NULL;
//					current = left;
//					break;
//				}
//				else{
////					printf("Noting is true, there must be an alternative path!\n");
////					printf("j: %i",j);
//				}
//			}
//			edge = edge->next;
//		}
//		if(!leftbool) break;
//		else leftbool = 0;
//	}
//	// Connect to matrix origin
//	k = j-1;
//	if(print_Message){
//		if(j>0) printf("Connect to matrix origin j (score: %i) (%c/%c): %i\n",current->ml[j],current->letter,seq[k],j);
//	}
//	while(j!=0){
//		k = j-1;
//		if(print_Message){
//			if(j>0) printf("%i -> %i = %i + %i\n",current->ml - alMatrix[0],current->ml[j], alMatrix[0][j-1],SM1[codes[(int)current->letter]][codes[(int)seq[k]]]);
//		}
//		if(j>0 && current->ml[j] == alMatrix[0][j-1] + SM1[codes[(int)current->letter]][codes[(int)seq[k]]]){
//			// Entry from Diagonal
//			if(print_Message) printf("Origin Diagonal\n");
//			readseq[readlen++] = seq[k];
//			refseq[reflen++] = current->letter;
//			if(current->letter == seq[k]){
//				current->counter++;
//			}
//			else{
//				// make new letter
//				// how to connect???
//			}
//			j--;
//			for(;j>0;j--){
//				readseq[readlen++] = seq[k];
//				refseq[reflen++] = '-';
//				if(print_Message) printf("go left till origin is reached\n");
//			}
//			break;
//		}
//		else if(j>0 && current->ml[j] == current->ml[j-1] + GAP_PENALTY){
//			// Entry from left -> Gap in ref -> Stay in current matrix line;
//			if(print_Message) printf("Origin left\n");
//			readseq[readlen++] = seq[k];
//			refseq[reflen++] = '-';
//			j--;
//			continue;
//		}
//		else if(current->ml[j] == alMatrix[0][j] + GAP_PENALTY){
//			// Entry from above -> Gap in seq
//			if(print_Message) printf("Origin top\n");
//			readseq[readlen++] = '-';
//			refseq[reflen++] = current->letter;
//			for(;j>0;j--){
//				readseq[readlen++] = seq[k];
//				refseq[reflen++] = '-';
//				if(print_Message) printf("go left till origin is reached\n");
//			}
//			break;
//		}
//		else{
//			printf("%c (j: %i)\n",seq[k],j);
//			printf("-\t");
//			for(i=0;i<=len;i++){
//				printf(" %i",alMatrix[0][i]);
//			}
//			printf("\n");
//			printf("%c\t",current->letter);
//			for(i=0;i<=len;i++){
//				printf(" %i",current->ml[i]);
//			}
//			printf("\n");
//			printf("Number of ends: %i\n",end_num);
//			printf("Depth: %i\n", depth);
//			printf("Matrixline: %i\n",alMatrix[0]-current->ml);
//			printf("Nothing is true: This Case should not happen\n");
//			poa_showAlignment(readseq,refseq,readlen);
//			exit(1);
//
//		}
//	}
//	if(print_Message) printf("j: %i\n",j);
//
//	clock_gettime(CLOCK_MONOTONIC, &ts_finish);
//	sumTrace += (((ts_finish.tv_sec * 1000000000) + ts_finish.tv_nsec) - ((ts_start.tv_sec * 1000000000) + ts_start.tv_nsec));
//
//
//	char print_align = 0;
//
//	// Show Alignment
//
//
//	if(print_align){
//		poa_showAlignment(readseq,refseq,readlen);
//	}
//
//	// free the letters and set matrix to 0
//	for(i=1;i<line;i++){
//		alMatrix_Letter[i]->ml = NULL;
//		for(j=1;j<=len;j++){
//			alMatrix[i][j] = (i+j) * GAP_PENALTY;
//		}
//	}
//
//	free(readseq);
//	free(refseq);
//
////	exit(1);
//}