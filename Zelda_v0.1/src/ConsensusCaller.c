/*
 ============================================================================
 Name        : ConsensusCaller.c
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
#include <time.h>
#include <string.h>
#include "ConsensusCaller.h"
#include "DBGraph_scaffold.h"
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
#define MIN_SCAFF_LEN 200
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
uint32_t          maxNumNodes =   1000000;	// Initial max size of LPOLetter_T array

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
//				Letters[j].source.next = NULL;
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
void poa_initBackbone2(struct Sequence* contig, char* seq){ //struct reads* read
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
//		current->source.ipos = i;
//		current->source.iseq = read->ID;
//		current->source.next = NULL;
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
			left->right->vFlag = 0;
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

void poa_catBackbone(struct Sequence* contig, struct myovlList *G, char* seq, int leftID, int rightID){ // Parameter :  struct reads* read,
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
	struct LetterEdge* oldEdge;

//	printf("LeftID: %i\n",leftID);

	for(;i<len;i++){
		current = &Letters[numNodes];
		current->letter = seq[i];
		current->align_ring = NULL;
		current->ml = NULL;
		current->counter=1;
//		// _________
//		// This things not until alignment
//		current->source.ipos = i;
//		current->source.iseq = read->ID;
//		current->source.next = NULL;
//		// _________
//		printf("NewLetter: %i\n",numNodes);
		current->left = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		current->left->counter = 1;
		current->left->dest = leftID;
		current->left->next = NULL;
//		printf("Set Right edge to leftID: %i\n",leftID);
		left = &Letters[leftID];
		oldEdge = left->right;
		left->right = (struct LetterEdge*)malloc(sizeof(struct LetterEdge));
		left->right->counter = 1;
		left->right->dest = numNodes;
		left->right->vFlag = 0;
		left->right->next = oldEdge;

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
	printf("Write POA-Graph to : %s\n",dotFile);
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

void poa_avgCov(struct Sequence* contig){
	int i;
	uint32_t totCov = 0;
	for(i=0;i<numNodes;i++){
		totCov += Letters[i].counter;
	}
	contig->avgCov = (float)totCov/contig->length;
	printf("Contig: %s \t\t AvgCov: %0.2f\n",contig->name,contig->avgCov);
}

static inline void resetLetterNum(int* letters){
	letters[0] = 0;
	letters[1] = 0;
	letters[2] = 0;
	letters[3] = 0;
}

static inline void resetLetterSt(struct LetterEdge** letters){
	letters[0] = NULL;
	letters[1] = NULL;
	letters[2] = NULL;
	letters[3] = NULL;
}

long double variant_qualaty(int n,int k){
	long double right = 0.995;
	long double wrong = 0.005;
    long long ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return (-10 * log10((ans*(pow(wrong,k)*pow(right,n-k)))));
}

#define maxAltLen 100

static inline unsigned char poa_makeCigar(char* cigar, char* ref, char* alt){
	unsigned char numM = 0;
	unsigned char numX = 0;
	unsigned char totX = 0;
	uint16_t rl = strlen(ref);
	uint16_t al = strlen(alt);
	cigar[0]='\0';
	// SNP
	if(al==1 && rl ==1){
		sprintf(cigar,"1X");
		return 0;
	}
	else{

	}
	int min = _min(rl,al);
	int i;
	for(i=0;i<min;i++){
		if(ref[i] == alt[i]){
			if(numM) numM++;
			else{
				if(numX){
					sprintf(cigar,"%s%iX",cigar,(int)numX);
					numX=0;
				}
				numM++;
			}
		}
		else{
			if(numX) numX++;
			else{
				if(numM){
					sprintf(cigar,"%s%iM",cigar,(int)numM);
					numM=0;
				}
				numX++;
			}
			totX++;
		}
	}
	if(numM) sprintf(cigar,"%s%iM",cigar,(int)numM);
	if(numX) sprintf(cigar,"%s%iX",cigar,(int)numX);
	if(rl == al) return 1;
	else{
		if(al > rl){
			sprintf(cigar,"%s%iI",cigar,al-rl);
			if(totX) return 4;
			else return 2;
		}
		else{
			sprintf(cigar,"%s%iD",cigar,rl-al);
			if(totX) return 4;
			else return 3;
		}
	}
	return 0;
}

void poa_deleteVariant(struct POG* pog){
	int i;
	struct Variation* var;
	for(i=0;i<pog->contigNum;i++){
		var = pog->contig[i].var;
		while(var){
			pog->contig[i].var = var->next;
			free(var->altSeq);
			free(var->refSeq);
			free(var);
			var = pog->contig[i].var;
		}
	}
}

void poa_reportVariant(struct POG* pog, char* vcfFile, char* ref){
	// Write Alternative to VCF
	FILE* vcf = fopen(vcfFile,"w");

	// ToDo: Write Header
    time_t current_time;
    char* cigar = (char*)malloc(2*maxAltLen);
    char* c_time_string;
    /* Obtain current time. */
    current_time = time(NULL);

    if (current_time == ((time_t)-1))
    {
        (void) fprintf(stderr, "Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL)
    {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }
    c_time_string[strlen(c_time_string)-1] = '\0';
    fprintf(vcf,"##fileformat=VCFv4.1\n");
    fprintf(vcf,"##fileDate=\"%s\"\n",c_time_string);
    fprintf(vcf,"##source=%s\n",version);
    fprintf(vcf,"##reference=%s\n",ref);
    fprintf(vcf,"##phasing=none\n");
    fprintf(vcf,"##filter=\"QUAL > 20\"\n");
    fprintf(vcf,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    fprintf(vcf,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n");
//    fprintf(vcf,"##INFO=<ID=DPB,Number=1,Type=Float,Description=\"Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype\">\n");
//    fprintf(vcf,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n");
//    fprintf(vcf,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
//    fprintf(vcf,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">\n");
    fprintf(vcf,"##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count, with partial observations recorded fractionally\">\n");
    fprintf(vcf,"##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=PRO,Number=1,Type=Float,Description=\"Reference allele observation count, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=PAO,Number=A,Type=Float,Description=\"Alternate allele observations, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=QR,Number=1,Type=Integer,Description=\"Reference allele quality sum in phred\">\n");
//    fprintf(vcf,"##INFO=<ID=QA,Number=A,Type=Integer,Description=\"Alternate allele quality sum in phred\">\n");
//    fprintf(vcf,"##INFO=<ID=PQR,Number=1,Type=Float,Description=\"Reference allele quality sum in phred for partial observations\">\n");
//    fprintf(vcf,"##INFO=<ID=PQA,Number=A,Type=Float,Description=\"Alternate allele quality sum in phred for partial observations\">\n");
//    fprintf(vcf,"##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">\
//    fprintf(vcf,"##INFO=<ID=SAF,Number=A,Type=Integer,Description=\"Number of alternate observations on the forward strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SAR,Number=A,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SRP,Number=1,Type=Float,Description=\"Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=SAP,Number=A,Type=Float,Description=\"Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=AB,Number=A,Type=Float,Description=\"Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">\n");
//    fprintf(vcf,"##INFO=<ID=ABP,Number=A,Type=Float,Description=\"Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RUN,Number=A,Type=Integer,Description=\"Run length: the number of consecutive repeats of the alternate allele in the reference genome\">\n");
//    fprintf(vcf,"##INFO=<ID=RPP,Number=A,Type=Float,Description=\"Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RPPR,Number=1,Type=Float,Description=\"Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RPL,Number=A,Type=Float,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">\n");
//    fprintf(vcf,"##INFO=<ID=RPR,Number=A,Type=Float,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">\n");
//    fprintf(vcf,"##INFO=<ID=EPP,Number=A,Type=Float,Description=\"End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=EPPR,Number=1,Type=Float,Description=\"End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=DPRA,Number=A,Type=Float,Description=\"Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.\">\n");
//    fprintf(vcf,"##INFO=<ID=ODDS,Number=1,Type=Float,Description=\"The log odds ratio of the best genotype combination to the second-best.\">\n");
//    fprintf(vcf,"##INFO=<ID=GTI,Number=1,Type=Integer,Description=\"Number of genotyping iterations required to reach convergence or bailout.\">\n");
    fprintf(vcf,"##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n");
    fprintf(vcf,"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.  Note that INDEL alleles do not have the first matched base (which is provided by default, per the spec) referred to by the CIGAR.\">\n");
//    fprintf(vcf,"##INFO=<ID=NUMALT,Number=1,Type=Integer,Description=\"Number of unique non-reference alleles in called genotypes at this position.\">\n");
//    fprintf(vcf,"##INFO=<ID=MEANALT,Number=A,Type=Float,Description=\"Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.\">\n");
    fprintf(vcf,"##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">\n");
//    fprintf(vcf,"##INFO=<ID=MQM,Number=A,Type=Float,Description=\"Mean mapping quality of observed alternate alleles\">\n");
//    fprintf(vcf,"##INFO=<ID=MQMR,Number=1,Type=Float,Description=\"Mean mapping quality of observed reference alleles\">\n");
//    fprintf(vcf,"##INFO=<ID=PAIRED,Number=A,Type=Float,Description=\"Proportion of observed alternate alleles which are supported by properly paired read fragments\">\n");
//    fprintf(vcf,"##INFO=<ID=PAIREDR,Number=1,Type=Float,Description=\"Proportion of observed reference alleles which are supported by properly paired read fragments\">\n");
    fprintf(vcf,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
//    fprintf(vcf,"##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">\n");
//    fprintf(vcf,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n");
    fprintf(vcf,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    fprintf(vcf,"##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n");
    fprintf(vcf,"##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">\n");
    fprintf(vcf,"##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n");
    fprintf(vcf,"##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">\n");
    fprintf(vcf,"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  unknown\n");

	int i;
	long double qual;
	struct Variation* var;
	for(i=0;i<pog->contigNum;i++){
		var = pog->contig[i].var;
		while(var){
			qual = variant_qualaty(var->dp,var->ao);
			if(qual>20){
				var->type = poa_makeCigar(cigar,var->refSeq,var->altSeq);
				fprintf(vcf,"%s\t",pog->contig[i].name);
				fprintf(vcf,"%i\t",var->pos);
				fprintf(vcf,".\t");
				fprintf(vcf,"%s\t",var->refSeq);
				fprintf(vcf,"%s\t",var->altSeq);
				fprintf(vcf,"%.3Lf\t",qual);
				fprintf(vcf,".\t");
				fprintf(vcf,"AO=%i;CIGAR=%s;DP=%i;LEN=%i;RO=%i;TYPE=%s\t",var->ao,cigar,var->dp,var->len,var->ro,varType[(int)var->type]);
				fprintf(vcf,"GT:DP:RO:QR:AO:QA\t");
				fprintf(vcf,"%i:%i:%i:%i:%i:%i\t",1,var->dp,var->ro,1,var->ao,1);
				fprintf(vcf,"\n");
			}
			var = var->next;
		}
	}

	free(cigar);

	fclose(vcf);
}

void poa_recMainPath(struct Letter_T* currentLetter, struct Letter_T* endLetter, int startPos, int altLen, char* altSeq, int altCov, struct Sequence* contig){
	char verbose = 0;
	if(verbose) printf("Checkpoint recMainPath");
	char* refSeq = (char*)malloc(maxAltLen);
	int refLen = 0;
	int refCov = 100000;
	struct LetterEdge* edge;
	struct Letter_T* startLetter = currentLetter;
	while(currentLetter != endLetter){
		if(refLen == maxAltLen){
			printf("No Proper Ref Path Found (Length Error)\n");
			return;
		}
		refSeq[refLen++] = currentLetter->letter;
		edge = currentLetter->right;
		if(edge){
			while(edge){
				if(edge->vFlag){
					refCov = _min(refCov,edge->counter);
					break;
				}
				edge = edge->next;
			}
			currentLetter = & Letters[edge->dest];
		}
	}
	char type = -1;
	if(altLen == 2 && refLen == 2) type = 0; 	// snp
	else if(altLen == refLen) type = 1; 		// mnp
	else if(altLen > refLen) type = 2;			// ins
	else if(altLen < refLen) type = 3;			// del

	struct Variation* var = (struct Variation*)malloc(sizeof(struct Variation));
	if(type == 0){
		var->altSeq = (char*)malloc(2);
		var->altSeq[0] = altSeq[1];
		var->altSeq[1] = '\0';
		var->refSeq = (char*)malloc(2);
		var->refSeq[0] = refSeq[1];
		var->refSeq[1] = '\0';
	}
	else{
		altSeq[altLen] = '\0';
		refSeq[refLen] = '\0';
		var->altSeq = (char*)malloc(strlen(altSeq)+1);
		var->refSeq = (char*)malloc(strlen(refSeq)+1);
		strcpy(var->altSeq,altSeq);
		strcpy(var->refSeq,refSeq);
	}
	var->dp = _min(startLetter->counter,endLetter->counter);
	var->ao = altCov;
	var->ro = refCov;
	var->pos = startPos;
	if(type == 0) var->pos++;
	var->next = NULL;
	var->type = type;
	// ToDo: rest entries down under
	var->len = strlen(var->altSeq);

	if(!contig->var){
		contig->var = var;
	}
	else contig->lastvar->next = var;
	contig->lastvar = var;
	if(verbose) printf("Write Variation:\n");
	if(verbose) printf("\t%s\t%i\t%s\t%s\tDP:%i;AO:%i;RO:%i\n",contig->name,var->pos,var->refSeq,var->altSeq,var->dp,var->ao,var->ro);

	free(refSeq);
	// walk main path to Letter is End of Alt Path
	// call poa_reportVariant()
}

/**
 * DFS-Search
 */
void poa_recVariantPath(struct Letter_T* startLetter, int startPos, int len, char* seq, struct LetterEdge* edge, int cov, struct Sequence* contig){
	// Include Path counting
	struct Letter_T* current = &Letters[edge->dest];
	// IF current is Consensus Path call poa_recMainPath()
	if(current->vFlag) poa_recMainPath(startLetter,current,startPos,len,seq,cov,contig);
	// ELSE next not flagged go deeper
	else{
		edge = current->right;
		if(edge){
			seq[len] = current->letter;
			while(edge){
				if(edge->counter > 2 && len < maxAltLen){
					poa_recVariantPath(startLetter,startPos,len+1,seq,edge,_min(cov,edge->counter),contig);
				}
				edge = edge->next;
			}
		}
	}
}

void poa_variantCalling(struct Sequence* contig){
//	int verbose = 0;
	int i=0;
	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	char* altPath = (char*)malloc(2*maxAltLen);
	contig->var = NULL;
	contig->lastvar = NULL;

	while(1){
		edge = current->right;
		if(edge){
			while(edge){
				if(!edge->vFlag && edge->counter > 2){
					altPath[0] = current->letter;
					poa_recVariantPath(current,i,1,altPath,edge,edge->counter,contig);
				}
				edge = edge->next;
			}
			edge = current->right;
			while(edge){
				if(edge->vFlag) break;
				edge = edge->next;
			}
			current = &Letters[edge->dest];
			i++;
		}
		else{
			break;
		}
	}
	free(altPath);

}


void poa_consensus2(struct Sequence* contig){
	printf("CHECKPOINT: PO-MSA to Contig\n");
	int verbose = 0;
//	if(strcmp(contig->name,"Scaffold_19_90086_186097_len:")==0) verbose = 1;
	char* seq = (char*)malloc(contig->length + 100);
	int i=0;

	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	struct LetterEdge* bestedge;
	struct Letter_T* ring;
	struct Letter_T* bestRing = NULL;

	printf("Contig: %s%i\n",contig->name,contig->length);

	struct LetterEdge* lettersSt[4];

	while(1){
		if(current->counter < 5) seq[i++] = current->letter+32;
		else seq[i++] = current->letter;
		if(current->vFlag) printf("Flag was set\n");
		current->vFlag = 1;
		edge = current->right;
		resetLetterSt(lettersSt);
		if(edge){
			bestedge = edge;
			while(edge){
				lettersSt[codes[(int)Letters[edge->dest].letter]] = edge;
				if(bestedge->counter < edge->counter){
					bestedge = edge;
				}
				edge = edge->next;
			}
			edge = bestedge;
			bestRing = &Letters[bestedge->dest];
			if(Letters[bestedge->dest].align_ring){
				ring = Letters[bestedge->dest].align_ring;
				if(verbose) printf("Pos: %i -> BestRing Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,bestRing->letter,(int)bestRing->letter,bestRing->counter);
				while(ring != &Letters[bestedge->dest]){
					if(verbose) printf("\tPos: %i -> Ring Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,ring->letter,codes[(int)ring->letter],ring->counter);
					if(lettersSt[codes[(int)ring->letter]] && ring->counter > bestRing->counter && ring == &Letters[lettersSt[codes[(int)ring->letter]]->dest]){
						bestRing = ring;
						edge = lettersSt[codes[(int)ring->letter]];
					}
					ring = ring->align_ring;
				}
				if(verbose) printf("Pos: %i -> BestRing Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,bestRing->letter,(int)bestRing->letter,bestRing->counter);
			}
			edge->vFlag = 1;
			current = bestRing;
		}
		else{
			current->vFlag = 1;
			break;
		}
		bestRing = NULL;
	}

	seq[i]='\0';
//	printf(">Correct_%s\n",contig->name);
//	int k;
//	for(k=0;k<i;k+=80){
//		printf("%.80s\n",&seq[k]);
//	}
//	printf("\n");
	contig->sequence = seq;
	contig->length = strlen(seq);
	sprintf(contig->name,"%s%i",contig->name,contig->length);
	poa_variantCalling(contig);
	poa_avgCov(contig);
	sprintf(contig->name,"%s_avgCov:%.2f",contig->name,contig->avgCov);
	if(verbose) exit(1);
}

void poa_consensus(struct Sequence* contig){
	printf("CHECKPOINT: PO-MSA to Contig\n");
	int verbose = 0;
	if(strcmp(contig->name,"Scaffold_19_90086_186097_len:")==0) verbose = 1;
	char* seq = (char*)malloc(contig->length + 100);
	int i=0;
	int j;

	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	struct LetterEdge* bestedge;
	struct Letter_T* ring;
	struct Letter_T* bestRing = NULL;

	printf("Contig: %s%i\n",contig->name,contig->length);

	int letters[4];
	int totLetters;

	while(1){
		if(current->counter < 5) seq[i++] = current->letter+32;
		else seq[i++] = current->letter;
		if(current->align_ring){
			resetLetterNum(letters);
			letters[codes[(int)current->letter]] = current->counter;
			totLetters = current->counter;
			ring = current->align_ring;
			bestRing = current;
			while(ring != current){
				if(verbose) printf("Pos: %i -> Ring Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,ring->letter,(int)ring->letter,ring->counter);
				letters[codes[(int)ring->letter]] = ring->counter;
				if(ring->counter > bestRing->counter) bestRing = ring;
				totLetters += ring->counter;
				ring = ring->align_ring;
			}
			if(current->counter < totLetters - 3){
				printf("Variation at Pos: %i (Bases: %i)\n",i-1,totLetters);
				printf("\tBestLetter: %c (%i)\n",current->letter,current->counter);
				printf("BestRingLetter: %c (%i)\n",bestRing->letter,bestRing->counter);
				for(j=0;j<4;j++){
					printf("\t\t%c: %i\n",rev_codes[j],letters[j]);
				}
			}
		}
		edge = current->right;
		if(edge){
			bestedge = edge;
			while(edge){
				if(bestRing && &Letters[edge->dest] == bestRing){
					bestedge = edge;
					break;
				}
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
		bestRing = NULL;
	}

	seq[i]='\0';
//	printf(">Correct_%s\n",contig->name);
//	int k;
//	for(k=0;k<i;k+=80){
//		printf("%.80s\n",&seq[k]);
//	}
//	printf("\n");
	contig->sequence = seq;
	if(verbose) exit(1);
}

void poa_printContigs(struct POG* pog, char* contigFile){
	printf("CHECKPOINT: Write CorrectContigs in fasta\n");
	FILE* correctContigs = fopen(contigFile,"w");
	char hideEnds = 0;

	int i,j,k;
	int len;
	int lenNew;

	int nextID;
	int insert;
	for(i=0;i<pog->contigNum;i++){
		if(!pog->contig[i].vflag){
			len = strlen(pog->contig[i].sequence);
			fprintf(correctContigs,">%s\n",pog->contig[i].name);
			j=0;
			if(hideEnds){
				if(len>480){
					j = 160;
					len -= 160;
				}
			}
			for(;j<len;j+=80){
				if(j+80 < len) fprintf(correctContigs,"%.80s\n",&pog->contig[i].sequence[j]);
				else fprintf(correctContigs,"%.80s",&pog->contig[i].sequence[j]);
			}
			nextID = i;
			while(pog->contig[nextID].seqEdge){
				printf("Write Fasta Over Bridge\n");
				j = len;
				insert = pog->contig[nextID].seqEdge->insertLen;
				if(insert<10) insert = 10;
				len += insert;
				nextID = pog->contig[nextID].seqEdge->nextScaff;
				while(j<len){
					j++;
					fprintf(correctContigs,"N");
					if(j%80==0) fprintf(correctContigs,"\n");
				}
				printf("Write Fasta Behind Bridge\n");
				lenNew = strlen(pog->contig[nextID].sequence);
				printf("Write Fasta Behind Bridge (len: %i)\n",lenNew);
				len += lenNew;
				k=0;
				if(insert < 0) k = insert*-1;
				for(;k<lenNew;k+=80){
					if(k+80 < len) fprintf(correctContigs,"%.80s\n",&pog->contig[nextID].sequence[k]);
					else fprintf(correctContigs,"%.80s",&pog->contig[nextID].sequence[k]);
				}
			}
			fprintf(correctContigs,"\n");
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
//    	Letters[i].source.next = NULL;
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
    				pog->contig[pog->contigNum].name = (char*)malloc(strlen(name)+10);
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
       				poa_initBackbone2(&pog->contig[pog->contigNum],readseq); // parameter 2 ,&reads[i]
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
        			poa_catBackbone(&pog->contig[pog->contigNum],G,readseq,i,breadID); // Parameter 3: &reads[breadID],
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
            						poa_catBackbone2(&pog->contig[pog->contigNum],G,readseq,oldbreadID,breadID); // Parameter 3 &reads[breadID]
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
	char verbose = 1;
    int gesLen = 0;													// Sum over all scaffold length
	int *nStat = (int*)malloc(sizeof(int)*aS->num);				// List of Scaffold length

	int i;

    int v;
    int k=0;
    int anzlen;
    for (v = 0; v < aS->num; v++){
    	if(aS->scaff[v].len >= 200) nStat[k++] = aS->scaff[v].len;
    }

	int temp;
	anzlen = k;
//	k = aS->num;
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

	for(i=0;i<anzlen;i++) gesLen += nStat[i];

	int sum=0;
	int ns=0;

	struct scaffEdge* scaffEdge;
	int startJunction;

	if(verbose){
		printf("Scaffold Paths -> num: %i\n",anzlen);
		printf("ScaffStat:\n");
		aS->numbridge = aS->num;
		int bridgeJunction;
		for(i=0;i<aS->numbridge;i++){
			if(aS->scaff[i].len < 200) continue;
			startJunction = aS->scaff[i].startJunction;
			printf("Scaffold: %i (len: %i bp) Type: %i\n",i,aS->scaff[i].len,aS->scaff[i].type);
			scaffEdge = aS->scaff[i].first;
			printf(KRED"%i"KNRM,startJunction);
			while(scaffEdge){
				if(scaffEdge->bridge){
					if(paths[scaffEdge->ID].leftJunction == scaffEdge->targetJunction) bridgeJunction = paths[scaffEdge->ID].rightJunction;
					else bridgeJunction = paths[scaffEdge->ID].leftJunction;
					printf(" ..(%i)..> "KYEL"%i"KNRM,scaffEdge->bridge->estLen,bridgeJunction);
					printf(" -> "KGRN"%i"KNRM,scaffEdge->ID);
					printf(" -> "KRED"%i"KNRM,scaffEdge->targetJunction);
				}
				else{
					printf(" -> "KGRN"%i"KNRM,scaffEdge->ID);
					printf(" -> "KRED"%i"KNRM,scaffEdge->targetJunction);
				}
				scaffEdge = scaffEdge->next;
			}
			printf("\n");
		}
	}

//	printf("ScaffStat 3 !!!	\n");
	int scaffID;
	int len;
	int bridgeJunction;
	struct scaffEdge* oldscaffEdge = NULL;
	for(i=0;i<aS->numbridge;i++){
		aS->scaff[i].next = NULL;
		if(aS->scaff[i].len < 200 && i<aS->num) continue;
		startJunction = aS->scaff[i].startJunction;
//		printf("Scaffold: %i (len: %i bp) Type: %i\n",i,aS->scaff[i].len,aS->scaff[i].type);
		scaffEdge = aS->scaff[i].first;
		len = scaffEdge->len;
//		printf(KRED"%i"KNRM,startJunction);
		scaffID = i;
		while(scaffEdge){
			if(scaffEdge != aS->scaff[i].first && scaffEdge->bridge){
				if(paths[scaffEdge->ID].leftJunction == scaffEdge->targetJunction) bridgeJunction = paths[scaffEdge->ID].rightJunction;
				else bridgeJunction = paths[scaffEdge->ID].leftJunction;
				// New ScaffEdge
				aS->scaff[aS->numbridge].ID = aS->num;
				aS->scaff[aS->numbridge].startJunction = bridgeJunction;
				aS->scaff[aS->numbridge].endJunction = scaffEdge->targetJunction;
				aS->scaff[aS->numbridge].len = aS->scaff[i].len - (len - scaffEdge->len);
				aS->scaff[aS->numbridge].type = aS->scaff[i].type;
				aS->scaff[aS->numbridge].first = scaffEdge;
				// Update old scaffEdge
				aS->scaff[i].len = len - scaffEdge->len;
				aS->scaff[i].next = &aS->scaff[aS->numbridge];
				oldscaffEdge->next = NULL;
//				printf(" ..(%i)..> "KYEL"%i"KNRM,scaffEdge->bridge->estLen,bridgeJunction);
//				printf(" -> "KGRN"%i"KNRM,scaffEdge->ID);
//				printf(" -> "KRED"%i"KNRM,scaffEdge->targetJunction);
				aS->numbridge++;
				break;
			}
			else{
//				printf(" -> "KGRN"%i"KNRM,scaffEdge->ID);
//				printf(" -> "KRED"%i"KNRM,scaffEdge->targetJunction);
			}
			if(scaffEdge->next && scaffEdge->next->bridge){
				oldscaffEdge = scaffEdge;
			}
			scaffEdge = scaffEdge->next;
			len += scaffEdge->len;
		}
//		printf("\n");
	}

	printf("\n");
	printf("Largest Scaffold: %i bp\n",nStat[0]);
	printf("Number of Scaffolds (>=200bp): %i\n",anzlen);
	printf("Total Length over all Scaffolds: %i\n",gesLen);
	for(i=0;i<anzlen;i++){
		sum += nStat[i];
		if(sum > (gesLen/10) && ns == 0){
			printf("N10: %i (L10: %i)\n",nStat[i],i+1);
			ns++;
		}
		if(sum > (gesLen/4) && ns == 1){
			printf("N25: %i (L25: %i)\n",nStat[i],i+1);
			ns++;
		}
		if(sum > (gesLen/2) && ns == 2){
			printf("N50: %i (L50: %i)\n",nStat[i],i+1);
			ns++;
		}
		if(sum > (gesLen/4)*3 && ns == 3){
			printf("N75: %i (L75: %i)\n",nStat[i],i+1);
			ns++;
		}
		if(sum > (gesLen/10)*9 && ns == 4){
			printf("N90: %i (L90: %i)\n\n",nStat[i],i+1);
			ns++;
		}
	}
	free(nStat);
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

struct scaffold_set* scaffold_init_old(){
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
	return aS;
}

struct scaffold_set* scaffold_init(){
    int i;
    struct scaffold_set* aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
    aS->num = 0;
    aS->nummax = 1000;
    aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*aS->nummax);

    char verbose = 0;

    int depth = 0;
    int len = 0;

    int lpos = 0;
    int rpos = 0;
    int lelem = 0;
    int relem = 0;
    int lmaxelem = 1000;
    int rmaxelem = 1000;
    int current;
    struct contigScaff* left = (struct contigScaff*)malloc(sizeof(struct contigScaff)*lmaxelem);
    struct contigScaff* right = (struct contigScaff*)malloc(sizeof(struct contigScaff)*rmaxelem);
    struct pathEdge* edge = (struct pathEdge*)malloc(sizeof(struct pathEdge));

    for(i=1;i<pathsNum;i++){
    	if(!paths[i].flag){
    		if((paths[i].scaffflag & 32) && (paths[i].scaffflag & 2)){
    			continue;
    		}
    		if(verbose) printf("\t NEW SCAFFOLD\n");
    		paths[i].flag++;
    		// initial left
    		lpos = 0;
    		rpos = 0;
    		lelem = 0;
    		relem = 0;
    		edge = paths[i].leftPath;
    		while(edge && !edge->sibl){
    			if(verbose) printf("(%i) Set left %i\n",i,edge->ID);
    			left[lelem].ID = edge->ID;
    			if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
    			else left[lelem].sameside = 0;
    			lelem++;
    			if(lelem == lmaxelem){
    				lmaxelem *= 2;
    				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
    				if(!left){
    					printf("Error in realloc lmaxelem in scaffold_init\n");
    					exit(1);
    				}
    			}
    			edge = edge->next;
    		}
    		// initial right
    		edge = paths[i].rightPath;
    		while(edge && !edge->sibl){
    			if(verbose) printf("(%i) Set right %i\n",i,edge->ID);
    			right[relem].ID = edge->ID;
    			if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
    			else right[relem].sameside = 0;
    			relem++;
    			if(relem == rmaxelem){
    				rmaxelem *= 2;
    				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
    				if(!right){
    					printf("Error in realloc rmaxelem in scaffold_init\n");
    					exit(1);
    				}
    			}
    			edge = edge->next;
    		}
    		while(lpos<lelem || rpos<relem){
    			if(verbose) printf("Somithing was set, go deeper\n");
    			// left
    			while(lpos<lelem){
    				current = left[lpos].ID;
    				paths[current].flag++;
    				// left -> left
    				if(left[lpos].sameside) edge = paths[current].leftPath;
    				else edge = paths[current].rightPath;
    				while(edge && !edge->sibl){
    					if(edge->depth + lpos == lelem){
    						if(verbose) printf("(%i) (lpos: %i) Set left left %i (target: %i)\n",i,lpos,edge->ID,edge->targetJunction);
    						left[lelem].ID = edge->ID;
    		    			if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
    		    			else left[lelem].sameside = 0;
    		    			lelem++;
    		    			if(lelem == lmaxelem){
    		    				lmaxelem *= 2;
    		    				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
    		    				if(!left){
    		    					printf("Error in realloc lmaxelem in scaffold_init\n");
    		    					exit(1);
    		    				}
    		    			}
    					}
    					edge = edge->next;
    				}
    				// left -> right
    				if(left[lpos].sameside) edge = paths[current].rightPath;
    				else edge = paths[current].leftPath;
    				while(edge && !edge->sibl){
    					if(edge->depth == lpos + relem +2){
    						if(verbose) printf("(%i) (lpos: %i) Set left right %i",i,lpos,edge->ID);
    						right[relem].ID = edge->ID;
    		    			if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
    		    			else right[relem].sameside = 0;
    		    			relem++;
    		    			if(relem == rmaxelem){
    		    				rmaxelem *= 2;
    		    				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
    		    				if(!right){
    		    					printf("Error in realloc rmaxelem in scaffold_init\n");
    		    					exit(1);
    		    				}
    		    			}
    					}
    					edge = edge->next;
    				}
    				lpos++;
    			}
    			// right
    			while(rpos < relem){
    				current = right[rpos].ID;
    				paths[current].flag++;
    				// right -> right
    				if(right[rpos].sameside) edge = paths[current].rightPath;
    				else edge = paths[current].leftPath;
    				while(edge && !edge->sibl){
    					if(edge->depth + rpos == relem){
    						right[relem].ID = edge->ID;
    						if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
    						else right[relem].sameside = 0;
    						relem++;
    		    			if(relem == rmaxelem){
    		    				rmaxelem *= 2;
    		    				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
    		    				if(!right){
    		    					printf("Error in realloc rmaxelem in scaffold_init\n");
    		    					exit(1);
    		    				}
    		    			}
    					}
    					edge = edge->next;
    				}
    				// right -> left
    				if(right[rpos].sameside) edge = paths[current].leftPath;
    				else edge = paths[current].rightPath;
    				while(edge && !edge->sibl){
    					if(edge->depth == rpos + lelem +2){
    						left[lelem].ID = edge->ID;
    						if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
    						else left[lelem].sameside = 0;
    						lelem++;
    		    			if(lelem == lmaxelem){
    		    				lmaxelem *= 2;
    		    				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
    		    				if(!left){
    		    					printf("Error in realloc lmaxelem in scaffold_init\n");
    		    					exit(1);
    		    				}
    		    			}
    					}
    					edge = edge->next;
    				}
    				rpos++;
    			}
    		}
        	if(verbose){
        		int a = lpos-1;
        		printf("Scaffold Nodes (%i):\n",i);
        		while(a>-1){
        			printf("%i (%i) ",left[a].ID,left[a].sameside);
        			a--;
        		}
        		printf("<- %i -> ",i);
        		a = 0;
        		while(a<rpos){
        			printf("%i (%i) ",right[a].ID,right[a].sameside);
        			a++;
        		}
        		printf("\n");
        	}
        	// save results in aS
        	depth = 0;
        	len = 0;
        	int ID;
        	struct scaffEdge* scaffedge;
        	struct scaffEdge* scaffedgenew;
        	aS->scaff[aS->num].ID = aS->num;
        	if(lpos){
    			aS->scaff[aS->num].type = 0;
        		lpos--;
        		scaffedge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        		ID = left[lpos].ID;
        		len += paths[ID].len;
        		scaffedge->ID = ID;
        		scaffedge->depth = depth;
        		scaffedge->len = len;
        		scaffedge->next = NULL;
        		if(left[lpos].sameside){
        			aS->scaff[aS->num].startJunction = paths[ID].leftJunction;
        			aS->scaff[aS->num].endJunction = paths[ID].rightJunction;
        			scaffedge->targetJunction = paths[ID].rightJunction;
        		}
        		else{
        			aS->scaff[aS->num].startJunction = paths[ID].rightJunction;
        			aS->scaff[aS->num].endJunction = paths[ID].leftJunction;
        			scaffedge->targetJunction = paths[ID].leftJunction;
        		}
        		depth ++;
        		aS->scaff[aS->num].type = 1;
        		aS->scaff[aS->num].first = scaffedge;
        		while(lpos){
        			lpos--;
        			scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        			ID = left[lpos].ID;
        			scaffedgenew->ID = ID;
        			scaffedgenew->len = paths[ID].len;
        			scaffedgenew->depth = depth;
        			scaffedgenew->next = NULL;
        			if(left[lpos].sameside) scaffedgenew->targetJunction = paths[ID].rightJunction;
        			else scaffedgenew->targetJunction = paths[ID].leftJunction;
        			depth ++;
        			len += scaffedgenew->len;
        			scaffedge->next = scaffedgenew;
        			scaffedge = scaffedge->next;
        		}

        		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        		scaffedgenew->ID = i;
        		scaffedgenew->next = NULL;
        		scaffedgenew->depth = depth;
        		scaffedgenew->targetJunction = paths[i].rightJunction;
        		scaffedgenew->len = paths[i].len;
        		scaffedge->next = scaffedgenew;
        		scaffedge = scaffedgenew;
        		len += paths[i].len;
        		depth++;
           	}
        	else{
        		len += paths[i].len;
        		scaffedge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        		scaffedge->ID = i;
        		scaffedge->next = NULL;
        		scaffedge->depth = depth;
        		scaffedge->targetJunction = paths[i].rightJunction;
        		scaffedge->len = paths[i].len;
        		aS->scaff[aS->num].type = 1;
        		aS->scaff[aS->num].endJunction = paths[i].rightJunction;
        		aS->scaff[aS->num].startJunction = paths[i].leftJunction;
        		aS->scaff[aS->num].first = scaffedge;
        		depth++;
        	}
        	rpos = 0;
        	while(rpos < relem){
    			aS->scaff[aS->num].type = 0;
        		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        		ID = right[rpos].ID;
        		scaffedgenew->ID = ID;
        		scaffedgenew->next = NULL;
        		scaffedgenew->depth = depth;
        		scaffedgenew->len = paths[ID].len;
        		if(right[rpos].sameside) scaffedgenew->targetJunction = paths[ID].rightJunction;
        		else scaffedgenew->targetJunction = paths[ID].leftJunction;
        		scaffedge->next = scaffedgenew;
        		scaffedge = scaffedgenew;
        		len += paths[ID].len;
        		depth++;
        		rpos++;
        	}
        	aS->scaff[aS->num].len = len;
        	aS->num++;
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

struct scaffold_set* scaffold_init2(){
    int i,j;
    struct scaffold_set* aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
    aS->num = 0;
    aS->nummax = 1000;
    aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*aS->nummax);

    char verbose = 0;

    int depth = 0;
    int len = 0;

    int lpos = 0;
    int rpos = 0;
    int lelem = 0;
    int relem = 0;
    int lmaxelem = 1000;
    int rmaxelem = 1000;
    int current;
    struct contigScaff* left = (struct contigScaff*)malloc(sizeof(struct contigScaff)*lmaxelem);
    struct contigScaff* right = (struct contigScaff*)malloc(sizeof(struct contigScaff)*rmaxelem);
    struct pathEdge* edge; // = (struct pathEdge*)malloc(sizeof(struct pathEdge));

    for(j=0;j<2;j++){
        for(i=1;i<pathsNum;i++){
        	if(!paths[i].flag){
        		if(j || ((paths[i].scaffflag & 32) && (paths[i].scaffflag & 2))){
        			continue;
        		}
        		if(verbose) printf("\t NEW SCAFFOLD\n");
        		paths[i].flag++;
        		// initial left
        		lpos = 0;
        		rpos = 0;
        		lelem = 0;
        		relem = 0;
        		edge = paths[i].leftPath;
        		while(edge && !paths[edge->ID].flag && !edge->sibl){
        			if(verbose) printf("(%i) Set left %i\n",i,edge->ID);
        			left[lelem].ID = edge->ID;
        			if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
        			else left[lelem].sameside = 0;
        			if(edge->junctionCon>=0) left[lelem].bridge = edge;
        			else left[lelem].bridge = NULL;
        			lelem++;
        			if(lelem == lmaxelem){
        				lmaxelem *= 2;
        				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
        				if(!left){
        					printf("Error in realloc lmaxelem in scaffold_init\n");
        					exit(1);
        				}
        			}
        			edge = edge->next;
        		}
        		// initial right
        		edge = paths[i].rightPath;
        		while(edge && !paths[edge->ID].flag && !edge->sibl){
        			if(verbose) printf("(%i) Set right %i\n",i,edge->ID);
        			right[relem].ID = edge->ID;
        			if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
        			else right[relem].sameside = 0;
        			if(edge->junctionCon>=0) right[relem].bridge = edge;
        			else right[relem].bridge = NULL;
        			relem++;
        			if(relem == rmaxelem){
        				rmaxelem *= 2;
        				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
        				if(!right){
        					printf("Error in realloc rmaxelem in scaffold_init\n");
        					exit(1);
        				}
        			}
        			edge = edge->next;
        		}
        		while(lpos<lelem || rpos<relem){
        			if(verbose) printf("Somithing was set, go deeper\n");
        			// left
        			while(lpos<lelem){
        				current = left[lpos].ID;
        				paths[current].flag++;
        				// left -> left
        				if(left[lpos].sameside) edge = paths[current].leftPath;
        				else edge = paths[current].rightPath;
        				while(edge && !paths[edge->ID].flag && !edge->sibl){
        					if(edge->depth + lpos == lelem){
        						if(verbose) printf("(%i) (lpos: %i) Set left left %i (target: %i)\n",i,lpos,edge->ID,edge->targetJunction);
        						left[lelem].ID = edge->ID;
        		    			if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
        		    			else left[lelem].sameside = 0;
        		    			if(edge->junctionCon>=0) left[lelem].bridge = edge;
        		    			else left[lelem].bridge = NULL;
        		    			lelem++;
        		    			if(lelem == lmaxelem){
        		    				lmaxelem *= 2;
        		    				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
        		    				if(!left){
        		    					printf("Error in realloc lmaxelem in scaffold_init\n");
        		    					exit(1);
        		    				}
        		    			}
        					}
        					edge = edge->next;
        				}
        				// left -> right
        				if(left[lpos].sameside) edge = paths[current].rightPath;
        				else edge = paths[current].leftPath;
        				while(edge && !paths[edge->ID].flag && !edge->sibl){
        					if(edge->depth == lpos + relem +2){
        						if(verbose) printf("(%i) (lpos: %i) Set left right %i",i,lpos,edge->ID);
        						right[relem].ID = edge->ID;
        		    			if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
        		    			else right[relem].sameside = 0;
        		    			if(edge->junctionCon>=0) right[relem].bridge = edge;
        		    			else right[relem].bridge = NULL;
        		    			relem++;
        		    			if(relem == rmaxelem){
        		    				rmaxelem *= 2;
        		    				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
        		    				if(!right){
        		    					printf("Error in realloc rmaxelem in scaffold_init\n");
        		    					exit(1);
        		    				}
        		    			}
        					}
        					edge = edge->next;
        				}
        				lpos++;
        			}
        			// right
        			while(rpos < relem){
        				current = right[rpos].ID;
        				paths[current].flag++;
        				// right -> right
        				if(right[rpos].sameside) edge = paths[current].rightPath;
        				else edge = paths[current].leftPath;
        				while(edge && !paths[edge->ID].flag && !edge->sibl){
        					if(edge->depth + rpos == relem){
        						right[relem].ID = edge->ID;
        						if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
        						else right[relem].sameside = 0;
        		    			if(edge->junctionCon>=0) right[relem].bridge = edge;
        		    			else right[relem].bridge = NULL;
        						relem++;
        		    			if(relem == rmaxelem){
        		    				rmaxelem *= 2;
        		    				right = (struct contigScaff*)realloc(right,sizeof(struct contigScaff)*rmaxelem);
        		    				if(!right){
        		    					printf("Error in realloc rmaxelem in scaffold_init\n");
        		    					exit(1);
        		    				}
        		    			}
        					}
        					edge = edge->next;
        				}
        				// right -> left
        				if(right[rpos].sameside) edge = paths[current].leftPath;
        				else edge = paths[current].rightPath;
        				while(edge && !paths[edge->ID].flag && !edge->sibl){
        					if(edge->depth == rpos + lelem +2){
        						left[lelem].ID = edge->ID;
        						if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
        						else left[lelem].sameside = 0;
        		    			if(edge->junctionCon>=0) left[lelem].bridge = edge;
        		    			else left[lelem].bridge = NULL;
        						lelem++;
        		    			if(lelem == lmaxelem){
        		    				lmaxelem *= 2;
        		    				left = (struct contigScaff*)realloc(left,sizeof(struct contigScaff)*lmaxelem);
        		    				if(!left){
        		    					printf("Error in realloc lmaxelem in scaffold_init\n");
        		    					exit(1);
        		    				}
        		    			}
        					}
        					edge = edge->next;
        				}
        				rpos++;
        			}
        		}
            	if(verbose){
            		int a = lpos-1;
            		printf("Scaffold Nodes (%i):\n",i);
            		while(a>-1){
            			printf("%i (%i) ",left[a].ID,left[a].sameside);
            			a--;
            		}
            		printf("<- %i -> ",i);
            		a = 0;
            		while(a<rpos){
            			printf("%i (%i) ",right[a].ID,right[a].sameside);
            			a++;
            		}
            		printf("\n");
            	}
            	// save results in aS
            	depth = 0;
            	len = 0;
            	int ID;
            	struct scaffEdge* scaffedge;
            	struct scaffEdge* scaffedgenew;
            	aS->scaff[aS->num].ID = aS->num;
        		aS->scaff[aS->num].type = 1;
        		struct pathEdge* pathedge;
            	if(lpos){
        			aS->scaff[aS->num].type = 0;
            		lpos--;
            		scaffedge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
            		ID = left[lpos].ID;
            		len += paths[ID].len;
            		scaffedge->ID = ID;
            		scaffedge->depth = depth;
            		scaffedge->len = len;
            		scaffedge->next = NULL;
            		scaffedge->bridge = NULL;
            		pathedge = left[lpos].bridge;
            		if(left[lpos].sameside){
            			aS->scaff[aS->num].startJunction = paths[ID].leftJunction;
            			aS->scaff[aS->num].endJunction = paths[ID].rightJunction;
            			scaffedge->targetJunction = paths[ID].rightJunction;
            		}
            		else{
            			aS->scaff[aS->num].startJunction = paths[ID].rightJunction;
            			aS->scaff[aS->num].endJunction = paths[ID].leftJunction;
            			scaffedge->targetJunction = paths[ID].leftJunction;
            		}
            		depth ++;
            		aS->scaff[aS->num].first = scaffedge;
            		while(lpos){
            			lpos--;
            			scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
            			ID = left[lpos].ID;
            			scaffedgenew->ID = ID;
            			scaffedgenew->len = paths[ID].len;
            			scaffedgenew->depth = depth;
            			scaffedgenew->next = NULL;
            			scaffedgenew->bridge = pathedge;
            			pathedge = left[lpos].bridge;
            			if(left[lpos].sameside) scaffedgenew->targetJunction = paths[ID].rightJunction;
            			else scaffedgenew->targetJunction = paths[ID].leftJunction;
            			depth ++;
            			len += scaffedgenew->len;
            			scaffedge->next = scaffedgenew;
            			scaffedge = scaffedge->next;
            		}
            		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
            		scaffedgenew->ID = i;
            		scaffedgenew->next = NULL;
            		scaffedgenew->depth = depth;
            		scaffedgenew->targetJunction = paths[i].rightJunction;
            		scaffedgenew->len = paths[i].len;
            		scaffedgenew->bridge = pathedge;
            		scaffedge->next = scaffedgenew;
            		scaffedge = scaffedgenew;
            		len += paths[i].len;
            		depth++;
               	}
            	else{
            		len += paths[i].len;
            		scaffedge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
            		scaffedge->ID = i;
            		scaffedge->next = NULL;
            		scaffedge->depth = depth;
            		scaffedge->targetJunction = paths[i].rightJunction;
            		scaffedge->len = paths[i].len;
            		scaffedge->bridge = NULL;
    //        		aS->scaff[aS->num].type = 1;
            		aS->scaff[aS->num].endJunction = paths[i].rightJunction;
            		aS->scaff[aS->num].startJunction = paths[i].leftJunction;
            		aS->scaff[aS->num].first = scaffedge;
            		depth++;
            	}
            	rpos = 0;
            	while(rpos < relem){
        			aS->scaff[aS->num].type = 0;
            		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
            		ID = right[rpos].ID;
            		scaffedgenew->ID = ID;
            		scaffedgenew->next = NULL;
            		scaffedgenew->depth = depth;
            		scaffedgenew->len = paths[ID].len;
            		// ____________________
            		if(right[rpos].bridge) scaffedgenew->bridge = right[rpos].bridge;
            		else scaffedgenew->bridge = NULL;
            		// ____________________
            		if(right[rpos].sameside) scaffedgenew->targetJunction = paths[ID].rightJunction;
            		else scaffedgenew->targetJunction = paths[ID].leftJunction;
            		scaffedge->next = scaffedgenew;
            		scaffedge = scaffedgenew;
            		len += paths[ID].len;
            		depth++;
            		rpos++;
            	}
            	aS->scaff[aS->num].len = len;
            	aS->num++;
            	if(aS->num == aS->nummax){
            		aS->nummax *= 2;
            		aS->scaff = (struct scaffold*)realloc(aS->scaff,sizeof(struct scaffold)*aS->nummax);
            		if(!aS->scaff){
            			printf("No realloc of struct scaffold possible\n Abort\n");
            			exit(1);
            		}
            	}

        	}
        	// Rest singletons
        }
    }

    free(left);
    free(right);
	return aS;
}

static inline void resetLetters(struct Letter_T* Letters){
	int i;
	struct Letter_T* current;
	struct LetterEdge* edge;
	struct LetterEdge* nextedge;
//	struct LetterSource_S* sedge;
//	struct LetterSource_S* nextsedge;
	for(i=0;i<numNodes;i++){
//				printf("Letter NUmber %i\n",i);
		current = &Letters[i];
		edge = current->left;
		while(edge){
//					printf("Free left edge\n");
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
		edge = current->right;
		while(edge){
//					printf("Free right edge\n");
//					printf("Free right edge at pos: %i\n",i);
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
//		sedge = &current->source;
//		if(sedge){
//			sedge = sedge->next;
//			while(sedge){
////						printf("Free source\n");
//				nextsedge = sedge->next;
//				free(sedge);
//				sedge = nextsedge;
//			}
//		}
		current->vFlag = 0;
		current->counter = 0;
//		current->source.next = NULL;
		current->align_ring = NULL;
		current->junction = 0;
		current->score = 0;
		current->left = NULL;
		current->right = NULL;
	}
}

/**
 * Same implementation of the pao algorithms but graph touring over the paths instead of the junctions. Spanning junctions of scaffolds span over
 * @return
 */
struct POG* make_poaScaff(struct myovlList* G, struct reads* reads, char scaffolding, struct para* para){
	char verbose = 0;
	char verbose2 = 0;

	struct scaffold_set* aS;
	if(scaffolding){
		printf("Checkpoint: Init Scaffold Correction\n");
		aS = scaffold_init2();
	}
	else{
		printf("Checkpoint: Init Contig Correction\n");
		aS = contigs_init(G); // ,reads
	}

//	if(verbose)
		scaffold_stats(aS);
//	exit(1);

    int i,j;
    int breadID;
    struct bread* bread;
    struct bread* internb;
    printf("MaxReadLen: %i\n",maxReadLen);
    char* readseq = (char*)malloc(sizeof(char)*(maxReadLen+1));
    char* revreadseq = (char*)malloc(maxReadLen+1);
    char* name = (char*)malloc(1000);
    char* dotPath = (char*)malloc(1000);

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
//    	Letters[i].source.next = NULL;
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
    int finJunction;
//    struct scaffEdge* scaffedge;
//    struct bread* counterbread;


    //__________ Scaffold Consensus Calling
    struct scaffEdge* scaffEdge;
    int startJunction;
    for(i=0;i<aS->numbridge;i++){
    	if(aS->scaff[i].len > MIN_SCAFF_LEN || i >= aS->num){
    		scaffEdge = aS->scaff[i].first;
    		printf("FirstEdge: %i\n",scaffEdge->ID);
    		startJunction = aS->scaff[i].startJunction;
    		dir = G->read[startJunction]->dir;
    		bread = G->read[startJunction]->first;
    		// Scaffold is a Singleton; Just on Contig (no scaffEdge)
    		if(aS->scaff[i].type == 1){
    			finJunction = aS->scaff[i].endJunction;
    			while(bread){
    				if(verbose && bread->dest) printf("destID: %i == endJunction: %i\n",bread->dest->ID, aS->scaff[i].endJunction);
    				if(bread->dest && bread->dest->ID == aS->scaff[i].endJunction) break;
    				else bread = bread->next;
    			}
//    			if(bread->dest->ID != 284801){
//    				continue;
//    			}
    		}
    		else{
    			while(scaffEdge){
    				finJunction = scaffEdge->targetJunction;
    				scaffEdge = scaffEdge->next;
    			}
    			scaffEdge = aS->scaff[i].first;
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
   				poa_initBackbone2(&pog->contig[pog->contigNum],readseq); // Parameter 2: &reads[aS->scaff[i].startJunction],
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
    			poa_catBackbone(&pog->contig[pog->contigNum],G,readseq,startJunction,breadID); // Parameter 3: &reads[breadID],
    			poa_heuristic_align2(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts,overhang);
    			inserts++;
    			internb = bread;
//				while(breadID != aS->scaff[i].endJunction){
				while(breadID != finJunction){
					if(G->read[breadID]->flag == JUNCTION){
						printf("bread is a JUNCTION\n");
						scaffEdge = scaffEdge->next;
						if(scaffEdge->bridge){

						}

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
	    						poa_catBackbone2(&pog->contig[pog->contigNum],G,readseq,oldbreadID,breadID); // Parameter 3: &reads[breadID],
	//        						poa_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
	    						if(verbose) printf("ALINING PROPER READ\n");
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
			        			if(verbose2) printf("ALIGING CONTAINED READ\n");
			        			if(verbose2) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose2) printf("Read: %s\n",readseq);
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
	    						poa_catBackbone2(&pog->contig[pog->contigNum],G,readseq,oldbreadID,breadID); // Parameter 3: &reads[breadID],
	//        						poa_align(&pog->contig[pog->contigNum],&reads[breadID],readseq,1,inserts);
	//            						printf("Insert Number: %i in contig %i\n",inserts,pog->contigNum);
	    						if(verbose2) printf("ALIGING PROPER READ\n");
	    						if(verbose2) printf("c %i (%i)\n",G->read[internb->ID]->dir,internb->ID);
			        			if(verbose2) printf("Read: %s\n",readseq);
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
			sprintf(name,"Scaffold_%i_%i_%i_len:",i+1,aS->scaff[i].startJunction,breadID);
			pog->contig[pog->contigNum].name = (char*)malloc(strlen(name)+100);
			strcpy(pog->contig[pog->contigNum].name,name);

    		poa_consensus2(&pog->contig[pog->contigNum]);
    		if(verbose){
    			sprintf(dotPath,"%s/%s.dot",para->asemblyFolder,pog->contig[pog->contigNum].name);
    			poa_toDot(dotPath);
    		}

    		resetLetters(Letters);
    		numNodes = 0;
    		aS->scaff[i].scaffoldID = pog->contigNum;
    		if(aS->scaff[i].next){
    			printf("Scaffold %i has a connection\n",pog->contigNum);
    			pog->contig[pog->contigNum].seqEdge = (struct sequenceEdge*)malloc(sizeof(struct sequenceEdge));
    			pog->contig[pog->contigNum].seqEdge->insertLen = aS->scaff[i].next->first->bridge->estLen;
    			pog->contig[pog->contigNum].seqEdge->ori = 0;
    		}
    		else pog->contig[pog->contigNum].seqEdge = NULL;
    		if(i >= aS->num) pog->contig[pog->contigNum].vflag = 1;
    		else pog->contig[pog->contigNum].vflag = 0;
    		pog->contigNum++;

    		printf("Matrix time:    %.3f s\n",(float)sumMatrix/1000000000);
    		printf("Backtrace time: %.3f s\n",(float)sumTrace/1000000000);
    	}
    }

    for(i=0;i<aS->numbridge;i++){
    	if(aS->scaff[i].len > MIN_SCAFF_LEN || i >= aS->num){
    		if(aS->scaff[i].next){
    			pog->contig[aS->scaff[i].scaffoldID].seqEdge->nextScaff = aS->scaff[i].next->scaffoldID;
    			printf("Connect Scaffold %i with %i bp to scaffold %i\n",aS->scaff[i].scaffoldID,pog->contig[aS->scaff[i].scaffoldID].seqEdge->insertLen,aS->scaff[i].next->scaffoldID);
    		}
    	}
    }


    // _________ Scaffold Consensus Calling
    if(verbose)printf("Free ScaffoldsSet\n");
    free_schaffoldSet(aS);
    if(verbose)printf("Free ScaffoldsSeqs\n");
    free(readseq);
    free(revreadseq);
    free(name);
    if(verbose)printf("Free Alignment Matrix\n");
    for(i=0;i<=maxReadLen*MATRIX_MAX_BR;i++){
    	free(alMatrix[i]);
    }
    free(alMatrix);
    if(verbose)printf("Free DotPath\n");
    free(dotPath);
    if(verbose)printf("Return\n");
    return pog;
}

void free_POG(struct POG* contigs_pog){
	int i;
	struct Letter_T* current;
	struct LetterEdge* edge;
	struct LetterEdge* nextedge;
//	struct LetterSource_S* sedge;
//	struct LetterSource_S* nextsedge;
	for(i=0;i<numNodes;i++){
//				printf("Letter NUmber %i\n",i);
		current = &Letters[i];
		edge = current->left;
		while(edge){
//					printf("Free left edge\n");
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
		edge = current->right;
		while(edge){
//					printf("Free right edge\n");
//					printf("Free right edge at pos: %i\n",i);
			nextedge = edge->next;
			free(edge);
			edge = nextedge;
		}
//		sedge = &current->source;
//		if(sedge){
//			sedge = sedge->next;
//			while(sedge){
////						printf("Free source\n");
//				nextsedge = sedge->next;
//				free(sedge);
//				sedge = nextsedge;
//			}
//		}

	}
	for(i=0;i<contigs_pog->contigNum;i++){
		free(contigs_pog->contig[i].name);
		free(contigs_pog->contig[i].sequence);
	}

	free(Letters);
	free(contigs_pog->contig);
	free(contigs_pog);
	Letters = NULL;
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
