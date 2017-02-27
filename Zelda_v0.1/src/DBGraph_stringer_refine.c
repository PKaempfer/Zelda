/*
 ============================================================================
 Name        : DBGraph_stringer_refine.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Heart of Zelda: Converts the dBG to a String Graph by calculating
 	 	 	   transitively irreducible overlaps from the dBG.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "DBGraph_stringer.h"
#include "DBGraph.h"

static char status_char[] = { 'W', 'C', 'S', 'P', 'J' };
static const char *arrowdir[2]={"normal","inv"};

void radixSort(int nodeNum, int* nStat){
	int i,j;
	int temp0N,temp1N;
	int *temp0 = (int*)malloc(sizeof(int)*nodeNum);
	int *temp1 = (int*)malloc(sizeof(int)*nodeNum);
	uint32_t mask = 1;
	for(i=0;i<24;i++){ // Max Length of largest Contig 16mbp
		temp0N = 0;
		temp1N = 0;
		for(j=0;j<nodeNum;j++){
			if(nStat[j] & mask) temp0[temp0N++] = nStat[j];
			else temp1[temp1N++] = nStat[j];
		}
		memcpy(&nStat[0],&temp0[0],sizeof(int)*temp0N);
		memcpy(&nStat[temp0N],&temp1[0],sizeof(int)*temp1N);
		mask = mask << 1;
	}
}

void countRemainingNodes(){
	int i;
	int nodeNum = 0;
	int nodeNum100 = 0;
	int nodeLen = 0;
	int nodeLen100 = 0;
	for(i=1;i<=redGraph->V;i++){
		if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
			nodeNum++;
			nodeLen += redGraph->array[i].len;
			if(redGraph->array[i].len > (100-nK)){
				nodeLen100 += redGraph->array[i].len;
				nodeNum100++;
			}
		}
	}

	printf("\nGraph stats:\n");
	printf("Number of contigs: %i\n",nodeNum);
	printf("Number of contigs > 100bp: %i\n",nodeNum100);
	printf("Avg contig length: %2.f\n",(float)nodeLen/nodeNum100);

	int *nStat = (int*)malloc(sizeof(int)*(nodeNum+1));
	int k=0;
	for(i=1;i<=redGraph->V;i++){
		if(redGraph->array[i].head || redGraph->array[i].tail || redGraph->array[i].headread){
			nStat[k++]=redGraph->array[i].len+k;
		}
	}

	radixSort(nodeNum,nStat);

	printf("Contig Length: %i\n",nodeLen);
	printf("Contig Length (>100): %i\n",nodeLen100);
	printf("Largest Contig: %i\n\n", nStat[0]);

	int sum=0;
	int ns=0;
	for(i=0;i<nodeNum;i++){
		sum += nStat[i];
		if(sum > (nodeLen/10) && ns == 0){
			printf("N10: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/4) && ns == 1){
			printf("N25: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/2) && ns == 2){
			printf("N50: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/4)*3 && ns == 3){
			printf("N75: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/10)*9 && ns == 4){
			printf("N90: %i\n\n",nStat[i]);
			ns++;
		}
	}

	free(nStat);
}

// recursive depth first search
void tagComponent(struct string_graph *S, int v, int *nStat, char *edgeflag){
	int  w, e, o;
	int  w2,e2,o2;
	int t;
	for(o = 0; o < 2; o++){
		w = v + o;
		for (e = S->side[w - 1].top + 1; e <= S->side[w].top; e++){
			if(!edgeflag[e]){
				edgeflag[e] = 1;
				nStat[e] = S->edge[e].length;
				t = S->edge[e].target;
				if(t%2==0) t-=1;

				for(o2 = 0; o2 < 2; o2++){
					w2 = t + o2;
					for (e2 = S->side[w2 - 1].top + 1; e2 <= S->side[w2].top; e2++){
						if(VERT(S->edge[e2].target) == VERT(v)){
							edgeflag[e2] = 1;
						}
					}
				}

				if(S->ID[VERT(t)] > 0){
					S->ID[VERT(t)] *=-1;
					tagComponent(S,t,nStat,edgeflag);
				}
			}
		}

	}
}

void stringGraphStats(struct string_graph *S){
    int comp = 0;														// Number of components
    int nodeLen = 0;													// Sum over all Contig length
	int *nStat = (int*)malloc(sizeof(int)*(2*S->nedges+2));			// List of Contig length
	char *edgeflag = (char*)malloc(sizeof(char)*(2*S->nedges)+3);
	int k,i;

	for(i=0;i<=(2*S->nedges);i++) edgeflag[i] = 0;
	for(i=0;i<=(2*S->nedges);i++) nStat[i] = 0;

    int nends  = S->nends;
    unsigned char*  status = S->status;
    int       v;
    for (v = 1; v <= nends; v += 2)
    {
        if ( status[VERT(v)] == WIDOWED ||
             status[VERT(v)] == CONTAINED )
        {
            continue;
        }

        // Use ID as visited-flag
        if(status[VERT(v)]==JUNCTION && S->ID[VERT(v)] > 0){
        	S->ID[VERT(v)] *=-1;
        	tagComponent(S,v,nStat,edgeflag);
        	comp++;
        }
    }

    // Restore real IDs
    for (v = 1; v <= nends; v += 2){
        if(status[VERT(v)]==JUNCTION){
        	S->ID[VERT(v)] *=-1;
        }
    }

	int temp;
	k = (2*S->nedges);
	int newk;
	// Sort Contigs (Pathlength between junctions) in String graph by length
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

	for(i=0;i<=(2*S->nedges);i++) nodeLen += nStat[i];

	int sum=0;
	int ns=0;

	printf("Largest Contig: %i bp\n",nStat[0]);
	printf("Number of Components: %i\n",comp);
	printf("Total Length over all edges: %i\n",nodeLen);
	for(i=0;i<(2*S->nedges);i++){
		sum += nStat[i];
		if(sum > (nodeLen/10) && ns == 0){
			printf("N10: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/4) && ns == 1){
			printf("N25: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/2) && ns == 2){
			printf("N50: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/4)*3 && ns == 3){
			printf("N75: %i\n",nStat[i]);
			ns++;
		}
		if(sum > (nodeLen/10)*9 && ns == 4){
			printf("N90: %i\n\n",nStat[i]);
			ns++;
		}
	}

	free(nStat);
	free(edgeflag);

}

struct myovlList* initOVLgraph(int numreads){
	struct myovlList* ovlgraph = (struct myovlList*)malloc(sizeof(struct myovlList));
	ovlgraph->V = numreads;
	ovlgraph->read = (struct aread**)malloc(sizeof(struct aread*)*(numreads+2)); // 0 ist NULL, last is \0->pointer
	printf("Number of all reads: %i\n",ovlgraph->V);
	int i;
	for(i=0;i<=numreads;i++){
		ovlgraph->read[i] = NULL;
	}
	return ovlgraph;
}

struct string_graph* initStringGraph(struct myovlList* ovlgraph,char* pathAssembly, char dotdump){
//	stringer2(ovlgraph);
	printf("FIND INITIAL CONTAINMENTS\n");
	tag_A_Contained(ovlgraph);
	printf("CALCULATE OVERLAPS\n");
	stringer3(ovlgraph);
//	printf("CATEGORIZE OVERLAPS\n");
	struct string_graph* S = catOVLgraph(ovlgraph, pathAssembly);
//	catOVLgraph(ovlgraph);
	if(dotdump){
		char* tempPath = (char*)malloc(200);
		printf("WRITE DOT-FILES\n");
		sprintf(tempPath,"%s/catOrgOVL.dot",pathAssembly);
		printOVLgraph(ovlgraph,1,tempPath);
		sprintf(tempPath,"%s/catRevOVL.dot",pathAssembly);
		printOVLgraph(ovlgraph,0,tempPath);
		sprintf(tempPath,"%s/catRedRevOVL.dot",pathAssembly);
		printReducedOVLgraph(ovlgraph,1,tempPath);
		sprintf(tempPath,"%s/string.dot",pathAssembly);
		printStringGraph(ovlgraph,tempPath);
		free(tempPath);
	}
	return S;
}


void printStringGraph(struct myovlList *ovlGraph, char *filename){
	FILE* string = fopen(filename,"w");
	int i, j;
	int dir;
	int len;
	int num;
	float dens = 0;
	struct bread* bread;
	struct bread* internb;

	fprintf(string,"digraph StringGraph {\n");
//	printf("Number of nodes in the StringGraph: %i\n",ovlGraph->V);

	for(i = 0; i <= ovlGraph->V; i++){
		if(ovlGraph->read[i] && ovlGraph->read[i]->flag == JUNCTION){
			fprintf(string,"%i [label=\"%i (%i)\"]\n",i,i,readStartList[i]);
			bread = ovlGraph->read[i]->first;
			while(bread){
				len = 0;
				if(ovlGraph->read[bread->ID]->flag != CONTAINED){
					num=0;
					dir = bread->sideflag;
//					printf("\t%i (%i)\n",bread->ID,bread->overhang);
					len+=bread->overhang;
					j = bread->ID;
					while(ovlGraph->read[j]->flag != JUNCTION){
						internb = ovlGraph->read[j]->first;
						while(internb){
							if(internb->sideflag == dir && ovlGraph->read[internb->ID]->flag != CONTAINED){
								len += internb->overhang;
								num++;
								j = internb->ID;
//								printf("\t%i (%i)\n",j,internb->overhang);
								break;
							}
							internb = internb->next;
						}
					}
//					printf("%i -> %i[label=\"%i\"]\n",i,j,len);
					dens = (float)(num+1)/len;
					if(num || j!=i) fprintf(string,"%i -> %i[label=\"%i(%i) -> %0.3f\", fontcolor=\"blue\"];\n",i,j,len,num,dens);
				}
				bread = bread->next;
			}
		}
	}
	fprintf(string,"}\n");
	fclose(string);
}

void printOVLgraph(struct myovlList *ovlGraph,int dir, char *ovlPath){
	int i;
	struct bread* bread;
	FILE *ovl = fopen(ovlPath,"w");
	fprintf(ovl,"digraph ovl {\n");

	for (i=0;i<=ovlGraph->V;i++){
		if(ovlGraph->read[i]){
			fprintf(ovl,"%i [label=\"%i (%i , %i)\"];\n",i,i,ovlGraph->read[i]->length,ovlGraph->read[i]->flag);
			bread = ovlGraph->read[i]->first;
			while(bread){
				if(bread->sideflag==dir){
					fprintf(ovl,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, label=\"%i\"];\n",i,bread->ID,arrowdir[(int)!ovlGraph->read[i]->dir],arrowdir[(int)ovlGraph->read[bread->ID]->dir],bread->overhang);
				}
				bread = bread->next;
			}
		}
	}
	fprintf(ovl,"}\n");
	fclose(ovl);
}

void printReducedOVLgraph(struct myovlList *ovlGraph,int dir, char *ovlPath){
	int i;
	struct bread* bread;
	FILE *ovl = fopen(ovlPath,"w");
	fprintf(ovl,"digraph ovl {\n");

	for (i=0;i<=ovlGraph->V;i++){
		if(ovlGraph->read[i] && ovlGraph->read[i]->flag!=CONTAINED){
			fprintf(ovl,"%i [label=\"%i (%i , %i)\"];\n",i,i,ovlGraph->read[i]->length,ovlGraph->read[i]->flag);
			bread = ovlGraph->read[i]->first;
			while(bread){
				if(bread->sideflag==dir && ovlGraph->read[bread->ID]->flag!=CONTAINED){
					fprintf(ovl,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, label=\"%i\"];\n",i,bread->ID,arrowdir[(int)!ovlGraph->read[i]->dir],arrowdir[(int)ovlGraph->read[bread->ID]->dir],bread->overhang);
				}
				bread = bread->next;
			}
		}
	}
	fprintf(ovl,"}\n");
	fclose(ovl);
}

void checkTransitivity(struct myovlList *ovlGraph){
	printf("CHECK OVERLAP-TRANSITIVETY\n");
	int i;
	int in;
	int target;
	struct bread *aread;
	struct bread *bread;

	for(i=0;i<=ovlGraph->V;i++){
		if(ovlGraph->read[i]){
			aread = ovlGraph->read[i]->first;
			while(aread){
				target = aread->ID;
				bread = ovlGraph->read[target]->first;
				in = 0;
				while(bread){
					if(bread->ID == i){
						in = 1;
						break;
					}
					bread = bread->next;
				}
				if(!in){
					printf("Leak of transitivity at %i -> %i\n",i,target);
					printf("\n");
					printf("Read: %i (%i)\n",i,ovlGraph->read[i]->flag);
					bread = ovlGraph->read[i]->first;
					while(bread){
						printf("\t-> %i\t%i\t%i  -> dir: %i  -> flag: %i\n",bread->ID,ovlGraph->read[i]->dir != ovlGraph->read[bread->ID]->dir,bread->overhang,(int)bread->sideflag,ovlGraph->read[bread->ID]->flag);
						bread = bread->next;
					}
					printf("\n");
					i=target;
					printf("Read: %i (%i)\n",i,ovlGraph->read[i]->flag);
					bread = ovlGraph->read[i]->first;
					while(bread){
						printf("\t-> %i\t%i\t%i  -> dir: %i  -> flag: %i\n",bread->ID,ovlGraph->read[i]->dir != ovlGraph->read[bread->ID]->dir,bread->overhang,(int)bread->sideflag,ovlGraph->read[bread->ID]->flag);
						bread = bread->next;
					}
					printf("\n");
					return;
				}
				aread = aread->next;
			}
			in = 0;

		}
	}
	printf("String graph is transitively closed\n");
}

struct string_graph* catOVLgraph(struct myovlList *ovlGraph, char* pathAssembly){
	char verbose = 0;
	int i;
	struct bread *bread;
//	struct KannteNode* kannte;

	// Delete double edges and prop-to-contained OVLs not being transitively supported
//	struct bread *revbread;
//	struct bread *tempread;
//	for(i=0;i<=ovlGraph->V;i++){
//		if(ovlGraph->read[i]){
//			if(ovlGraph->read[i]->flag==PROPER){
//				if((bread = ovlGraph->read[i]->first)){
//					while(bread && ovlGraph->read[bread->ID]->flag == CONTAINED && ovlGraph->read[bread->ID]->first->ID != i){
////						printf("first delete Edge %i -> %i (flag: %i)\n",i,bread->ID,ovlGraph->read[bread->ID]->flag);
//						ovlGraph->read[i]->first = bread->next;
//						free(bread);
//						bread = ovlGraph->read[i]->first;
//					}
//					while(bread && bread->next){
//						if(ovlGraph->read[bread->next->ID]->flag == CONTAINED && ovlGraph->read[bread->next->ID]->first->ID != i){
////							printf("delete Edge %i -> %i (flag: %i)\n",i,bread->next->ID,ovlGraph->read[bread->next->ID]->flag);
//							tempread = bread->next->next;
//							free(bread->next);
//							bread->next = tempread;
//							continue;
//						}
//						bread = bread->next;
//					}
//				}
//			}
//			if((bread = ovlGraph->read[i]->first)){
//				while(bread){
//					revbread = bread;
//					while(revbread->next){
//						if(bread->ID == revbread->next->ID && bread->overhang == revbread->next->overhang){
////							printf("Doubled edge (%i -> %i)\n",i,bread->ID);
//							tempread = revbread->next->next;
//							free(revbread->next);
//							revbread->next = tempread;
//						}
//						else revbread = revbread->next;
//					}
//					bread = bread->next;
//				}
//			}
//		}
//	}

	checkTransitivity(ovlGraph);

	printf("CHECKPOINT: Read categorizer\n");
	int wid = 0;
	int cont = 0;
	int str = 0;
	int prop = 0;
	int junc = 0;
	int all=0;
	int in,out;
	int totin = 0, totout = 0;

	for(i=0;i<=ovlGraph->V;i++){
		if(ovlGraph->read[i] && ovlGraph->read[i]->flag != CONTAINED){
			in = 0; out = 0;
			if((bread = ovlGraph->read[i]->first)){

				while(bread){
					if(ovlGraph->read[bread->ID]->flag != CONTAINED){
						if(bread->sideflag) in++;
						else out++;
					}

					bread = bread->next;
				}
			}
			// No overlaps
			if(!ovlGraph->read[i]->first){
//				printf("WIDOWED Read: %i\n",i);
				ovlGraph->read[i]->flag = WIDOWED;
				wid++;
			}
			else if(in == 0 && out == 0){
				ovlGraph->read[i]->flag = STRAY;
				str++;
			}
			else if(in != 1 || out != 1){
				totin += in;
				totout += out;
				ovlGraph->read[i]->flag = JUNCTION;
				junc++;
			}
			else if(in == 1 && out == 1) prop++;
		}
		else if(ovlGraph->read[i]){
			cont++;
		}
	}

	all = wid+cont+str+prop+junc;

	if(verbose){
		printf("Sting graph as Table:\n");
		printf("Correct memopos of ovlgraph: %p\n",ovlGraph);

		for(i=0;i<=ovlGraph->V;i++){
			if(ovlGraph->read[i] && ovlGraph->read[i]->flag == JUNCTION){// && i == 126342){
				printf("Read: %i (%i) -> Dir: %i\n",i,ovlGraph->read[i]->flag,ovlGraph->read[i]->dir);
				bread = ovlGraph->read[i]->first;
				while(bread){
					if(bread->dest){
						printf("\t-> %i (destJ: %i)\t%i\t%i  -> dir: %i  -> flag: %i\n",bread->ID,bread->dest->ID,ovlGraph->read[i]->dir != ovlGraph->read[bread->ID]->dir,bread->overhang,(int)bread->sideflag,ovlGraph->read[bread->ID]->flag);
					}
					else{
						printf("\t-> %i\t%i\t%i  -> dir: %i  -> flag: %i\n",bread->ID,ovlGraph->read[i]->dir != ovlGraph->read[bread->ID]->dir,bread->overhang,(int)bread->sideflag,ovlGraph->read[bread->ID]->flag);
					}
					bread = bread->next;
				}
			}
		}
	}


	printf("StringGraph Stats: \n");
	printf("Widowed:\t\t\t%i\n",wid);
	printf("Contained:\t\t\t%i\n",cont);
	printf("Stray:\t\t\t\t%i\n",str);
	printf("Proper:\t\t\t\t%i\n",prop);
	printf("Junction:\t\t\t%i\n",junc);
	printf("\t-> Incomming edges:\t%i\n",totin);
	printf("\t-> Outgoing edges: \t%i\n",totout);
	printf("All:\t\t\t\t%i\n",all);

	// Init Stringgraph
	struct string_graph *S = (struct string_graph*)malloc(sizeof(struct string_graph));

	S->nverts = ovlGraph->V;
	S->nends = 2 * S->nverts;
	S->nedges = totout;
	S->edge = (struct overhang*)malloc(sizeof(struct overhang) * (totout*2+1));
	struct readend* side = (struct readend*)malloc(sizeof(struct readend) * (2*S->nverts+2));
	S->length = (int*)malloc(sizeof(int) * (S->nverts+2));
	S->ID = (int*)malloc(sizeof(int) * (S->nverts+2));
	S->status = (unsigned char*)malloc(sizeof(unsigned char) * (S->nverts+2));
	side[0].top=-1;
    for (i = 1; i <= (2*S->nverts); i++){
        side[i].top = 0;
    }


    S->side = side;

	// Fill String graph
	ovlToString(ovlGraph,S, pathAssembly);

	return S;
}

void ovlToString(struct myovlList *G, struct string_graph *S, char* pathAssembly){

	char verbose = 0;

	int i;
	int rc = 1;
	int si=1;
	int overpos = 0;
	int dir;
	int len;

	int num;
	int bID;

	struct bread* internb;
	struct bread* bread;

	int nends = S->nends;
	int* length = S->length;
	int* ID = S->ID;
	unsigned char* status = S->status;
	struct readend *side = S->side;
	struct overhang* edge = S->edge;
	struct aread** read = G->read;

//	for(i=1;i<=O->V;i++){
//		if(i == 243898){
//			bID = i;
//			bread = read[bID]->first;
//			while(bread){
//				printf("to: %i (dir: %i) (flag: %i / %i)\n",bread->ID,bread->sideflag,read[bread->ID]->flag);
//				bread = bread->next;
//			}
//		}
//	}
//	printf("RAUS\n");
//
//	return;

	struct bread* firstbread;
	int lastbreadID;

	for(i=1;i<=G->V;i++){
		rc=i;
		ID[rc] = i;
		length[rc] = readLenList[i];
		if(!G->read[i]){
			status[rc] = WIDOWED;
			si+=2;
			continue;
		}
		else if((status[rc] = G->read[i]->flag) == JUNCTION){
			if(verbose) printf("Overlaps of JUNCTION read: %i (dir: %i) (%i -> %i)\n",i,read[i]->dir,G->read[i]->stk,G->read[i]->endk);
//			bread = read[i]->first;
//			while(bread){
//				printf("to: %i (dir: %i) (flag: %i / %i)\n",bread->ID,bread->sideflag,bread->flag,read[bread->ID]->flag);
//				bread = bread->next;
//			}

			// Print nearest proper overlaps of junction reads


//			status[rc] = JUNCTION;
			if(verbose) printf("\nSIDE ONE:\n");
			dir=!read[i]->dir;
			bread = read[i]->first;
			while(bread){
				lastbreadID = i;
				len = read[i]->length;
				firstbread = bread;
				bID = bread->ID;
				if(read[bID]->flag != CONTAINED && bread->sideflag==dir){
					num=0;
					len+=bread->overhang;
					while(read[bID]->flag != JUNCTION){
						internb = read[bID]->first;
						if(verbose) {
							if(!internb) printf("Infinite loop\n");
						}
						while(internb){
							if(internb->sideflag == dir && read[internb->ID]->flag != CONTAINED){
								len += internb->overhang;
								num++;
								lastbreadID = bID;
								bID = internb->ID;
//								if(internb == internb->next) read[bID]->flag = JUNCTION;
								break;
							}
//							if(internb == internb->next) read[bID]->flag = JUNCTION;

							internb = internb->next;
						}
					}
//					dens = (float)(num+1)/len;
					side[si].top ++;
					if(verbose) printf("ovl: %i -> %i Overpos: %i (si: %i)\n",i,bID,overpos,si);
					if(read[bID]->dir != read[i]->dir) edge[overpos].target = bID*2;
					else edge[overpos].target = bID*2-1;
					edge[overpos].length = len;
					edge[overpos].nreads = num;
					// Set the destination junction of a path to the first bread in the overlap graph
					if(verbose) printf("Set Path: %i -> %i (len: %i)\n",i,bID,len);
					if(!firstbread->dest){
						firstbread->dest = (struct dest*)malloc(sizeof(struct dest));
					}
					internb = read[bID]->first;
					while(internb){
						if(internb->ID == lastbreadID && internb != bread){
							if(!internb->dest) internb->dest = (struct dest*)malloc(sizeof(struct dest));
							internb->dest->counterbread = firstbread;
							firstbread->dest->counterbread = internb;
							break;
						}
						internb = internb->next;
					}
					firstbread->dest->len = len;
					firstbread->dest->ID = bID;
					firstbread->dest->flag = 0;
					overpos++;
					if(verbose) printf("Dir 0 JUNCTION: %i (%i -> %i) OVL: %i  (%i -> %i)  /  (len: %i, AnzToJunction: %i)\n",bID,G->read[bID]->stk,G->read[bID]->endk,bread->ID,G->read[bread->ID]->stk,G->read[bread->ID]->endk,bread->overhang,num+1);
				}
				bread = bread->next;
			}
			if(verbose) printf("\nOTHER SIDE\n");
			si++;
			dir=(read[i]->dir);
			bread = read[i]->first;
			while(bread){
				lastbreadID = i;
				len = 0;
				bID = bread->ID;
				firstbread = bread;
//				printf("While bread BID  %i (flag: %i / %i):\n",bID,read[bID]->flag,read[bID]->first->ID);
				if(read[bID]->flag != CONTAINED && bread->sideflag==dir){
					num=0;
					len+=bread->overhang;
					while(read[bID]->flag != JUNCTION){
//						printf("inner While bread BID  %i (flag: %i / %i):\n",bID,read[bID]->flag,read[bID]->first->ID);
						internb = read[bID]->first;
						if(verbose)
							if(!internb) printf("Infinite loop\n");
						while(internb){
//							printf("InnerWhile\n");
							if(internb->sideflag == dir && read[internb->ID]->flag != CONTAINED){
//								printf("InnerWhileIf\n");
								len += internb->overhang;
								num++;
								lastbreadID = bID;
								bID = internb->ID;
//								if(internb == internb->next) read[bID]->flag = JUNCTION;
								break;
							}
//							if(memcmp internb == internb->next) read[bID]->flag = JUNCTION;
							internb = internb->next;
						}
//						printf("OuterWhile %i\n",internb->ID);
					}
//					dens = (float)(num+1)/len;
					side[si].top ++;
					if(verbose) printf("ovl: %i -> %i Overpos: %i (si: %i)\n",i,bID,overpos,si);
					if(read[bID]->dir != read[i]->dir) edge[overpos].target = bID*2-1;
					else edge[overpos].target = bID*2;
//					edge[overpos].target = bID*2;
					edge[overpos].length = len;
					edge[overpos].nreads = num;
					// Set the destination junction of a path to the first bread in the overlap graph
					if(verbose) printf("Set Path: %i -> %i (len: %i)\n",i,bID,len);
					if(!firstbread->dest){
						firstbread->dest = (struct dest*)malloc(sizeof(struct dest));
					}
					internb = read[bID]->first;
					while(internb){
						if(internb->ID == lastbreadID && internb != bread){
							if(!internb->dest) internb->dest = (struct dest*)malloc(sizeof(struct dest));
							internb->dest->counterbread = firstbread;
							firstbread->dest->counterbread = internb;
							break;
						}
						internb = internb->next;
					}
					firstbread->dest->len = len;
					firstbread->dest->ID = bID;
					firstbread->dest->flag = 0;
					overpos++;
				}
//				printf("Next bread\n");
				bread = bread->next;
			}
			si++;
			rc++;
//			printf("\n");
		}
		else si += 2;
	}

//    printf("Side: %i",side[0].top);
//    for (i = 1; i <= nends; i++){
//        printf(", %i",side[i].top);
//    }
//    printf("\n");

    for (i = 1; i <= nends; i++){
        side[i].top += side[i - 1].top;
    }


//    for (i = nends; i >= 1; i--){
//        side[i].top = side[i - 1].top;
//    }

//    printf("Side: %i",side[0].top);
//    for (i = 1; i <= nends; i++){
//        printf(", %i",side[i].top);
//    }
//    printf("\n");

    if(verbose) printf("CHECKPOINT: StringGraphList\n");
	if(verbose) print_string_graph_list(S,"Test_label");
	char* tempPath = (char*) malloc(200);
	sprintf(tempPath,"%s/stringGraph.dot",pathAssembly);
	print_overlap_graph_dot(S,tempPath,"StringGraph");
	free(tempPath);
}

int backedge_length(struct string_graph* G, int v1, int v2){
    // given v1 -> v2 , check if v2 -> v1 exists

    struct readend*  side   = G->side;
    struct overhang* edge   = G->edge;

    v1 = MATE(v1);
    v2 = MATE(v2);

    int e;
    for (e = side[v2 - 1].top + 1; e <= side[v2].top; e++)
    {
        if (edge[e].target == v1)
        {
            return edge[e].length;
        }
    }

    return 0;
}

char* vtx_name(int v, int width)
{
    static char name[100];

    sprintf(name, "%*d%c", width, VERT(v), (v%2) ? 'b' : 'e');

    /*
    if (v % 2)
    {
        sprintf(name, "%*db", width, VERT(v));
    }
    else
    {
        sprintf(name, "%*de", width, VERT(v));
    }
    */

    return name;
}

void print_string_graph_list(struct string_graph* G, char* label){
    int       nends  = G->nends;
    unsigned char*     status = G->status;
    struct readend*  side   = G->side;
    struct overhang* edge   = G->edge;

    int       v, w, e, o;

    printf("\n%s\n", label);

//    printf("Sides:\n");
//    for (v = 0; v <= nends; v ++)
//    {
//    	if(status[VERT(v)] == JUNCTION)
//    	printf("v: %i -> top: %i (mark %i)\n",VERT(v),side[v].top,side[v].mark);
//    }
//    printf("\n");

    for (v = 1; v <= nends; v += 2)
    {
        if ( status[VERT(v)] == WIDOWED ||
             status[VERT(v)] == CONTAINED )
        {
            continue;
        }

        if(status[VERT(v)]==JUNCTION)
        printf("%5d(%i): %c\n", VERT(v),G->ID[VERT(v)], status_char[status[VERT(v)]]);

        for (o = 0; o < 2; o++){
        	w = v + o;

            for (e = side[w - 1].top + 1; e <= side[w].top; e++){
//            	if(status[VERT(v)]==JUNCTION) printf("     %c -> %i[%3d]\n", (o ? 'e' : 'b'), edge[e].target, edge[e].length);
                printf("     %c -> %s [len: %3d]\n",
                       (o ? 'e' : 'b'), vtx_name(edge[e].target, 5), edge[e].length);
            }
        }
    }

    fflush(stdout);
}

void print_overlap_graph_dot(struct string_graph* G, const char* pcFile, const char* title){
	char verbose = 0;
    int       nends  = G->nends;
    unsigned char*     status = G->status;
    struct readend*  side   = G->side;
    struct overhang* edge   = G->edge;

    if(verbose){
        printf("Create overlap graph:\n");
        printf("nverts: %i\n",G->nverts);
        printf("nends: %i\n",G->nends);
        printf("edges: %i\n",G->nedges);
        printf("StatusSize: %i\n",G->nverts+1);
        printf("sideSize: %i\n",G->nends+1);
    }

    stringGraphStats(G);

    // v = overhang of vertex v/2
    // w = overhang left and right (w = v, w = v+1)
    // e = edge (w -> t)
    // o = 0 or 1 (left or right overhang)
    // t = target overhang (target v * 2)

    int       v, w, e, o, t, dashed, len_back;
    int bHasEdges = 0;

    FILE* fileOut = fopen(pcFile, "w");

    fprintf(fileOut, "digraph %s {\n", title);

    for (v = 1; v <= nends; v += 2)
    {
        bHasEdges = 0;

        for (o = 0; o < 2; o++)
        {
            w = v + o;

            for (e = side[w - 1].top + 1; e <= side[w].top; e++)
            {
                bHasEdges = 1;
                t = edge[e].target;

                len_back = backedge_length(G, w, t);

                if (len_back)
                {
                    if (w < t) continue;

                    dashed = 0;
                }
                else
                {
                    dashed = 1;
                }

                fprintf(fileOut, " %d -> %d [headlabel=\"%d\", color=%s, weight=%d, dir=both, divergence=%d",
                        VERT(v), VERT(t),
                        edge[e].length,
                        edge[e].reduced ? "gray" : "darkgreen",
                        edge[e].reduced ? 1 : 5,
                        edge[e].divergence);

                if (dashed)
                {
                    fprintf(fileOut, ", style=dashed");
                }
                else
                {
                    fprintf(fileOut, ", taillabel=%d", len_back);

                }

                if (w % 2)      // b
                {
                    fprintf(fileOut, ", arrowtail=normal");
                }
                else            // e
                {
                    fprintf(fileOut, ", arrowtail=inv");
                }

                if (t % 2)      // b
                {
                    fprintf(fileOut, ", arrowhead=inv");
                }
                else            // e
                {
                    fprintf(fileOut, ", arrowhead=normal");
                }

                fprintf(fileOut, "];\n");
            }
        }

        if (bHasEdges)
        {
            fprintf(fileOut, "  %d [label=%d,status=%c];\n", VERT(v), VERT(v), status_char[status[VERT(v)]]);
        }
    }

    fprintf(fileOut, "}\n");

    fclose(fileOut);
}

