/*
 ============================================================================
 Name        : DBGraph_scaffold.c
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Scaffolding unit of the DBSGA
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "DBGraph_scaffold.h"
#include "readDB.h"

/** Macro defines minimum length of paths to consider in scaffolding (should be zero I guess)*/
#define MIN_CONTIG_LEN 0

static const char* peOri[3]={"FR","RF","FF"};
static char status_char[] = { 'W', 'C', 'S', 'P', 'J' };

/** List of all paths */
struct path* paths = NULL;
/** Number of all paths*/
int pathsNum = 1;
/** Threshold for path list reallocation */
int maxPathNum = 100;
/** Number of components in String Graph*/
int componentNum;

/**
 * This Function initializes all needed data structures for paired read annotation.
 * @param G			The Overlap graph. At this point already used as a string graph, because it contains the same transitively reduces information
 * @param reads		Annotated reads.
 * @param numread	Number of all reads of all kinds of Libraries
 */
void initAnnotation(struct myovlList* G, struct reads* reads){
	char verbose = 0;
	int i;
	paths = (struct path*)malloc(sizeof(struct path) * maxPathNum);
	struct pc_anno* pc_anno;
	struct j_anno* j_anno;

	for(i=1; i <= G->V; i++){

		if(G->read[i] && G->read[i]->flag == JUNCTION){
			j_anno = (struct j_anno*)malloc(sizeof(struct j_anno));
			j_anno->vFlag = 0;
			if(verbose) printf("InitJunction Annotation for read: %i (flag:%i)\n",i,j_anno->vFlag);
			reads[i].annotation = (void*)j_anno;
		}
		else if(G->read[i]){
//			printf("Init annotation pc: %i\n",i);
			pc_anno = (struct pc_anno*)malloc(sizeof(struct pc_anno));
			pc_anno->pathID = 0;
			pc_anno->lJunctionDist = 0;
			pc_anno->rJunctionDist = 0;
			reads[i].annotation = (void*)pc_anno;
		}
	}
}

/**
 * DFS touring through the overlap graph -> recursive function call.
 * @param readID	Is the ID of the Junction read previously reached in touring by path traversing from ancient Junction. Only called if Junction was not reached by function call before.
 * @param G			The Overlap Graph (String Graph).
 * @param reads		Annotated reads of all kinds of libraries.
 */
void junctionDFS(int readID, struct myovlList* G, struct reads* reads){
	int verbose = 0;
	static int recursionDepth = 0;
	recursionDepth++;
	int i = readID;

	struct j_anno* j_anno;
	struct pc_anno* pc_anno;

	j_anno = (struct j_anno*)reads[i].annotation;
	j_anno->vFlag = 1;

	int breadID;
	struct bread* bread = G->read[i]->first;
	struct bread* internb;
	struct bread* counterbread;
	int dir = G->read[i]->dir;
	int bdir;
	int overhang;
	int counteroverhang = 0;
	int pathlength;

	if(verbose) printf("Start Annotation at Junction node: %i jDir: %i\n",i,dir);

	while(bread){
		breadID = bread->ID;
		// first edge of junction is proper or junction, has a junction destination, is not flagged as already visited, and the path has a minimum length
		if(bread && G->read[breadID]->flag != CONTAINED && bread->dest && !bread->dest->flag && bread->dest->len >= MIN_CONTIG_LEN){
			if(verbose) printf("-> Path Number: %i (len: %i) in comp: %i -> ljDir: %i (side: %i) (IDs: %i -> %i)\n",pathsNum,bread->dest->len,componentNum,dir,bread->sideflag,i,bread->dest->ID);
			// Flag path as visited
			bread->dest->flag = 1;
			pathlength = bread->dest->len;
			// Flag the other side of the unique path as already used, to not report the reverse path
			counterbread =  G->read[bread->dest->ID]->first;
			bread->dest->pathID = pathsNum;

			counterbread = bread->dest->counterbread;
			counterbread->dest->flag = 1;
			counterbread->dest->pathID = pathsNum;

			// Walk through the path, mark the reads (with their path membership and relative path position) and paths with their component membership
			bdir = bread->sideflag;
			overhang = bread->overhang;
			counteroverhang = 0;
			internb = bread;

			if(G->read[breadID]->flag == PROPER){
				counterbread = G->read[bread->ID]->first;
				while(counterbread){
					if(counterbread->ID == breadID){
						counteroverhang += counterbread->overhang;
						break;
					}
					counterbread = counterbread->next;
				}
				pc_anno = (struct pc_anno*)reads[internb->ID].annotation;
				pc_anno->pathID = pathsNum;
				pc_anno->lJunctionDist = overhang;
				pc_anno->rJunctionDist = pathlength - counteroverhang;
				if(verbose) printf("Set annotation p for : %i (path: %i / distlJ: %i / distrJ: %i)\n",bread->ID,pc_anno->pathID,pc_anno->lJunctionDist,pc_anno->rJunctionDist);
			}
			while(G->read[breadID]->flag != JUNCTION){
				internb = G->read[breadID]->first;
				while(internb){
					if(G->read[internb->ID]->flag == CONTAINED){
						// Read is not contained in the junction read, but in a proper overlap
						pc_anno = (struct pc_anno*)reads[internb->ID].annotation;
						if(!pc_anno) printf("No Annotation pointer allocated\n");
						pc_anno->pathID = pathsNum;
						pc_anno->lJunctionDist = overhang;
						counterbread = G->read[internb->ID]->first;
						while(counterbread){
							if(counterbread->ID == breadID){
								pc_anno->rJunctionDist = pathlength - (counteroverhang + counterbread->overhang);
								break;
							}
							counterbread = counterbread->next;
						}
						if(verbose && pathsNum == 51 && overhang < 300)
						printf("Set annotation c for : %i (path: %i / distlJ: %i / distrJ: %i) -> Side: %i\n",internb->ID,pc_anno->pathID,pc_anno->lJunctionDist,pc_anno->rJunctionDist,internb->sideflag);
					}
					internb = internb->next;
				}
				internb = G->read[breadID]->first;
				while(internb){
					if(internb->sideflag == bdir && G->read[internb->ID]->flag == PROPER){
						counterbread = G->read[internb->ID]->first;
						while(counterbread){
							if(counterbread->ID == breadID){
								counteroverhang += counterbread->overhang;
								break;
							}
							counterbread = counterbread->next;
						}
		    			overhang += internb->overhang;
						pc_anno = (struct pc_anno*)reads[internb->ID].annotation;
						pc_anno->pathID = pathsNum;
						pc_anno->lJunctionDist = overhang;
						pc_anno->rJunctionDist = pathlength - counteroverhang;
						if(verbose && pathsNum == 51 && overhang < 300)
						printf("Set annotation p for : %i (path: %i / distlJ: %i / distrJ: %i) -> Side: %i\n",internb->ID,pc_anno->pathID,pc_anno->lJunctionDist,pc_anno->rJunctionDist,internb->sideflag);
						breadID = internb->ID;
						break;
					}
					if(internb->sideflag == bdir && G->read[internb->ID]->flag == JUNCTION){
						breadID = internb->ID;
					}
					internb = internb->next;
				}
			}
			if(G->read[breadID]->flag == JUNCTION){
				j_anno = (struct j_anno*)reads[breadID].annotation;
				if(verbose) printf("End Junction: %i (flag: %i)\n",breadID,j_anno->vFlag);
				paths[pathsNum].ID = pathsNum;
				paths[pathsNum].component = componentNum;
				paths[pathsNum].leftJunction = i;
				paths[pathsNum].leftPath = NULL;
				paths[pathsNum].rightJunction = breadID;
				paths[pathsNum].rightPath = NULL;
				paths[pathsNum].len = pathlength;
				paths[pathsNum].flag = 0;
				paths[pathsNum].pathdir = bdir;
				pathsNum++;
				if(pathsNum==maxPathNum){
					maxPathNum *= 2;
					struct path* temp = (struct path*)realloc(paths,sizeof(struct path)*maxPathNum);
					if(temp) paths = temp;
					else{
						printf("Path realloc failed. PROCESS KILLED\n");
						exit(1);
					}
				}
				// End junction was not fully visited, call DFS touring on this junction
				if(!j_anno->vFlag){
					if(verbose) printf("Recursive DFS call on read: %i\n",breadID);
					junctionDFS(breadID,G,reads);
				}
			}
		}
		if(verbose) printf("-> End of recursive PATH (%i -> %i)\n",i,breadID);
		if(bread && G->read[breadID]->flag == CONTAINED){
			if(verbose) printf("Read %i is contained in junction read %i\n",breadID,i);
			pc_anno = (struct pc_anno*)reads[breadID].annotation;
			pc_anno->pathID = 0;
			if(verbose) printf("Read fin\n");
		}
		bread = bread->next;
	}
	if(verbose) printf("\t\t\t\tRead: %i tagged as fully visited (RecDepth: %i)\n",i,recursionDepth);
	recursionDepth--;
}

/**
 * Reports all paths with its associated components
 */
void statsScaff(){
	int i;
	printf("Paths:\n");
	for(i=1;i<pathsNum;i++){
		printf("%i (comp: %i): J: %i -> %i / len: %i\n",i,paths[i].component,paths[i].leftJunction,paths[i].rightJunction,paths[i].len);
	}
}

/**
 * External called init function for scaffolding
 * @param G			The Overlap Graph (String Graph).
 * @param reads		Annotated reads.
 * @param numread	Number of all reads of all kinds of Libraries
 */
void initScaff(struct myovlList* G, struct reads* reads){
	char verbose = 0;
	int i;
	componentNum = 1;

    initAnnotation(G, reads);

    struct j_anno* j_anno;
//    paths = (struct path*)malloc(sizeof(struct path)*maxPathNum);

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		j_anno = (struct j_anno*)reads[i].annotation;
    		if(!j_anno->vFlag){
    			if(verbose) printf("New component found, start recursive DFS call\n");
    			junctionDFS(i,G,reads);
    			componentNum++;
    		}
    	}
    }

    if(verbose) statsScaff();

    struct bread* bread;

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		bread = G->read[i]->first;
    		while(bread){
    			if(bread->dest && bread->dest->flag) bread->dest->flag = 0;
    			bread = bread->next;
    		}
    	}
    }

    printf("Number of components: %i\n",componentNum-1);
    printf("Initial Read-Annotation finished\n");
}

/**
 * If Scaffolding was not switched on, use contigs as Scaffolds
 * @param G
 * @param reads
 * @return
 */
struct scaffold_set* contigs_init(struct myovlList* G){ //, struct reads* reads
	int i;

    struct scaffold_set* aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
    aS->num = pathsNum - 1;
    aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*pathsNum);

    for(i = 1; i < pathsNum; i++){
    	aS->scaff[i-1].ID = paths[i].ID;
    	aS->scaff[i-1].startJunction = paths[i].leftJunction;
    	aS->scaff[i-1].endJunction = paths[i].rightJunction;
    	aS->scaff[i-1].len = paths[i].len;
    	aS->scaff[i-1].type = 1; // Singletons
    	aS->scaff[i-1].first = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
    	aS->scaff[i-1].first->ID = paths[i].ID;
    	aS->scaff[i-1].first->targetJunction = paths[i].rightJunction;
    	aS->scaff[i-1].first->next = NULL;
    	aS->scaff[i-1].first->bridge = NULL;
    	aS->scaff[i-1].next = -1;

    }

//    statsScaff();

    struct bread* bread;

    // Reset flag to 0
    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		bread = G->read[i]->first;
    		while(bread){
    			if(bread->dest && bread->dest->flag) bread->dest->flag = 0;
    			bread = bread->next;
    		}
    	}
    }

    printf("Number of components: %i\n",componentNum-1);
//    printf("Contigs treated as Scaffold Singletons\n");
    return aS;
}

void free_schaffoldSet(struct scaffold_set* aS){
	struct scaffEdge* edge;
    for(int i = 0; i < pathsNum-1; i++){
    	while(aS->scaff[i].first){
    		edge = aS->scaff[i].first->next;
    		if(aS->scaff[i].first->bridge) free(aS->scaff[i].first->bridge);
    		free(aS->scaff[i].first);
    		aS->scaff[i].first = edge;
    	}
    }

	free(aS->scaff);
	free(aS);
}


/**
 * Same Dot output like readStringGraph, but based on the Overlap Graph
 * @param G			The Overlap Graph
 * @param reads		Annotated Reads of all Libraries
 * @param dotfile	Output dot file
 */
void scaffGraphDot(struct myovlList* G, struct reads* reads, char* dotfile){
	printf("CHECKPOINT: Write ScaffGraph.dot\n");
	int verbose = 0;
	FILE* dot = fopen(dotfile,"w");
	int i;
	int dir;
	int breadID;
	int breadDestID;
//	int counterDir;
	int pathID;
	int pathlen;
//	int junctionID;
	char found;
	struct bread* bread;
	struct bread* counterbread;
	struct j_anno* j_anno;
	struct pc_anno* pc_anno;

	fprintf(dot,"digraph scaff {\n");

	for(i=1;i<=G->V;i++){
		if(G->read[i] && G->read[i]->flag == JUNCTION){
			fprintf(dot,"%i [label=\"%i\" ]\n",i,i);
			j_anno = (struct j_anno*)reads[i].annotation;
			j_anno->vFlag = 0;
			dir = G->read[i]->dir;
			bread = G->read[i]->first;
			while(bread){
				if(bread->dest && !bread->dest->flag){
					bread->dest->flag = 1;
					breadDestID = bread->dest->ID;
					breadID = bread->ID;
					if(verbose) printf("i: %i -> BreadDestID %i first: %i\n",i,breadDestID,breadID);
					// bread is PROPER
					if(G->read[breadID]->flag == PROPER){
						if(verbose) printf("PROPER\n");
						if(verbose) printf("breadID: %i\n",breadID);
						pc_anno = (struct pc_anno*)reads[breadID].annotation;
						pathID = pc_anno->pathID;
						if(verbose) printf("PathID: %i\n",pathID);
						pathlen = paths[pathID].len;
					}
					// bread is a JUNCTION
					else{
						pathID = bread->dest->pathID * -1;
						pathlen = paths[pathID*-1].len;
					}
					found = 0;
					// ______
					counterbread = bread->dest->counterbread;
//					printf("Found: PID vs. CPID: %i == %i ? \n",pathID,counterbread->dest->pathID);
					counterbread->dest->flag = 1;
					found = 1;
					// ______
//					counterbread = G->read[breadDestID]->first;
//					while(counterbread){
//						printf("PID vs. CPID: %i == %i ? \n",pathID,counterbread->dest->pathID);
//						if(counterbread->dest && counterbread->dest->ID == i && bread != counterbread && abs(pathID) == counterbread->dest->pathID){
//							printf("Found: PID vs. CPID: %i == %i ? \n",pathID,counterbread->dest->pathID);
//							counterbread->dest->flag = 1;
//							found = 1;
//							break;
//						}
//						counterbread = counterbread->next;
//					}
					if(found){
						if(pathID > 0){
							fprintf(dot,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color=\"darkgreen\",label=\"%i (%i)\" ]\n",
									i,
									breadDestID,
									dir!=bread->sideflag ? "normal" : "inv",
									G->read[breadDestID]->dir != counterbread->sideflag ? "normal" : "inv",
									pathID,
									pathlen);
						}
						else if(pathID < 0){
							pathID *= -1;
							fprintf(dot,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color= \"gray\",label=\"%i (%i)\"]\n",
									i,
									breadDestID,
									dir!=bread->sideflag ? "normal" : "inv",
									G->read[breadDestID]->dir != counterbread->sideflag ? "normal" : "inv",
									pathID,
									pathlen);
						}

					}
					else{
						printf("Error at: %i -> %i (path: %i)\n",i,breadDestID,pathID);
						printf("Shit Happens\nAbort!!!\n");
						exit(1);
					}

				}
				bread = bread->next;
			}
		}
	}

	// set bridges as dashed paths
	printf("CHECKPOINT: Write ScaffGraph.dot (bridges)\n");
	verbose = 0;
	struct pathEdge* pathEdge;
	int r1junc,r2junc;
	int r1rev, r2rev;
	for(i=1;i<pathsNum;i++){
		if(paths[i].leftPath && paths[i].leftPath->estLen != -1 && paths[i].leftPath->ID > i){
			if(verbose) printf("left Test (%i)\n",i);
			pathEdge = paths[i].leftPath;
			r1junc = paths[i].leftJunction;
			if(pathEdge->junctionCon & 4) r2junc = paths[pathEdge->ID].rightJunction;
			else r2junc = paths[pathEdge->ID].leftJunction;
			r1rev = pathEdge->junctionCon & 2;
			r2rev = pathEdge->junctionCon & 1;
			if(verbose)
			printf("%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color= \"red\",style=dashed, label=\"~%i (c:%i)\"]\n",
					r1junc,
					r2junc,
					r1rev ? "normal" : "inv",
					r2rev ? "normal" : "inv",
					pathEdge->estLen,
					pathEdge->counter);
			fprintf(dot,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color= \"red\",style=dashed, label=\"~%i (c:%i)\"]\n",
					r1junc,
					r2junc,
					r1rev ? "normal" : "inv",
					r2rev ? "normal" : "inv",
					pathEdge->estLen,
					pathEdge->counter);
		}
		if(paths[i].rightPath && paths[i].rightPath->estLen != -1 && paths[i].rightPath->ID > i){
			if(verbose) printf("right Test (%i)\n",i);
			pathEdge = paths[i].rightPath;
			r1junc = paths[i].rightJunction;
			if(pathEdge->junctionCon & 4) r2junc = paths[pathEdge->ID].rightJunction;
			else r2junc = paths[pathEdge->ID].leftJunction;
			r1rev = pathEdge->junctionCon & 2;
			r2rev = pathEdge->junctionCon & 1;
			if(verbose) printf("%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color= \"red\",style=dashed, label=\"~%i (c:%i)\"]\n",
					r1junc,
					r2junc,
					r1rev ? "normal" : "inv",
					r2rev ? "normal" : "inv",
					pathEdge->estLen,
					pathEdge->counter);
			fprintf(dot,"%i -> %i [dir=both, arrowtail=%s, arrowhead=%s, color= \"red\",style=dashed, label=\"~%i (c:%i)\"]\n",
					r1junc,
					r2junc,
					r1rev ? "normal" : "inv",
					r2rev ? "normal" : "inv",
					pathEdge->estLen,
					pathEdge->counter);
		}
	}

	int scaffPaths = 0;

	if(scaffPaths){
		struct pathEdge* counterpathEdge;

		int startJunction;
		int endJunction;
		int targetEdge;

		int lbool = 0;

		for(i=1;i<pathsNum;i++){
			if(paths[i].leftPath){
				pathEdge = paths[i].leftPath;
				while(pathEdge){
					targetEdge = pathEdge->ID;

					counterpathEdge = paths[targetEdge].leftPath;
					while(counterpathEdge){
						if(counterpathEdge->ID == i) break;
						counterpathEdge = counterpathEdge->next;
					}
					if(counterpathEdge) lbool = 1;
					else{
						counterpathEdge = paths[targetEdge].rightPath;
						while(counterpathEdge){
							if(counterpathEdge->ID == i) break;
							counterpathEdge = counterpathEdge->next;
						}
						if(counterpathEdge) lbool = 0;
						else lbool = -1;
					}
					if(lbool == 1)			endJunction = paths[targetEdge].rightJunction;
					else if(lbool == 0) 	endJunction = paths[targetEdge].leftJunction;
					else{
						printf("No reverse Edge found, but have to be there\n Abort!\n");
						exit(1);
					}
					startJunction = paths[i].rightJunction;
					if(startJunction < endJunction){
						fprintf(dot,"%i -> %i [dir=both, arrowtail=normal, arrowhead=normal,color= \"gray\",style=dashed,label=\"c:%i\"]\n",startJunction,endJunction,pathEdge->counter);
					}
					pathEdge = pathEdge->next;
				}
			}
			if(paths[i].rightPath){
				pathEdge = paths[i].rightPath;
				while(pathEdge){
					targetEdge = pathEdge->ID;

					counterpathEdge = paths[targetEdge].leftPath;
					while(counterpathEdge){
						if(counterpathEdge->ID == i) break;
						counterpathEdge = counterpathEdge->next;
					}
					if(counterpathEdge) lbool = 1;
					else{
						counterpathEdge = paths[targetEdge].rightPath;
						while(counterpathEdge){
							if(counterpathEdge->ID == i) break;
							counterpathEdge = counterpathEdge->next;
						}
						if(counterpathEdge) lbool = 0;
						else lbool = -1;
					}
					if(lbool == 1)			endJunction = paths[targetEdge].rightJunction;
					else if(lbool == 0) 	endJunction = paths[targetEdge].leftJunction;
					else{
						printf("No reverse Edge found, but have to be there\n Abort!\n");
						exit(1);
					}
					startJunction = paths[i].leftJunction;
					if(startJunction < endJunction){
						fprintf(dot,"%i -> %i [dir=both, arrowtail=normal, arrowhead=normal,color= \"gray\",style=dashed,label=\"c:%i\"]\n",startJunction,endJunction,pathEdge->counter);
					}
					pathEdge = pathEdge->next;
				}
			}
		}
	}

	fprintf(dot,"}");
	fclose(dot);

    for(i=1; i < G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		bread = G->read[i]->first;
    		while(bread){
    			if(bread->dest && bread->dest->flag) bread->dest->flag = 0;
    			bread = bread->next;
    		}
    		j_anno = (struct j_anno*)reads[i].annotation;
    		j_anno->vFlag = 0;
    	}
    }

}

// pathstack (LiFo struct) stores visited paths and their end distance relatively from the origin read
int maxp = 100;		// Maximum elements on stack
int p=0; 			// temporary elements on stack
//int** pathstack = NULL;
#define PE 3		// Definition of the pathstack fields for each path
int* pathstack = NULL;
int* finpathStack = NULL;
int findist;
int finp;

void printPathStack(){
//	int i;
//	printf("Pathstack:\n");
//	for(i=0;i<p;i++){
//		printf("Path: %i -> Dist: %i -> TargetJunction: %i\n",pathstack[i*PE],pathstack[i*PE+1],pathstack[i*PE+2]);
//	}
}

void printfinPathStack(){
	int i;
	printf("FinPathstack:\n");
	for(i=0;i<finp;i++){
		printf("Path: %i -> Dist: %i -> TargetJunction: %i\n",finpathStack[i*PE],finpathStack[i*PE+1],finpathStack[i*PE+2]);
	}
}

void reportSib(struct pathEdge* pathEdge, int depth, int verbose){
	if(verbose) printf("Sibling Found\n");
	while(pathEdge){
		if(pathEdge->sibl){
			struct pathEdge* siblPath = pathEdge->sibl;
			while(siblPath){
				reportSib(siblPath,depth+1,verbose);
				siblPath = siblPath->sibl;
			}
		}
		int i;
		if(verbose){
			for(i=0;i<depth;i++){
				printf("\t");
			}
			printf("-> Path: %i, Depth: %i, Target: %i, Counts: %i (Part of Target: %i)\n",pathEdge->ID,pathEdge->depth,pathEdge->targetJunction,pathEdge->counter, pathEdge->targetcounter);
		}
		pathEdge = pathEdge->next;
	}
}

void connectPathStats(){
	int verbose = 0;
	int i;
	struct pathEdge* pathEdge;
	for(i=1;i<pathsNum;i++){
		pathEdge = paths[i].leftPath;
		if(pathEdge){
			if(verbose){
				if(paths[i].scaffflag & 32) printf("Ambi   -> ");
				if(paths[i].scaffflag & 16) printf("Proper -> ");
				else  printf("Start  -> ");

			}
			if(verbose) printf("Path %i (l: %i) -> Edges over left junction: %i:\n",i,paths[i].len,paths[i].leftJunction);
		}
		while(pathEdge){
			if(pathEdge->sibl){
				struct pathEdge* siblPath = pathEdge->sibl;
				while(siblPath){
					reportSib(siblPath,2,verbose);
					siblPath = siblPath->sibl;
				}
			}
			if(verbose){
				printf("\t -> Path: %i, Depth: %i, Target: %i, Counts: %i (Part of Target: %i)\n",pathEdge->ID,pathEdge->depth,pathEdge->targetJunction,pathEdge->counter,pathEdge->targetcounter);
			}
			pathEdge = pathEdge->next;
		}
		pathEdge = paths[i].rightPath;
		if(pathEdge){
			if(verbose){
				if(paths[i].scaffflag & 2) printf("Ambi   -> ");
				if(paths[i].scaffflag & 1) printf("Proper -> ");
				else  printf("Start  -> ");
			}
			if(verbose) printf("Path %i (l:%i) -> Edges over right junction: %i:\n",i,paths[i].len,paths[i].rightJunction);
		}
		while(pathEdge){
			if(pathEdge->sibl){
				struct pathEdge* siblPath = pathEdge->sibl;
				while(siblPath){
					reportSib(siblPath,2,verbose);
					siblPath = siblPath->sibl;
				}
			}
			if(verbose){
				printf("\t -> Path: %i, Depth: %i, Target: %i, Counts: %i (Part of Target: %i)\n",pathEdge->ID,pathEdge->depth,pathEdge->targetJunction,pathEdge->counter,pathEdge->targetcounter);
			}
			pathEdge =pathEdge->next;
		}
	}
}

void annotatePaths(){
	int i;
	struct pathEdge* edge;
	for(i=1;i<pathsNum;i++){
		// Leftside
		paths[i].scaffflag = 0;
		edge = paths[i].leftPath;
		if(edge){
			if(edge->sibl) paths[i].scaffflag = 32;
			else paths[i].scaffflag = 16;

		}
		// rightside
		edge = paths[i].rightPath;
		if(edge){
			if(edge->sibl) paths[i].scaffflag |= 2;
			else paths[i].scaffflag |= 1;
		}
	}
}

// Concordance test for reads on different paths
static char dirPathConcordance(struct myovlList* G, int rID){
	char verbose = 0;
	int j = rID -1;
	char r1 = -1;
	char r2 = -1;

	if(verbose) printf("j: %i j+1: %i\n",rID-1,rID);
	printPathStack();

	if(paths[pathstack[0]].leftJunction == pathstack[2]) r1 = 1;
	else if(paths[pathstack[0]].rightJunction == pathstack[2])	r1 = 0;
	else{
		printf("No junction found for R1: Abort\n");
		exit(1);
	}

	if(paths[pathstack[PE*(p-1)]].leftJunction == pathstack[PE*(p-1)+2]) r2 = 1;
	else if(paths[pathstack[PE*(p-1)]].rightJunction == pathstack[PE*(p-1)+2])	r2 = 0;
	else{
		printf("No junction found for R2: Abort\n");
		exit(1);
	}
	if(paths[pathstack[0]].pathdir) r1 = !r1;
	if(paths[pathstack[PE*(p-1)]].pathdir) r2 = !r2;

	if(!G->read[j]->dir) r1 = !r1;
	if(!G->read[j+1]->dir) r2 = !r2;

	if(G->read[j]->flag == CONTAINED) r1 = !r1;
	if(G->read[j+1]->flag == CONTAINED) r2 = !r2;

	if(verbose) printf("R1: %i R2: %i\n",r1,r2);
	if(r1 == -1 || r2 == -1) return -1;
	else if(!r1 && r2) return FR;
	else if(r1 && !r2) return RF;
	else if(r1 && r2) return FF;
	else return RR;
}

void connectPathsDFS2(struct myovlList* G, struct readFiles* files, struct reads* reads, int pathR2, int newjunctionID, int side, int rID,int expOri){
	if(p >= 5) return;
	int verbose = 0;
//	if(p>=1) verbose = 0;
	int verbose2 = 0;
	if(p>=10) verbose2 = 0;
	// if path found return 1 and stop all further recursive function calls
	if(verbose2) printf("Call with p: %i -> (path: %i , dist: %i) junction: %i side: %i\n",p-1,pathstack[(p-1)*PE],pathstack[(p-1)*PE+1],newjunctionID,side);
	struct bread* bread = G->read[newjunctionID]->first;
	while(bread){
		if(bread->dest && bread->sideflag == side){
			int newPath = bread->dest->pathID;
			int newSide = -1;
			if(verbose) printf("Search in breaddest: %i (path: %i)\n",bread->dest->ID,bread->dest->pathID);
			struct bread* tempbread = G->read[bread->dest->ID]->first;
			while(tempbread){
				if(tempbread->dest){
//					if(verbose) printf("Search in tempreadJunction: %i\n",tempbread->dest->ID);
				}
				if(tempbread->dest && tempbread->dest->pathID == newPath){
					newSide = tempbread->sideflag;
//					if(verbose) printf("Newside = %i\n",newSide);
					break;
				}
				tempbread = tempbread->next;
			}
			if(newSide > -1){
				pathstack[p*PE] = newPath;
				pathstack[p*PE+1] = paths[newPath].len + pathstack[(p-1)*PE+1];
				pathstack[p*PE+2] = bread->dest->ID;
				if(verbose) printf("p=%i -> Pop on: p:%i d:%i j:%i\n",p,pathstack[p*PE],pathstack[p*PE+1],pathstack[p*PE+2]);
				p++;

				if(p==maxp){
					maxp *= 2;
					pathstack = (int*)realloc(pathstack,sizeof(int)*maxp*PE);
					if(!pathstack){
						printf("Realloc Error in pathstack! Abort\n");
						exit(1);
					}
					finpathStack = (int*)realloc(finpathStack,sizeof(int)*maxp*PE);
					if(!finpathStack){
						printf("Realloc Error in finpathStack!\nAbort\n");
						exit(1);
					}
				}
				if(pathstack[(p-1)*PE] == pathR2){
					int dist;
					struct pc_anno* pc_anno = (struct pc_anno*)reads[rID].annotation;
					if(paths[pathR2].leftJunction == bread->dest->ID){
						if(verbose) printf("Junction: %i dist: %i (%i + %i )\n",bread->dest->ID,pc_anno->rJunctionDist,pathstack[(p-2)*PE+1], pc_anno->rJunctionDist);
						dist = pathstack[(p-2)*PE+1] + pc_anno->rJunctionDist;
					}
					else if(paths[pathR2].rightJunction == bread->dest->ID){
						if(verbose) printf("Junction: %i dist: %i (%i + %i)\n",bread->dest->ID,pc_anno->lJunctionDist,pathstack[(p-2)*PE+1], pc_anno->lJunctionDist);
						dist = pathstack[(p-2)*PE+1] + pc_anno->lJunctionDist;
					}
					else{
						printf("This case should never happens: Abort\n");
						exit(1);
					}
					if(verbose) printf("----------------------------------> DIST: %i\n",dist);
					if(dist > files->minInsert && dist < files->maxInsert){
						int dirCon = dirPathConcordance(G,rID);
						if(dirCon == expOri){
							// Make copy of the temporary stack and continue searching. If a second tour provides a concordant match, return and report no concordant match
							if(findist == 0){
								// ToDo: Test for correct orientations!
								// Call static function !
								if(verbose) printf("SUCCSESS -> CONCORDANCE ;-) \n");
								if(verbose) printf("--> Memcpy <-- Findist = %i\n",dist);
								if(verbose)	printf("Right Direction: %s, expected: %s\n",peOri[dirCon],peOri[files->oriPE]);
								memcpy(finpathStack,pathstack,sizeof(int)*PE*p);
								finp = p;
								findist = dist + reads[rID + 1].len;
								// TODO return is critical, hits could be ambiguous
	//							return;
								p--;
							}
							else{
								if(verbose) printf("--->> No unique concordant solution. Break search!!! <<---\n");
								findist = -1;
								return;
							}
						}
						else{
							if(verbose) printf("Wrong Direction: %s, expected: %s\n",peOri[dirCon],peOri[files->oriPE]);
							if(verbose)	printf("SUCCSESS -> But DISCORDANT Orientation :-(  \n");
							p--;
							return;
						}

					}
					else{
						if(verbose)
							printf("SUCCSESS -> But DISCORDANT Distance :-(  \n");
						p--;
						return;
					}
				}
				else{
					if(paths[newPath].len + pathstack[(p-1)*PE+1] < files->maxInsert){
						// Casus is 0 if no match, 1 if Concordant match, 2 if Discordant;
						connectPathsDFS2(G,files,reads,pathR2,bread->dest->ID,newSide ? 0 : 1, rID, expOri);
						p--;
						if(findist != 0) return;
					}
					else{
						p--;
						if(verbose2) printf("----> Distance (%i) is already to high. No further depth search\n",paths[newPath].len + pathstack[(p-1)*PE+1]);
					}
				}
			}
			else{
				printf("Case should never happens (newSide: %i): Abort\n",newSide);
				exit(1);
			}
		}
		bread = bread->next;
//		if(verbose) printf("Next bread\n");
	}
//	if(verbose) printf("End of recursion\n");
	return;
}


int connectPathsInit(struct myovlList* G, struct readFiles* files, struct reads* reads, int pathR1, int pathR2, struct pc_anno* pc_annoR1, int rID, int expOri){
	int verbose = 0;
	if(!pathstack){
		pathstack = (int *)malloc(sizeof(int)*maxp*PE); // Maximum depth 100 at this position.
		finpathStack = (int *)malloc(sizeof(int)*maxp*PE); // Maximum depth 100 at this position.
		// First Element is the Path ID, second the total distance from read R1 to the end of the currently visited path
	}

	struct pc_anno* intern_pc_anno;
//	struct j_anno* intern_j_anno;
	findist = 0;

	// go right
	// go right if the distance to the right junction is smaller than the maximum insert size
	// out source to new function
	p=0;
	if(verbose) printf("Search on the right side\n");
	if(pc_annoR1->rJunctionDist < files->maxInsert){
		int orgSide = -1;
		int rJunction = paths[pathR1].rightJunction;
		struct bread* orgbread;
		orgbread = G->read[rJunction]->first;
		if(verbose) printf("Right Junction : %i\n",rJunction);
		while(orgbread){
			if(orgbread->dest){
				if(verbose) printf("Search in read: %i\n",orgbread->dest->ID);
				if(G->read[orgbread->ID]->flag != JUNCTION){
					intern_pc_anno = (struct pc_anno*)reads[orgbread->ID].annotation;
					if(intern_pc_anno->pathID == pathR1){
						if(verbose) printf("Right Path found\n");
						orgSide = orgbread->sideflag;
						break;
					}
				}
				else{
					// Is never the path origin
				}
			}
			orgbread = orgbread->next;
		}
		if(orgSide>=0){
			pathstack[p*PE] = pathR1;
			pathstack[p*PE+1] = pc_annoR1->rJunctionDist;
			pathstack[p*PE+2] = rJunction;
			p++;
			findist = 0;
			connectPathsDFS2(G,files,reads,pathR2,rJunction,orgSide ? 0 : 1,rID,expOri);
			if(findist == -1) return 0;
		}
	}
	// TODO: Maybe there is a discordant match in the one direction, but a concordant match in the other one. Therefore I have always to calculate both directions
	// On the other hand, only concordant matches are of interest -> overwrite the pathstack if the match was discordant

	// go left if the distance to the left junction is smaller than the maximum insert size
	p=0;
	if(verbose) printf("Search on the left side\n");
	if(pc_annoR1->lJunctionDist < files->maxInsert){
		int orgSide = -1;
		int lJunction = paths[pathR1].leftJunction;
		struct bread* orgbread;
		orgbread = G->read[lJunction]->first;
		if(verbose) printf("Left Junction : %i\n",lJunction);
		while(orgbread){
			if(verbose) printf("Search in read: %i\n",orgbread->ID);
			if(orgbread->dest){
				if(G->read[orgbread->ID]->flag != JUNCTION){
					intern_pc_anno = (struct pc_anno*)reads[orgbread->ID].annotation;
					if(intern_pc_anno->pathID == pathR1){
						if(verbose) printf("left Path found\n");
						orgSide = orgbread->sideflag;
						break;
					}
				}
				else{
					// Is never the path origin
				}
			}
			orgbread = orgbread->next;
		}
		if(orgSide>=0){
			pathstack[p*PE] = pathR1;
			pathstack[p*PE+1] = pc_annoR1->lJunctionDist;
			pathstack[p*PE+2] = lJunction;
			p++;
			connectPathsDFS2(G,files,reads,pathR2,lJunction,orgSide ? 0 : 1,rID,expOri);
		}

	}

//	if(findist != 0){
//		free(finpathStack);
//	}
	// Nothing found in both directions: Interpreted as discordance;
	if(findist == 0) return 0; 		// No solutions
	else if(findist == -1) return 0; 	// Multiple solutions
	else return findist; 				// Exact concordant match
}


void connectPathsPairwise3(){
	int verbose = 0;
	int verbosenew = 0;
	int printFinPath = 0;
	int k,l;
//	int depth;
	struct pathEdge* pathEdge;
	struct pathEdge* lastEdge;

	if(printFinPath){
		printf("FinPath:\n");
		for(k=0;k<finp;k++){
			printf("P: %i -> dist: %i -> targetjunction: %i\n",finpathStack[k*PE],finpathStack[k*PE+1],finpathStack[k*PE+2]);
		}
	}

	if(verbose) printf("\t\tPaths\n");

	// Forward direction read R1 -> R2
	for(l=0;l<finp-1;l++){
		if(verbose) printf("\t\tPaths (l= %i)\n",l);
		k = l+1;
		if(finpathStack[l*PE+2] == paths[finpathStack[l*PE]].leftJunction){
			// i == 1, k == l+1
			// no first element -> create
			if(!paths[finpathStack[l*PE]].leftPath){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				pathEdge->targetJunction = finpathStack[k*PE+2];
				pathEdge->depth = k-l; // == 1
				pathEdge->counter = 1;
				if(k-l == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->sibl = NULL;
				pathEdge->next = NULL;
				pathEdge->junctionCon = -1;
				paths[finpathStack[l*PE]].leftPath = pathEdge;
				if(verbosenew) printf("NEW PATH 1: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
				// pathEdge set to first element
			}
			else{
				pathEdge = paths[finpathStack[l*PE]].leftPath;
				// First element found and counted up
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(k-l == finp-1) pathEdge->targetcounter++;
				}
				// First element not found, search in siblings
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(k-l == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						pathEdge->targetJunction = finpathStack[k*PE+2];
						pathEdge->counter = 1;
						pathEdge->depth = k-l; // == 1
						if(k-l == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						lastEdge->sibl = pathEdge;
						pathEdge->junctionCon = -1;
						if(verbosenew) printf("NEW PATH 2: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
		}
		else if(finpathStack[l*PE+2] == paths[finpathStack[l*PE]].rightJunction){
			// no first element -> create
			if(!paths[finpathStack[l*PE]].rightPath){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				pathEdge->targetJunction = finpathStack[k*PE+2];
				pathEdge->counter = 1;
				pathEdge->depth = k-l; // == 1
				if(k-l == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->next = NULL;
				pathEdge->sibl = NULL;
				pathEdge->junctionCon = -1;
				paths[finpathStack[l*PE]].rightPath = pathEdge;
				if(verbosenew) printf("NEW PATH 3: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
				// pathEdge set to first element
			}
			else{
				pathEdge = paths[finpathStack[l*PE]].rightPath;
				// First element found and counted up
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(k-l == finp-1) pathEdge->targetcounter++;
				}
				// First element not found, search in siblings
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(k-l == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						pathEdge->targetJunction = finpathStack[k*PE+2];
						pathEdge->counter = 1;
						pathEdge->depth = k-l; // == 1
						if(k-l == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						lastEdge->sibl = pathEdge;
						pathEdge->junctionCon = -1;
						if(verbosenew) printf("NEW PATH 4: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
		}
		else{
			printfinPathStack();
			printf("3. Both junctions are not registered as end junctions of this path\nAbort\n");
			exit(1);
		}
		k++;
		if(verbose) printf("-> Path3: First Link: %i -> %i (count: %i)\n",finpathStack[l*PE],pathEdge->ID,pathEdge->counter);
		while(k<finp){
			lastEdge = pathEdge;
			if(!pathEdge->next){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				pathEdge->targetJunction = finpathStack[k*PE+2];
				pathEdge->counter = 1;
				pathEdge->depth = k-l;
				if(k-l == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->next = NULL;
				pathEdge->sibl = NULL;
				pathEdge->junctionCon = -1;
				lastEdge->next = pathEdge;
				if(verbosenew) printf("NEW PATH 5: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
			}
			else{
				pathEdge = pathEdge->next;
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(k-l == finp-1) pathEdge->targetcounter++;
				}
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(k-l == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						pathEdge->targetJunction = finpathStack[k*PE+2];
						pathEdge->counter = 1;
						pathEdge->depth = k-l;
						if(k-l == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						pathEdge->junctionCon = -1;
						lastEdge->sibl = pathEdge;
						if(verbosenew) printf("NEW PATH 6: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
			if(verbose) printf(">%i / %i (path: %i -> %i, c: %i, depth: %i)\n",k,finp,finpathStack[l*PE],pathEdge->ID,pathEdge->counter,pathEdge->depth);
			k++;
		}

	}
	// Backward direction R2 -> R1
	if(verbose) printf("Tag Reverse direction\n");
//	verbosenew = 1;
	for(l=finp-1;l>0;l--){
		k = l-1;
		if(finpathStack[(l-1)*PE+2] == paths[finpathStack[l*PE]].leftJunction){
			// i == 1, k == l+1
			// no first element -> create
			if(!paths[finpathStack[l*PE]].leftPath){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				//_________________________
//				pathEdge->targetJunction = finpathStack[k*PE+2];
				if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
				else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
				//_________________________
				pathEdge->counter = 1;
				pathEdge->depth = l-k; // == 1
				if(l-k == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->sibl = NULL;
				pathEdge->next = NULL;
				pathEdge->junctionCon = -1;
				paths[finpathStack[l*PE]].leftPath = pathEdge;
				if(verbosenew) printf("NEW PATH 7: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
				// pathEdge set to first element
			}
			else{
				pathEdge = paths[finpathStack[l*PE]].leftPath;
				// First element found and counted up
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(l-k == finp-1) pathEdge->targetcounter++;
				}
				// First element not found, search in siblings
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(l-k == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
						else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
						pathEdge->counter = 1;
						pathEdge->depth = l-k; // == 1
						if(l-k == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						pathEdge->junctionCon = -1;
						lastEdge->sibl = pathEdge;
						if(verbosenew) printf("NEW PATH 8: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
		}
		else if(finpathStack[(l-1)*PE+2] == paths[finpathStack[l*PE]].rightJunction){
			// i == 1, k == l+1
			// no first element -> create
			if(!paths[finpathStack[l*PE]].rightPath){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				//__________________________
//				pathEdge->targetJunction = finpathStack[k*PE+2];
				if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
				else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
				//__________________________
				pathEdge->counter = 1;
				pathEdge->depth = l-k; // == 1
				if(l-k == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->sibl = NULL;
				pathEdge->next = NULL;
				pathEdge->junctionCon = -1;
				paths[finpathStack[l*PE]].rightPath = pathEdge;
				if(verbosenew) printf("NEW PATH 9: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
				// pathEdge set to first element
			}
			else{
				pathEdge = paths[finpathStack[l*PE]].rightPath;
				// First element found and counted up
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(l-k == finp-1) pathEdge->targetcounter++;
				}
				// First element not found, search in siblings
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(l-k == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
						else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
						pathEdge->counter = 1;
						pathEdge->depth = l-k; // == 1
						if(l-k == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						pathEdge->junctionCon = -1;
						lastEdge->sibl = pathEdge;
						if(verbosenew) printf("NEW PATH 10: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
		}
		else{
			printfinPathStack();
			printf("4. Both junctions are not registered as end junctions of this path\nAbort\n");
			printf("revPath: %i -> Junction: %i r: %i, l: %i\n",finpathStack[l*PE],finpathStack[(l-1)*PE+2],paths[finpathStack[l*PE]].rightJunction,paths[finpathStack[l*PE]].leftJunction);
			exit(1);
		}
		if(verbose) printf("-> rev Path3: First Link: %i -> %i (count: %i)\n",finpathStack[l*PE],pathEdge->ID,pathEdge->counter);
		k--;
		while(k>=0){
			lastEdge = pathEdge;
			if(!pathEdge->next){
				pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
				pathEdge->ID = finpathStack[k*PE];
				if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
				else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
				pathEdge->counter = 1;
				pathEdge->depth = l-k;
				if(l-k == finp-1) pathEdge->targetcounter = 1;
				else pathEdge->targetcounter = 0;
				pathEdge->estLen = -1;
				pathEdge->next = NULL;
				pathEdge->sibl = NULL;
				pathEdge->junctionCon = -1;
				lastEdge->next = pathEdge;
				if(verbosenew) printf("NEW PATH 11: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
//				printfinPathStack();
			}
			else{
				pathEdge = pathEdge->next;
				if(pathEdge->ID == finpathStack[k*PE]){
					pathEdge->counter++;
					if(l-k == finp-1) pathEdge->targetcounter++;
				}
				else{
					lastEdge = pathEdge;
					pathEdge = pathEdge->sibl;
					while(pathEdge){
						// Sibling found, count up
						if(pathEdge->ID == finpathStack[k*PE]){
							pathEdge->counter++;
							if(l-k == finp-1) pathEdge->targetcounter++;
							break;
						}
						if(!pathEdge->sibl) lastEdge = pathEdge;
						pathEdge = pathEdge->sibl;
					}
					// sibling not found, create new
					if(!pathEdge){
						pathEdge = (struct pathEdge*)malloc(sizeof(struct pathEdge));
						pathEdge->ID = finpathStack[k*PE];
						if(paths[finpathStack[k*PE]].leftJunction == finpathStack[k*PE+2]) pathEdge->targetJunction = paths[finpathStack[k*PE]].rightJunction;
						else pathEdge->targetJunction = paths[finpathStack[k*PE]].leftJunction;
						pathEdge->counter = 1;
						pathEdge->depth = l-k;
						if(l-k == finp-1) pathEdge->targetcounter = 1;
						else pathEdge->targetcounter = 0;
						pathEdge->estLen = -1;
						pathEdge->next = NULL;
						pathEdge->sibl = NULL;
						pathEdge->junctionCon = -1;
						lastEdge->sibl = pathEdge;
						if(verbosenew) printf("NEW PATH 12: \t -> %i -> %i :path: %i (depth: %i), target junction: %i\n",l,k,pathEdge->ID,pathEdge->depth,pathEdge->targetJunction);
					}
				}
			}
			if(verbose) printf(">%i / %i (path: %i -> %i, c: %i, depth: %i)\n",k,finp,finpathStack[l*PE],pathEdge->ID,pathEdge->counter,pathEdge->depth);
			k--;
		}
	}
}


inline static char readside(struct myovlList* G, int id){
	int verbose = 0;
	struct bread* sidebread;
	if(G->read[id]->flag == CONTAINED){
		sidebread = G->read[id]->first;
		if(sidebread->next){
			printf("More then one overlap: Should not be\n Abort");
			exit(1);
		}
		if(verbose) printf("Contained -> breaddir = %i\n",sidebread->dir);
		if(sidebread->dir){
			return !G->read[sidebread->ID]->dir;
		}
		else{
			return G->read[sidebread->ID]->dir;
		}
	}
	else{
		return G->read[id]->dir;
	}
	return -1;
}

char bridgeReadsSide(struct myovlList* G, struct reads* reads, int readID, char leftEnd, char proper){
	char verbose = 0;
	struct pc_anno* pc_anno = (struct pc_anno*)reads[readID].annotation;
	int leftJID;
	int rightJID;
	if(leftEnd){
		leftJID = paths[pc_anno->pathID].leftJunction;
		rightJID = paths[pc_anno->pathID].rightJunction;
	}
	else{
		rightJID = paths[pc_anno->pathID].leftJunction;
		leftJID = paths[pc_anno->pathID].rightJunction;
	}

	int lJside = -1;

	struct bread* bread = G->read[leftJID]->first;
	while(bread){
		if(bread->dest && bread->dest->ID == rightJID){
			lJside = bread->sideflag;
			break;
		}
		bread = bread->next;
	}

	if(verbose) printf("Proper?: %i\n",(int)proper);
	if(!proper) lJside = (lJside+1)%2;

	if(verbose) printf("ReadDir: %i SpurJunctionDir: %i (left?:%i)\n",G->read[readID]->dir,lJside,(int)leftEnd);

	return G->read[readID]->dir == lJside;
}

inline static char dirConcordance(int r1pathID, int r2pathID, char r1dir, char r2dir, int r1lJDis, int r2lJDis){
	static char verbose = 0;
	int r1 = -1;	// 0 if foreward 1 if reverse;
	int r2 = -1;
	if(r1pathID == r2pathID){
		// outside of lJunction is 1
		if(paths[r1pathID].pathdir){
			if(r1lJDis < r2lJDis){
				if(r1dir == 0) r1 = 0;
				else r1 = 1;
				if(r2dir == 0) r2 = 0;
				else r2 = 1;
			}
			else{
				if(r2dir == 0) r1 = 0;
				else r1 = 1;
				if(r1dir == 0) r2 = 0;
				else r2 = 1;
			}
		}
		else{
			if(r1lJDis < r2lJDis){
				if(r1dir == 0) r1 = 1;
				else r1 = 0;
				if(r2dir == 0) r2 = 1;
				else r2 = 0;
			}
			else{
				if(r2dir == 0) r1 = 1;
				else r1 = 0;
				if(r1dir == 0) r2 = 1;
				else r2 = 0;
			}
		}
	}
	if(verbose) printf("r1: %i r2: %i\n",r1,r2);
	if(r1 == -1 || r2 == -1) return -1;
	else if(r1 && r2) return 2;
	else if(!r1 && !r2) return 2;
	else if(r1 && !r2) return 1;
	else return 0;

}

char isSpur(struct myovlList* G, int pathID){
	char verbose = 0;
	char verbose2 = 0;
	char bool1 = -1;
	char bool2 = 0;
	char bool3 = 0;
	struct bread* bread;
	bread = G->read[paths[pathID].rightJunction]->first;
	if(verbose) printf("Look rigth: %i (pathID: %i)\n",paths[pathID].rightJunction, pathID);
	while(bread){
		if(verbose) printf("Search right: %i\n",bread->sideflag);
		if(bool1 == -1) bool1 = bread->sideflag;
		else if(bool1 != bread->sideflag){
			if(verbose) printf("Both sides found from the right Junction\n");
			bool2 = 1;
			break;
		}
		bread = bread->next;
	}
	if(!bool2 && verbose2) printf("RightSide (Read: %i of Path: %i) is dead (spur) end\n",paths[pathID].rightJunction,pathID);

	bool1 = -1;
	bread = G->read[paths[pathID].leftJunction]->first;
	if(verbose) printf("Look left: %i (pathID: %i)\n",paths[pathID].leftJunction, pathID);
	while(bread){
		if(verbose) printf("\tSearch left: %i\n",bread->sideflag);
		if(bool1 == -1) bool1 = bread->sideflag;
		else if(bool1 != bread->sideflag){
			if(verbose) printf("Both sides found from the left Junction\n");
			bool3 = 1;
			break;
		}
		bread = bread->next;
	}
	if(!bool3 && verbose2) printf("LeftSide (Read: %i of Path: %i) is dead (spur) end\n",paths[pathID].leftJunction,pathID);

	return (bool3 << 1) | bool2;
}

void printPath(struct myovlList* G, struct reads* reads, int pathID){
	int readID = paths[pathID].leftJunction;
	int breadID;
	printf("LeftJunction: %i\n",readID);
	struct bread* bread = G->read[readID]->first;
	struct bread* internbread;
	struct pc_anno* pc_anno;
	int sideflag;

	printf("List of Connected Junctions:\n");
	while(bread){
		if(bread->dest){
			printf(" --> %i\n",bread->dest->ID);
		}
		bread = bread->next;
	}
	printf("\n\n");

	bread = G->read[readID]->first;
	while(bread){
		if(bread->dest && bread->dest->ID == paths[pathID].rightJunction){
			sideflag = bread->sideflag;
			break;
		}
		bread = bread->next;
	}
	if(bread){
		breadID = bread->ID;
		while(breadID != paths[pathID].rightJunction){
			pc_anno = (struct pc_anno*)reads[breadID].annotation;
			printf(" P (ID: %i) --> lDist: %i -- rDist: %i\n",breadID,pc_anno->lJunctionDist, pc_anno->rJunctionDist);

			internbread = G->read[breadID]->first;
			while(internbread){
				if(G->read[internbread->ID]->flag == CONTAINED){
					pc_anno = (struct pc_anno*)reads[internbread->ID].annotation;
					printf(" C (ID: %i) --> lDist: %i -- rDist: %i\n",internbread->ID,pc_anno->lJunctionDist, pc_anno->rJunctionDist);
				}
				internbread = internbread->next;
			}
			internbread = G->read[breadID]->first;
			while(internbread){
				if(G->read[internbread->ID]->flag != CONTAINED && internbread->sideflag == sideflag){
					breadID = internbread->ID;
					break;
				}
				internbread = internbread->next;
			}
		}
		printf("RightJunction Reached: %i\n",breadID);
		printf("PathEnD\n");
		exit(1);
	}
	else{
		printf("No correct bread found\n");
		exit(1);
	}
}

//inline static char bridgeJunctionOri(struct myovlList* G, int r1J, int r2J ,char r1Right, char r2right){
//	char flag;
//	if(G->read[])
//
//	return 0;
//}

void setVirtualBridge(struct myovlList* G, int r1path, int r2path, char r1right, char r2right, int dist){
//	printf("Set Virtual Bridge\n");
	char verbose = 0;
	struct pathEdge* pathedgeR1;
	struct pathEdge* pathedgeR2;
	struct pathEdge* tempEdge;
	int avgDist;

	if(r1right)	pathedgeR1 = paths[r1path].rightPath;
	else pathedgeR1 = paths[r1path].leftPath;
	while(pathedgeR1){
		if(pathedgeR1->ID == r2path){
			// ToDo: Seek in Siblings
			// Set new Parameter
			avgDist = (pathedgeR1->counter * pathedgeR1->estLen) + dist;
			pathedgeR1->counter++;
			pathedgeR1->estLen = avgDist / pathedgeR1->counter;
			if(verbose) printf("Update Bridge: %i -> %i (dist: %i, avg: %i) Count: %i\n",r1path,r2path,dist,pathedgeR1->estLen,pathedgeR1->counter);
			break;
		}
		if(pathedgeR1->sibl){
			tempEdge = pathedgeR1->sibl;
			while(tempEdge != pathedgeR1){
				if(tempEdge->ID == r2path){
					pathedgeR1 = tempEdge;
					avgDist = (pathedgeR1->counter * pathedgeR1->estLen) + dist;
					pathedgeR1->counter++;
					pathedgeR1->estLen = avgDist / pathedgeR1->counter;
					if(verbose) printf("Update Bridge (Sibl!: BAD CASE): %i -> %i (dist: %i, avg: %i) Count: %i\n",r1path,r2path,dist,pathedgeR1->estLen,pathedgeR1->counter);
					break;
				}
				tempEdge = tempEdge->sibl;
			}
		}
		pathedgeR1 = pathedgeR1->next;
	}

	if(r2right)	pathedgeR2 = paths[r2path].rightPath;
	else pathedgeR2 = paths[r2path].leftPath;
	while(pathedgeR2){
		if(pathedgeR2->ID == r1path){
			// ToDo: Seek in Siblings
			// Set new Parameter
			avgDist = (pathedgeR2->counter * pathedgeR2->estLen) + dist;
			pathedgeR2->counter++;
			pathedgeR2->estLen = avgDist / pathedgeR2->counter;
			if(verbose)printf("Update Bridge: %i -> %i (dist: %i, avg: %i) Count: %i\n",r2path,r1path,dist,pathedgeR2->estLen,pathedgeR2->counter);
			break;
		}
		if(pathedgeR2->sibl){
			tempEdge = pathedgeR2->sibl;
			while(tempEdge != pathedgeR2){
				if(tempEdge->ID == r1path){
					pathedgeR2 = tempEdge;
					avgDist = (pathedgeR2->counter * pathedgeR2->estLen) + dist;
					pathedgeR2->counter++;
					pathedgeR2->estLen = avgDist / pathedgeR2->counter;
					if(verbose) printf("Update Bridge (Sibl!: BAD CASE): %i -> %i (dist: %i, avg: %i) Count: %i\n",r2path,r1path,dist,pathedgeR2->estLen,pathedgeR2->counter);
					break;
				}
				tempEdge = tempEdge->sibl;
			}
		}
		pathedgeR2 = pathedgeR2->next;
	}

	if(pathedgeR1){
		return;
	}

	// pathedgeR1 == NULL; -> Set new edge
	// Path not found, det new one!!!
	pathedgeR1 = (struct pathEdge*)malloc(sizeof(struct pathEdge));
	pathedgeR2 = (struct pathEdge*)malloc(sizeof(struct pathEdge));
	pathedgeR1->sibl = NULL;
	pathedgeR2->sibl = NULL;

	pathedgeR1->ID = r2path;
	pathedgeR1->counter = 1;
	pathedgeR1->depth = 1;
	pathedgeR1->estLen = dist;

	pathedgeR2->ID = r1path;
	pathedgeR2->counter = 1;
	pathedgeR2->depth = 1;
	pathedgeR2->estLen = dist;

	//TODO: Left or right junction: the target should not be the corresponding junction to be consistent with the other pathedged???

	if(r1right){
		pathedgeR2->targetJunction = paths[r1path].leftJunction;
//		pathedgeR2->targetJunction = paths[r1path].rightJunction;
		if(paths[r1path].rightPath){
			if(verbose) printf("R1 right: Multiple Virtual Edges\n");
			return;
			exit(1);
		}
		else{
			paths[r1path].rightPath = pathedgeR1;
			pathedgeR1->next = NULL;
			pathedgeR1->next = NULL;
		}
	}
	else{
		pathedgeR2->targetJunction = paths[r1path].rightJunction;
		if(paths[r1path].leftPath){
			if(verbose) printf("R1 left: Multiple Virtual Edges\n");
			return;
			exit(1);
		}
		else{
			paths[r1path].leftPath = pathedgeR1;
			pathedgeR1->next = NULL;
			pathedgeR1->next = NULL;
		}
	}

	if(r2right){
		pathedgeR1->targetJunction = paths[r2path].leftJunction;
		if(paths[r2path].rightPath){
			if(verbose) printf("R2 right: Multiple Virtual Edges\n");
			return;
			exit(1);
		}
		else{
			paths[r2path].rightPath = pathedgeR2;
			pathedgeR2->next = NULL;
			pathedgeR2->next = NULL;
		}
	}
	else{
		pathedgeR1->targetJunction = paths[r2path].rightJunction;
		if(paths[r2path].leftPath){
			if(verbose) printf("R2 left: Multiple Virtual Edges\n");
			return;
			exit(1);
		}
		else{
			paths[r2path].leftPath = pathedgeR2;
			pathedgeR2->next = NULL;
			pathedgeR2->next = NULL;
		}
	}

	// Junction Connection:
	// Junction bit: 0 if left, 1 if right
	if(r2right) pathedgeR1->junctionCon = (1 << 2);
	else pathedgeR1->junctionCon = 0;
	if(r1right) pathedgeR2->junctionCon = (1 << 2);
	else pathedgeR2->junctionCon = 0;
	// Is this correct?????
	// Ori bit incoming: 0 if forward to original read, 1 if reverse
//	printf("First Flag in r1: %i\n",G->read[G->read[paths[r1path].rightJunction]->first->ID]->flag);
//	printf("First Flag in r2: %i\n",G->read[G->read[paths[r2path].rightJunction]->first->ID]->flag);


	if(r1right)	pathedgeR1->junctionCon |= ((!G->read[paths[r1path].rightJunction]->first->sideflag) << 1);
	else pathedgeR1->junctionCon |= ((!G->read[paths[r1path].leftJunction]->first->sideflag) << 1);
	if(r2right) pathedgeR2->junctionCon |= ((!G->read[paths[r2path].rightJunction]->first->sideflag) << 1);
	else pathedgeR2->junctionCon |= ((!G->read[paths[r2path].leftJunction]->first->sideflag) << 1);
	// Ori bit outgoing 0 if forward to original read, 1 if reverse
	if(r2right) pathedgeR1->junctionCon |= (!G->read[paths[r2path].rightJunction]->first->sideflag);
	else pathedgeR1->junctionCon |= (!G->read[paths[r2path].leftJunction]->first->sideflag);
	if(r1right) pathedgeR2->junctionCon |= (!G->read[paths[r1path].rightJunction]->first->sideflag);
	else pathedgeR2->junctionCon |= (!G->read[paths[r1path].leftJunction]->first->sideflag);

	if(verbose) printf("Set New Bridge: %i -> %i (dist: %i) Count: %i (TargetJunction: %i)\n",r1path,r2path,dist,pathedgeR1->counter,pathedgeR1->targetJunction);
	if(verbose) printf("Set New Bridge: %i -> %i (dist: %i) Count: %i (TargetJunction: %i)\n",r2path,r1path,dist,pathedgeR2->counter,pathedgeR2->targetJunction);

	// Print the connection in scaffGraph.dot (with dashed lines);
	// Scaffold stats: the touring is essetial;

}

void buildBridge(struct myovlList* G, struct readFiles lib, struct reads* reads, int r1ID, char spurs, int r1path, int r2path, int oriPE){
	// Look which side is spur? Try all of bit type is 0
	char verbose = 0;
	int dist = 0;
	int hitnum = 0;
	struct pc_anno* r1anno = (struct pc_anno*)reads[r1ID].annotation;
	struct pc_anno* r2anno = (struct pc_anno*)reads[r1ID+1].annotation;
	char r1dir, r2dir; // is 1 if the direction is forward in regard of the dead-end junction, 0 otherwise
//	printf("Spur: %i\n",spurs);
	if(!((spurs >> 3) & 1)){
		// LL
		if(!((spurs >> 1) & 1)){
			if(verbose) printf("LL (DIR: %i (FLAG: %i) -> %i (FLAG: %i))\n",G->read[r1ID]->dir,G->read[r1ID]->flag, G->read[r1ID + 1]->dir,G->read[r1ID+1]->flag);
			r1dir = bridgeReadsSide(G,reads,r1ID,1,G->read[r1ID ]->flag == PROPER ? 1 : 0);
			r2dir = bridgeReadsSide(G,reads,r1ID+1,1,G->read[r1ID +1]->flag == PROPER ? 1 : 0);
			dist = r1anno->lJunctionDist + r2anno->lJunctionDist + G->read[r1ID]->length + G->read[r1ID+1]->length;
			if(dist < lib.maxInsert && ((oriPE == 0 && r1dir && r2dir) || (oriPE == 1 && !r1dir && !r2dir) || (oriPE == 2 && (r1dir != r2dir)))){
				if(verbose) printf("Set virtual Bridge\n");
				if(verbose) printf("Found Bridge over R1Left: %i R2Left: %i (Paths: %i : %i)\n",r1anno->lJunctionDist,r2anno->lJunctionDist,r1path,r2path);
				setVirtualBridge(G,r1path,r2path,0,0,dist-lib.avgInsert);
				hitnum++;
			}
		}
		// LR
		if(!(spurs & 1)){
			if(verbose) printf("LR (DIR: %i (FLAG: %i) -> %i (FLAG: %i))\n",G->read[r1ID]->dir,G->read[r1ID]->flag, G->read[r1ID + 1]->dir,G->read[r1ID+1]->flag);
			r1dir = bridgeReadsSide(G,reads,r1ID,1,G->read[r1ID ]->flag == PROPER ? 1 : 0);
			r2dir = bridgeReadsSide(G,reads,r1ID+1,0,G->read[r1ID +1]->flag == PROPER ? 1 : 0);
			dist = r1anno->lJunctionDist + r2anno->rJunctionDist + G->read[r1ID]->length + G->read[r1ID+1]->length;
			if(dist < lib.maxInsert && ((oriPE == 0 && r1dir && r2dir) || (oriPE == 1 && !r1dir && !r2dir) || (oriPE == 2 && (r1dir != r2dir)))){
				if(verbose) printf("Set virtual Bridge\n");
				if(verbose) printf("Found Bridge over R1Left: %i R2Right: %i (Paths: %i : %i)\n",r1anno->lJunctionDist,r2anno->rJunctionDist,r1path,r2path);
				setVirtualBridge(G,r1path,r2path,0,1,dist-lib.avgInsert);
				hitnum++;
			}
		}
	}
	if(!((spurs >> 2) & 1)){
		// RL
		if(!((spurs >> 1) & 1)){
			if(verbose) printf("RL (DIR: %i (FLAG: %i) -> %i (FLAG: %i))\n",G->read[r1ID]->dir,G->read[r1ID]->flag, G->read[r1ID + 1]->dir,G->read[r1ID+1]->flag);
			r1dir = bridgeReadsSide(G,reads,r1ID,0,G->read[r1ID ]->flag == PROPER ? 1 : 0);
			r2dir = bridgeReadsSide(G,reads,r1ID+1,1,G->read[r1ID +1]->flag == PROPER ? 1 : 0);
			dist = r1anno->rJunctionDist + r2anno->lJunctionDist + G->read[r1ID]->length + G->read[r1ID+1]->length;
			if(dist < lib.maxInsert && ((oriPE == 0 && r1dir && r2dir) || (oriPE == 1 && !r1dir && !r2dir) || (oriPE == 2 && (r1dir != r2dir)))){
				if(verbose) printf("Set virtual Bridge\n");
				if(verbose) printf("Found Bridge over R1Right: %i R2Left: %i (Paths: %i : %i)\n",r1anno->rJunctionDist,r2anno->lJunctionDist,r1path,r2path);
				setVirtualBridge(G,r1path,r2path,1,0,dist-lib.avgInsert);
				hitnum++;
			}
		}
		// RR
		if(!(spurs & 1)){
			r1dir = bridgeReadsSide(G,reads,r1ID,0,G->read[r1ID ]->flag == PROPER ? 1 : 0);
			r2dir = bridgeReadsSide(G,reads,r1ID+1,0,G->read[r1ID +1]->flag == PROPER ? 1 : 0);
			dist = r1anno->rJunctionDist + r2anno->rJunctionDist + G->read[r1ID]->length + G->read[r1ID+1]->length;
			if(verbose) printf("RR (DIR: %i (oriPE: %i) (FLAG: %i) -> %i (FLAG: %i)) dist<maxsize?(%i<%i)\n",G->read[r1ID]->dir,oriPE,G->read[r1ID]->flag, G->read[r1ID + 1]->dir,G->read[r1ID+1]->flag,dist,lib.maxInsert);
//			printf("%i < %i && (oriPE(%i) == 0 && r1dir(%i) && r2dir(%i))",dist,lib.maxInsert,oriPE,r1dir,r2dir);
			if(dist < lib.maxInsert && ((oriPE == 0 && r1dir && r2dir) || (oriPE == 1 && !r1dir && !r2dir) || (oriPE == 2 && (r1dir != r2dir)))){
				if(verbose) printf("Set virtual Bridge\n");
				if(verbose) printf("Found Bridge over R1Right: %i R2Right: %i (Paths: %i : %i)\n",r1anno->rJunctionDist,r2anno->rJunctionDist,r1path,r2path);
				setVirtualBridge(G,r1path,r2path,1,1,dist-lib.avgInsert);
				hitnum++;
			}
		}
	}
	if(hitnum != 1 && verbose){
		printf("Found NOTHING or to MUCH\n");
		printf("Print path 45\n");
//		if(r1path == 45) printPath(G,reads,r1path);
//		if(r2path == 45) printPath(G,reads,r2path);

	}

}

char isBridge(struct myovlList* G, int r1pathID, int r2pathID){
	char r1side = isSpur(G,r1pathID);
	char r2side = isSpur(G,r2pathID);
	if(r1side != 3 && r2side != 3){
		return (r1side << 2) | r2side;
	}
	return -1;
}

// The function will be called recursively if a path have to proper outgoing edges and so the first one is marked as having a sibling.
// Go straight in this function and call the function recursively with the found siblings
/**
 * Heart of the scaffolding module. Iteration over all read libraries and the reads. If the reads are paired the function categorizes the in order
 * of their position to each other, and stores the order of the contigs in the way the are toured to reach one read end from the other.
 * Reads may be DENIED, if they are not part of the graph anymore; SINGLE, if the are part of a SE library;
 * WIDOWED, if their corresponding read is denied; DISCORDANT, if both ends are part of the graph, but the are not in the expected distance or orientation
 * to each other; or CONCORDANT if both ends are in correct distance and orientation
 *
 * @param G			The String graph, containing all overlap informations
 * @param files		The Meta struct of all read libraries. Reads can be SE, PE, or MP libs. Information of the Min and Max distance accepted for a concordant match
 * @param reads		The Read database, including meta information of the reads: The contig and component the belongs to. The function in the graph (JUNCTION, PROPER, WIDOWED, CONTAINED)
 */
void readTouring(struct myovlList* G, struct readFiles* files, struct reads* reads){
	printf("Start ReadTouring\n");
	char verbose = 0;
	char verbose2 = 0;
	char verbose3 = 0;
	char verboseBridge = 0;
	char verboseBridgeYes = 0;
	char debug = 0;

	int i,j;
	int startID, endID;
	int pe = 0;

	// PE relation info
	int pathR1, pathR2;
	int compR1 = -1, compR2 = -1;
	int ldistR1 = -1, ldistR2 = -1;
	int rdistR1 = -1, rdistR2 = -1;

	// Concordance Stats
	int insert;
	int insert_avg;
	int insert_sum = 0;
	int insert_stddev;
	int insert_stddev_sum = 0;
	int insert_num;
	int insert_var_sum;
	float insert_var = 0;
	int* insert_array;

	// Read Annotation
	struct j_anno* j_annoR1 = NULL;
	struct j_anno* j_annoR2 = NULL;
	struct pc_anno* pc_annoR1 = NULL;
	struct pc_anno* pc_annoR2 = NULL;

	int minIns;
	int maxIns;
	int wi_num_tot = 0;			// Widowed Reads
	int di_num_tot = 0 ;		// Discordant Reads
	int co_num_tot = 0;			// Concordant Reads
	int	de_num_tot = 0;			// Rejected Reads
	int un_num_tot = 0;			// Unknown Reads
	int si_num_tot = 0;			// Single End Reads
	int dicomp_num_tot = 0;		// Reads on different components

	int interpathdist = 0;

	struct bread* bread;
// _____________________________________
	// Test for the right orientation of both reads (fr, ff or rr, rf)
	int r1side, r2side;
// _____________________________________
	int foundsum = 0;
	int nfoundsum = 0;


	for(i=0;i<files->libNum;i++){
		startID = files[i].startId;
		endID = files[i].endId;
		printf("Reads: %i -> %i\n",startID,endID);
		if(files[i].rightReads) pe = 1;
		if(pe){
			minIns = files[i].minInsert;
			maxIns = files[i].maxInsert;
			insert_array = (int*)malloc(sizeof(int)*((endID-startID)+1));
			printf("\n##### MP/PE LIB #####\n");
			insert_num = 0;
			insert_sum = 0;
			insert_stddev_sum = 0;
			int wi_num = 0;
			int di_num = 0 ;
			int mu_num = 0;
			int co_num = 0;
			int re_num = 0;
			int	de_num = 0;
			int un_num = 0;
			int br_num = 0;
			int dicomp_num = 0;		// Reads on different components
			// Annotate read Positions
			for(j=startID;j<=endID;j+=2){
//				if(j%10000==0 || j+1%10000==0)
					if(debug) printf("ReadPair: %i, %i\n",j,j+1);
				reads[j].flag = CONCORDANT;
				reads[j+1].flag = CONCORDANT;
				if(!G->read[j]){
//					if(verbose) printf("\t -> READ IS NOT PART OF THE GRAPH\n");
					reads[j].flag = REJECTED;
					pathR1 = 0;
				}
				else if(G->read[j]->flag == CONTAINED && G->read[j]->first){
					if(verbose) printf("Read %i was rejected in Stringer because of not found start and end connection\n",j);
					reads[j].flag = REJECTED;
					pathR1 = 0;
				}
				else if(G->read[j]->flag == CONTAINED || G->read[j]->flag == PROPER){
					pc_annoR1 = (struct pc_anno*)reads[j].annotation;
					pathR1 = pc_annoR1->pathID;
					compR1 = paths[pathR1].component;
					ldistR1 = pc_annoR1->lJunctionDist;
					rdistR1 = pc_annoR1->rJunctionDist;
//					if(verbose) printf("\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
				}
				else if(G->read[j]->flag == JUNCTION){
					j_annoR1 = (struct j_anno*)reads[j].annotation;
//					if(verbose) printf("\t -> ID: %i -> Read is a JUNCTION\n",j);
					pathR1 = 0;
				}
				else{
//					if(verbose) printf("\t -> READ IS NOT PART OF THE GRAPH\n");
					reads[j].flag = REJECTED;
					pathR1 = 0;
				}
				if(!G->read[j+1]){
//					if(verbose) printf("\t -> READ IS NOT PART OF THE GRAPH\n");
					reads[j+1].flag = REJECTED;
					pathR2 = 0;
				}
				else if(G->read[j+1]->flag == CONTAINED && G->read[j+1]->first){
					if(verbose) printf("Read %i was rejected in Stringer because of not found start and end connection\n",j+1);
					reads[j+1].flag = REJECTED;
					pathR2 = 0;
				}
				else if(G->read[j+1]->flag == CONTAINED || G->read[j+1]->flag == PROPER){
					pc_annoR2 = (struct pc_anno*)reads[j+1].annotation;
					pathR2 = pc_annoR2->pathID;
					compR2 = paths[pathR2].component;
					ldistR2 = pc_annoR2->lJunctionDist;
					rdistR2 = pc_annoR2->rJunctionDist;
//					if(verbose) printf("\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
				}
				else if(G->read[j+1]->flag == JUNCTION){
					j_annoR2 = (struct j_anno*)reads[j+1].annotation;
//					if(verbose) printf("\t -> ID: %i -> Read is a JUNCTION\n",j+1);
					pathR2 = 0;
				}
				else{
//					if(verbose) printf("\t -> READ IS NOT PART OF THE GRAPH (WIDOWED)\n");
					reads[j+1].flag = REJECTED;
					pathR2 = 0;
				}


				// 1. Both reads are path of the graph
				if(reads[j].flag != REJECTED && reads[j+1].flag != REJECTED){
					// 1.1. Both reads are part of a path and not a junction
					if(pathR1 && pathR2){
						// 1.1.1. Both reads are part of the same path
						if(pathR1 == pathR2){
							insert = ABS((ldistR1-ldistR2)) + reads[j+1].len;
							insert_array[insert_num] = insert;
							insert_num++;
							insert_sum += insert;
							files[i].avgInsert = insert_sum/insert_num;
							// 1.1.1.1. Correct insert size
							if(insert >= minIns && insert <= maxIns){
								if(debug) printf("1.1.1.1. ReadPair: %i, %i\n",j,j+1);
								// Not really interesting -> Should be the standard case!
								r1side = readside(G,j);
								r2side = readside(G,j+1);
								int dirCon = dirConcordance(pathR1,pathR2,r1side,r2side,pc_annoR1->lJunctionDist,pc_annoR2->lJunctionDist);

								if(dirCon != -1){
									// 1.1.1.1.1. correct directions
									if(files[i].oriPE == dirCon){
										if(debug) printf("1.1.1.1.1. ReadPair: %i, %i\n",j,j+1);
										if(verbose) printf("CON:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j,status_char[(int)G->read[j]->flag],compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist,r1side);
										if(verbose) printf("CON:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j+1,status_char[(int)G->read[j+1]->flag],compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist,r2side);
										if(verbose)	printf("Right Direction: %s, expected: %s\n",peOri[dirCon],peOri[files[i].oriPE]);
										reads[j].flag = CONCORDANT;
										reads[j+1].flag = CONCORDANT;
										co_num+=2;
									}
									// 1.1.1.1.2. wrong directions
									else{
										if(debug) printf("1.1.1.1.2. ReadPair: %i, %i\n",j,j+1);
										if(verbose) printf("DIS:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j,status_char[(int)G->read[j]->flag],compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist,r1side);
										if(verbose) printf("DIS:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j+1,status_char[(int)G->read[j+1]->flag],compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist,r2side);
										if(verbose) printf("Wrong Direction: %s, expected: %s\n",peOri[dirCon],peOri[files[i].oriPE]);
										reads[j].flag = DISCORDANT;
										reads[j+1].flag = DISCORDANT;
										di_num+=2;
									}

								}
								else{
									printf("Direction could not be detected\nAbort\n");
									exit(1);

								}
								if(verbose) printf("\t --> CONcordant Read Pair: insert size = %i\n",insert);
							}
							// 1.1.1.2. Wrong insert size
							else{
								if(debug) printf("1.1.1.2. ReadPair: %i, %i\n",j,j+1);
								reads[j].flag = DISCORDANT;
								reads[j+1].flag = DISCORDANT;
								di_num += 2;
								if(verbose) printf("DIS:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
								if(verbose) printf("DIS:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
								if(verbose) printf("\t --> DIScordant Read Pair: insert size = %i\n",insert);
							}
						}
						// 1.1.2. Both reads are part of different paths on same component
						// Real Scaffolding here (look for concordance)
						else if(paths[pathR1].component == paths[pathR2].component){
							if(debug) printf("1.1.2. ReadPair: %i, %i\n",j,j+1);
							// Tour to the corresponding path, collect insert size!!!
							interpathdist = 0;
							// walk from R1 read to the right if R2 was not found (concordantly), walk to the left
							if(verbose2) printf("SCAFFOLD reads:\n");
							if(verbose2) printf("UNK:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist,G->read[j]->dir);
							if(verbose2) printf("UNK:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist,G->read[j+1]->dir);
							if(verbose2) printf("TEST for Concordance\n");
							bread = G->read[j]->first;
							insert = connectPathsInit(G,files,reads,pathR1,pathR2,pc_annoR1,j+1,files[i].oriPE);
							// 1.1.2.1.1. Correct insert size and read directions
							if(insert > 0){
								if(debug) printf("1.1.2.1.1. ReadPair: %i, %i\n",j,j+1);
								r1side = readside(G,j);
								r2side = readside(G,j+1);

								if(verbose3) printf("CON:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j,status_char[(int)G->read[j]->flag],compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist,r1side);
								if(verbose3) printf("CON:\t -> ID: %i (flag: %c)\t comp: %i\t path: %i\t lDist: %i\t rDist: %i Dir: %i\n",j+1,status_char[(int)G->read[j+1]->flag],compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist,r2side);
								if(verbose2) printf("\t --> CONcordant Read Pair: insert size = %i\n",insert);
								if(verbose2) printf("Connect the path and count the observations:\n");
								connectPathsPairwise3();
								reads[j].flag = CONCORDANT;
								reads[j+1].flag = CONCORDANT;
								re_num += 2;
							}
							// 1.1.2.1.2. Wrong insert size and or wrong direction
							// -> Same treatment as 1.1.3 may usable for contig bridging
							else if(insert == 0){
								if(debug) printf("1.1.2.1.2. ReadPair: %i, %i\n",j,j+1);
								// TODO: write Function to catch bridging pairs and save bridge-information somehow
								if(verbose2) printf("DIS:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
								if(verbose2) printf("DIS:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
								char spurs = isBridge(G,pathR1,pathR2);
								if(spurs != -1){
									if(verboseBridgeYes) printf("YES 1 DISBRIDGE:\t-> ID: %i\t comp: %i\t path: %i\t lDist: %i (ID: %i)\t rDist: %i (ID: %i)\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,paths[pc_annoR1->pathID].leftJunction,pc_annoR1->rJunctionDist,paths[pc_annoR1->pathID].rightJunction);
									if(verboseBridgeYes) printf("YES 2 DISBRIDGE:\t-> ID: %i\t comp: %i\t path: %i\t lDist: %i (ID: %i)\t rDist: %i (ID: %i)\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,paths[pc_annoR2->pathID].leftJunction,pc_annoR2->rJunctionDist,paths[pc_annoR2->pathID].rightJunction);
									reads[j].flag = CONCORDANT;
									reads[j+1].flag = CONCORDANT;
									br_num += 2;
									buildBridge(G,files[i],reads,j,spurs,pathR1,pathR2,files[i].oriPE);
								}
								else{
									if(verboseBridge) printf("NO 1 DISBRIDGE:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
									if(verboseBridge) printf("NO 2 DISBRIDGE:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
									reads[j].flag = DISCORDANT;
									reads[j+1].flag = DISCORDANT;
									di_num += 2;
								}
							}
							// 1.1.2.1.3.
							else if(insert == -1){
								if(debug) printf("1.1.2.1.3. ReadPair: %i, %i\n",j,j+1);
								mu_num++;
							}
						}
						// 1.1.3. Both reads are part of different paths in different components.
						// ToDo: Look for orientation and distance to Component end -> Connect Components // Estimate bridge length by average insert size and distance to the contig end
						else{
							if(debug) printf("1.1.3. ReadPair: %i, %i\n",j,j+1);
							// TODO: write Function to catch bridging pairs and save bridge-information somehow
							char spurs = isBridge(G,pathR1,pathR2);
							if(spurs != -1){
								if(verboseBridgeYes) printf("YES 1 DIFFCOMP:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
								if(verboseBridgeYes) printf("YES 2 DIFFCOMP:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
								reads[j].flag = CONCORDANT;
								reads[j+1].flag = CONCORDANT;
								br_num += 2;
								buildBridge(G,files[i],reads,j,spurs,pathR1,pathR2,files[i].oriPE);
							}
							else{
								if(verboseBridge) printf("NO 1 DIFFCOMP:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
								if(verboseBridge) printf("NO 2 DIFFCOMP:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
								reads[j].flag = DISCORDANT;
								reads[j+1].flag = DISCORDANT;
								dicomp_num += 2;
							}
						}

					}
					// 1.2. At least one Read is a Junction
					else{
						if(debug) printf("1.2. ReadPair: %i, %i\n",j,j+1);
						un_num += 2;
						if(G->read[j]->flag == JUNCTION){
							if(verbose) printf("JUN:\t -> ID: %i\n",j);
						}
						else{
							if(verbose) printf("CON:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
						}
						if(G->read[j+1]->flag == JUNCTION){
							if(verbose) printf("JUN:\t -> ID: %i\n",j+1);
						}
						else{
							if(verbose) printf("CON:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
						}
					}
				}
				// 2. One read is not part of the graph
				// 2.1.
				else if(reads[j].flag == REJECTED && reads[j+1].flag == CONCORDANT){
					if(debug) printf("2.1. ReadPair: %i, %i\n",j,j+1);
					reads[j+1].flag = WIDOW;
					de_num ++;
					wi_num ++;
					if(verbose) printf("DEN:\t -> ID: %i\n",j);
					if(verbose) printf("WID:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j+1,compR2,pathR2,pc_annoR2->lJunctionDist,pc_annoR2->rJunctionDist);
				}
				// 2.2.
				else if(reads[j].flag == CONCORDANT && reads[j+1].flag == REJECTED){
					if(debug) printf("2.2. ReadPair: %i, %i\n",j,j+1);
					reads[j].flag = WIDOW;
					de_num ++;
					wi_num ++;
					if(verbose) printf("WID:\t -> ID: %i\t comp: %i\t path: %i\t lDist: %i\t rDist: %i\n",j,compR1,pathR1,pc_annoR1->lJunctionDist,pc_annoR1->rJunctionDist);
					if(verbose) printf("DEN:\t -> ID: %i\n",j+1);
				}
				// 3. Both reads are not part of the graph
				else{
					if(debug) printf("3. ReadPair: %i, %i\n",j,j+1);
					de_num += 2;
					if(verbose) printf("DEN:\t -> ID: %i\n",j);
					if(verbose) printf("DEN:\t -> ID: %i\n",j+1);
				}
			}


			insert_avg = insert_sum / insert_num;
			insert_var_sum = 0;
			for(j=0;j<insert_num;j++){
				insert_stddev_sum += ABS((insert_avg - insert_array[j]));
				insert_var_sum +=  ABS((insert_avg - insert_array[j])) * ABS((insert_avg - insert_array[j]));
			}
			insert_var = (float)insert_var_sum/insert_num;
			insert_stddev = insert_stddev_sum/insert_num;
			printf("Library Stats:\n");
			printf("Concordant Reads pairs: %i\n",insert_num);
			printf("Avg insert size:        %i\n",insert_avg);
			printf("Stddev of insert size:  %i\n",insert_stddev);
			printf("Var of insert size:     %.2f\n",insert_var);
			printf("\n");
			printf("CONCORDANT READS : %i\n",co_num);
			printf("BRIDGING   READS : %i\n",br_num);
			printf("REPRESOLVE READS : %i\n",re_num);
			printf("DISCORDANT READS : %i\n",di_num);
			printf("DIFFCOMP   READS : %i\n",dicomp_num);
			printf("MULTIPLE P READS : %i\n",mu_num);
			printf("REJECTED   READS : %i\n",de_num);
			printf("WIDOWED    READS : %i\n",wi_num);
			printf("UNKNOWN    READS : %i\n",un_num);
			printf("FOUND      READS : %i\n",foundsum);
			printf("NOT FOUND  READS : %i\n",nfoundsum);
			co_num_tot += co_num;
			co_num_tot += br_num;
			co_num_tot += re_num;
			di_num_tot += di_num;
			de_num_tot += de_num;
			wi_num_tot += wi_num;
			un_num_tot += un_num;
			dicomp_num_tot += dicomp_num;
			free(insert_array);
		}
		else{
			printf("\n##### SE LIB #####\n");
			printf("From ID: %i to %i\n",startID,endID);
			int si_num = 0;
			int de_num = 0 ;
			for(i=startID;i<=endID;i++){
				if(G->read[i]){
					reads[i].flag = SINGLE;
					si_num ++;
				}
				else{
					reads[i].flag = REJECTED;
					de_num++;
				}
			}
			printf("REJECTED   READS : %i\n",de_num);
			printf("SINGLE     READS : %i\n",si_num);
			de_num_tot += de_num;
			si_num_tot += si_num;
		}
		pe = 0;
	}

	printf("\n");
	printf("Collective Statistica:\n");
	printf("CONCORDANT READS : %i\n",co_num_tot);
	printf("DISCORDANT READS : %i\n",di_num_tot);
	printf("DIFFCOMP   READS : %i\n",dicomp_num_tot);
	printf("REJECTED   READS : %i\n",de_num_tot);
	printf("WIDOWED    READS : %i\n",wi_num_tot);
	printf("UNKNOWN    READS : %i\n",un_num_tot);
	printf("SINGLE     READS : %i\n",si_num_tot);
	printf("\n");

	annotatePaths();
	connectPathStats();

}

