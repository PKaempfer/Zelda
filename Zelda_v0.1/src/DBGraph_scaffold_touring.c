#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>
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


inline static void hierholzer_scaffConnection(struct pathEdge* pathedge, int ID1, int ID2, char* found1, char* found2){
	struct pathEdge* sibedge;
	while(pathedge){
		if(pathedge->ID == ID1){
			(*found1)=1;
		}
		if(pathedge->ID == ID2){
			(*found2)=1;
		}
		if((*found1) && (*found2)) return;
		sibedge = pathedge->sibl;
		while(sibedge){
			hierholzer_scaffConnection(sibedge,ID1,ID2,found1,found2);
			sibedge = sibedge->sibl;
		}
		pathedge = pathedge->next;
	}
}

inline static void hierholzer_altPath(struct jPath** jPath, int* parents, int new, char* flags){
	char verbose = 0;
	char verbose2 = 0;
	int i,j,k,l,m;
	char found1 = 0;
	char found2 = 0;
	struct pathEdge* pathedge;
	for(i=1;i<new;i++){
		l=i;
		while(l+1<new && parents[l+1]==parents[i]){
			l++;
		}
		if(l>i){
			for(j=i;j<l;j++){
				for(k=j+1;k<=l;k++){
					if(jPath[j]->eJID == jPath[k]->eJID){
						if(verbose) printf("Found Altarnative Parts (Paths: %i,%i)\n",jPath[j]->pathID,jPath[k]->pathID);
						if(verbose) printf("Check with Circle Path (l: %i i: %i)\n",l,i);
						found1 = 0;
						found2 = 0;
						for(m=1;m<new;m++){
							found1 = 0;
							found2 = 0;
							pathedge = paths[jPath[m]->pathID].leftPath;
							hierholzer_scaffConnection(pathedge,jPath[j]->pathID,jPath[k]->pathID,&found1,&found2);
							if(found1 && found2) continue;
							pathedge = paths[jPath[m]->pathID].rightPath;
							hierholzer_scaffConnection(pathedge,jPath[j]->pathID,jPath[k]->pathID,&found1,&found2);
							// ToDo: Maybe to dangerous / don't use unflagged nodes
							if(!flags[m]){
								found1 = !found1;
								found2 = !found2;
							}
							if(found1 && !found2){
								if(verbose2) printf("Path %i, supports %i but not %i\n",jPath[m]->pathID,jPath[j]->pathID,jPath[k]->pathID);
								flags[k]=0;
							}
							if(!found1 && found2){
								if(verbose2) printf("Path %i, supports %i but not %i\n",jPath[m]->pathID,jPath[k]->pathID,jPath[j]->pathID);
							}
						}
					}
				}
			}
		}
	}
	for(i=1;i<new;i++){
		if(!flags[parents[i]]){
			flags[i]=0;
		}
	}
}

inline static char hierholzer_uniqCircles(int* parents, int new, char* flags){
	int i,j,k,l;
	for(i=1;i<new;i++){
		l=i;
		while(l+1<new && parents[l+1]==parents[i]){
			l++;
		}
		if(l>i){
			for(j=i;j<l;j++){
				for(k=j+1;k<=l;k++){
					if(flags[j] && flags[k]) return 0;
				}
			}
		}
	}
	return 1;
}

inline static void hierholzer_cleanCircle(int* parents, int* goals, int goals_i, char* flags){
	int i;
	int parentPos;
	for(i=0;i<goals_i;i++){
		parentPos = goals[i];
		while(parentPos){
			if(flags[parentPos]) break;
			flags[parentPos] = 1;
			parentPos = parents[parentPos];
		}
	}
}

inline static void hierholzer_backtracking(struct jPath** jPath, int* parents, int* goals, int goals_i, char* flags, struct reads* reads, int* circleID){
	int i;
	int parent;
	struct j_anno* j_anno;
	struct contigCircle* circle;
	for(i=0;i<goals_i;i++){
		if(flags[goals[i]]){
			parent = goals[i];
			while(parent){
				// EndJuncion
				j_anno = (struct j_anno*)reads[jPath[parent]->eJID].annotation;
				circle = (struct contigCircle*)malloc(sizeof(struct contigCircle));
				circle->ID = (*circleID);
				circle->dest = jPath[parent]->sJID;
				circle->pathID = jPath[parent]->pathID;
				circle->next = j_anno->circle;
				j_anno->circle = circle;

				// StartJunction
				j_anno = (struct j_anno*)reads[jPath[parent]->sJID].annotation;
				circle = (struct contigCircle*)malloc(sizeof(struct contigCircle));
				circle->ID = (*circleID);
				circle->dest = jPath[parent]->eJID;
				circle->pathID = jPath[parent]->pathID;
				circle->next = j_anno->circle;
				j_anno->circle = circle;

				// Path
				paths[jPath[parent]->pathID].circfreq++;
				paths[jPath[parent]->pathID].freq--;

				parent = parents[parent];
			}
			(*circleID)++;
		}
	}
}

inline static void hierholzerTourPath(int i, char inDir , struct reads* reads){
	char verbose = 0;
	if(verbose) printf("Checkpoint Hierholzer Tour\n");
	int j;
	struct j_anno* j_anno;

	int nodesize = 1000;
	int goalsize = 10;
	int goal_i = 0;
	int goalNum;
	int start, end, new = 0;
	int tempparent;
	char side;
	char set;
	static int circleID = 0;

	struct jPath* edge;
	struct jPath* newedge;
	struct jPath** jPath = (struct jPath**)malloc(sizeof(struct jPath*)*nodesize);
	int* parents = (int*)malloc(sizeof(int)*nodesize);
	char* flags = (char*)malloc(nodesize);
	int* goals = (int*)malloc(sizeof(int)*goalsize);

	goal_i = 0;
	j_anno = (struct j_anno*)reads[i].annotation;
	if(j_anno->inDegree == j_anno->outDegree && j_anno->inDegree > 1){
		printf("Node: %i is a potential startpoint of a contig circle.\n",i);
		new = 1;
		start = 1;

		if(inDir) edge = j_anno->outEdge;
		else edge = j_anno->inEdge;
		while(edge && new < nodesize && goal_i < goalsize){
			if(paths[edge->pathID].freq){
				jPath[new] = edge;
				parents[new] = 0;
				flags[new] = 0;
				if(edge->eJID == i) goals[goal_i] = new;
				new++;
			}
			edge = edge->next;
		}

		// Start touring: BFS
		end = new-1;
		while(start<=end){
			while(start<=end){
				if(jPath[start]->eJID == i){
					start++;
					continue;
				}
				j_anno = (struct j_anno*)reads[jPath[start]->eJID].annotation;
				side = 0;
				edge = jPath[start];
				newedge = j_anno->inEdge;
				while(newedge){
					if(newedge->pathID == edge->pathID){
						side = 1;
						break;
					}
					newedge = newedge->next;
				}
				if(side) newedge = j_anno->outEdge;
				else newedge = j_anno->inEdge;
				while(newedge && new < nodesize && goal_i < goalsize){
					set = 1;
					tempparent = start;
					while(tempparent){
						if(jPath[tempparent]->pathID == newedge->pathID){
							if(verbose){
								printf("Found Inner circle at %i\n",newedge->pathID);
							}
							set = 0;
							break;
						}
						tempparent = parents[tempparent];
					}
					if(set && paths[newedge->pathID].freq){
						if(newedge->eJID != newedge->sJID){
							jPath[new] = newedge;
							parents[new] = start;
							flags[new] = 0;
							if(i==newedge->eJID) goals[goal_i++] = new;
							new++;
						}
					}
					newedge = newedge->next;
				}
				start++;
			}
			end = new -1;
		}
		if(verbose){
			printf("Goals:");
			for(j=0;j<goal_i;j++){
				printf("\t%i",goals[j]);
			}
			printf("\n");
			printf("Circle Tree:\n");
			for(j=1;j<new;j++){
				printf("\t%i (%i)",jPath[j]->pathID,parents[j]);
			}
			printf("\n");
		}
		// 2. Tour: identify alternative parts in the circle and validate their possibility
		hierholzer_cleanCircle(parents,goals,goal_i,flags);
		hierholzer_altPath(jPath,parents,new,flags);
		if(verbose){
			printf("Circle Tree:\n");
			for(j=1;j<new;j++){
				if(flags[j]) printf("\t%i (%i)",jPath[j]->pathID,parents[j]);
				else printf("\t%i (%i)",jPath[j]->pathID*-1,parents[j]);
			}
			printf("\n");
		}
		// Check if circle is unique
		goalNum = 0;
		for(j=0;j<goal_i;j++){
			if(flags[goals[j]]) goalNum++;
		}
		// ToDo: Look if any circle is fully part of another one
		if(hierholzer_uniqCircles(parents,new,flags)){
			hierholzer_backtracking(jPath,parents,goals,goal_i,flags,reads,&circleID);
			printf("Found %i unique circle paths (Circles: %i)\n",goalNum,circleID);
			// Backtracking - Find correct circle
		}
	}
}

/**
 * Searches for circles in the String graph and tag them. This circles are not investigated in the
 * Scaffolding tour but substituted when one of them nodes are firstly visited.
 */
void hierholzerTourAll(struct myovlList* G, struct reads* reads){
	printf("Checkpoint Hierholzer Tour\n");
	int i;
//	struct j_anno* j_anno;

//	int nodesize = 10000;
//	int goalsize = 100;
//	int goal_i = 0;
//	int goalNum;
//	int start, end, new = 0;
//	int tempparent;
//	char side;
//	char set;
//	int circleID = 0;

//	struct jPath* edge;
//	struct jPath* newedge;
//	struct jPath** jPath = (struct jPath**)malloc(sizeof(struct jPath*)*nodesize);
//	int* parents = (int*)malloc(sizeof(int)*nodesize);
//	char* flags = (char*)malloc(nodesize);
//	int* goals = (int*)malloc(sizeof(int)*goalsize);

//	hierholzerTourPath(178515,0,reads);
//	hierholzerTourPath(178515,1,reads);

//	exit(1);

	for(i=1;i<=G->V;i++){
		if(G->read[i] && G->read[i]->flag == JUNCTION){
			hierholzerTourPath(i,0,reads);
			hierholzerTourPath(i,1,reads);
//			exit(1);
		}
	}
}

static inline void prepPathsFlag(){
	int i;
	for(i=1;i<pathsNum;i++){
		paths[i].flag = paths[i].freq  + paths[i].circfreq;
	}
}

static inline int scaffold_lookForward(char right, int current, struct contigScaff* uniPath, int pos, int len){
	char verbose = 1;
	char verbose2 = 1;
	if(verbose) printf("Path look forward\n");
//	int dist = len - pos;
	struct pathEdge* edge;
	struct pathEdge* sibl;
	char found = 0;
	if(right) edge = paths[current].rightPath;
	else edge = paths[current].leftPath;
	pos++;
	while(edge){
		if(pos < len){
			if(verbose) printf("Pos %i < Len %i\n",pos,len);
			if(edge->sibl){
				if(verbose) printf("Edge has Sibls\n");
				sibl = edge;
				found = 0;
				while(sibl){
					if(sibl->ID == uniPath[pos].ID){
						if(verbose) printf("Correct Sibl found: %i\n",sibl->ID);
						found = 1;
						pos++;
						edge = sibl->next;
						break;
					}
					sibl = sibl->sibl;
				}
				if(found) continue;
			}
			else{
				if(edge->ID == uniPath[pos].ID){
					// Go ahead -> Continue
					if(verbose) printf("Same ID: %i\n",edge->ID);
				}
				else{
					if(verbose) printf("Wrong ID: %i != %i\n",edge->ID,uniPath[pos].ID);
					break;
				}
			}
		}
		else{
			if(verbose) printf("NewPos");
			if(edge->sibl){
				if(verbose) printf("Siblings, break\n");
				break;
			}
			else{
				// ToDo: Implement double check: lock backward: if not possible to reach the path: set node on an index, to the next
				// folk, then go not this way and delete the node from the index
				if(verbose) printf("Unique, elongate path: %i\n",edge->ID);
				uniPath[len].ID = edge->ID;
				len++;
			}
		}
		pos++;
		edge = edge->next;
	}
	if(verbose2){
		printf("List of the current path: \n");
		for(int i=0;i<len;i++){
			printf("\t%i",uniPath[i].ID);
		}
		printf("\n");
	}

	return len;
}

//static inline void scaffold_gocirclePath(struct scaffold_set* aS, struct contigCircle* circle, int circleID, char right, int current, struct contigScaff* uniPath, int pos, int uniPathLen){
//	printf("CHECKPOINT: Start running the circle with ID: %i over path: %i\n",circleID,circle->pathID);
//}

static inline void scaffold_deleteJunctionEdge(int pathID, struct reads* reads, int circleID){
	printf("CHECKPOINT: Delete junction edges to paths\n");
	// delete leftside edge
	char del = 0;
	char deledge = 0;
	if(paths[pathID].circfreq + paths[pathID].freq == 1) deledge = 1;
	struct contigCircle* circle;
	struct contigCircle* pcircle;
	int jun = paths[pathID].leftJunction;
	struct j_anno* j_anno = (struct j_anno*)reads[jun].annotation;
	struct jPath* jPath = j_anno->inEdge;
	struct jPath* pjPath = NULL;
	while(jPath){
		if(jPath->pathID == pathID){
			if(deledge){
				if(pjPath) pjPath->next = jPath->next;
				else j_anno->inEdge = jPath->next;
				free(jPath);
			}
			j_anno->inDegree--;
			del = 1;
			break;
		}
		pjPath = jPath;
		jPath = jPath->next;
	}
	if(!del){
		jPath = j_anno->outEdge;
		pjPath = NULL;
		while(jPath){
			if(jPath->pathID == pathID){
				if(deledge){
					if(pjPath) pjPath->next = jPath->next;
					else j_anno->outEdge = jPath->next;
					free(jPath);
				}
				j_anno->outDegree--;
				del = 1;
				break;
			}
			pjPath = jPath;
			jPath = jPath->next;
		}
	}
	if(!del){
		printf("Path of the junctions (%i) not found (InDegree: %i / Outdegree: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);
		exit(1);
	}

	if(circleID>=0){
		circle = j_anno->circle;
		pcircle = NULL;
		while(circle){
			if(circle->ID == circleID && circle->pathID == pathID){
				if(pcircle) pcircle->next  = circle->next;
				else j_anno->circle = circle->next;
				free(circle);
				break;
			}
			pcircle = circle;
			circle = circle->next;
		}
	}
	printf("Junction: %i (in: %i, out: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);

	// rigth side
	jun = paths[pathID].rightJunction;
	j_anno = (struct j_anno*)reads[jun].annotation;
	jPath = j_anno->outEdge;
	pjPath = NULL;
	del = 0;
	while(jPath){
		if(jPath->pathID == pathID){
			if(deledge){
				if(pjPath) pjPath->next = jPath->next;
				else j_anno->outEdge = jPath->next;
				free(jPath);
			}
			j_anno->outDegree--;
			del = 1;
			break;
		}
		pjPath = jPath;
		jPath = jPath->next;
	}
	if(!del){
		jPath = j_anno->inEdge;
		pjPath = NULL;
		while(jPath){
			if(jPath->pathID == pathID){
				if(deledge){
					if(pjPath) pjPath->next = jPath->next;
					else j_anno->inEdge = jPath->next;
					free(jPath);
				}
				j_anno->inDegree--;
				del = 1;
				break;
			}
			pjPath = jPath;
			jPath = jPath->next;
		}
	}
	if(circleID>=0){
		circle = j_anno->circle;
		pcircle = NULL;
		while(circle){
			if(circle->ID == circleID && circle->pathID == pathID){
				if(pcircle) pcircle->next  = circle->next;
				else j_anno->circle = circle->next;
				free(circle);
				break;
			}
			pcircle = circle;
			circle = circle->next;
		}
	}
	if(circleID>=0) paths[pathID].circfreq--;
	else paths[pathID].freq--;
	printf("Junction: %i (in: %i, out: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);
}

int* circleList = NULL;
int circleListLen = 0;
int circleListMaxLen;

int* scaffold_gocirclePath(struct scaffEdge* scaffedge, struct scaffold_set* aS, struct jPath* jPath, int circleID, char right, int jpos, struct contigScaff* uniPath, int* pos_len, int uniPathMaxLen, struct reads* reads){
	printf("CHECKPOINT: Start running the circle with ID: %i over path: %i\n",circleID,jPath->pathID);
	// Follow the circle

	int njpos = jpos;
	int lastpath;
	struct jPath* tempJpath = jPath;
	struct scaffEdge* scaffedgenew;
	struct j_anno* j_anno;
	int circleID1;
	int circleID2;
	struct contigCircle* correctCircle = NULL;
	char out;

	do{
		if(paths[jPath->pathID].leftJunction == njpos){
			njpos = paths[jPath->pathID].rightJunction;
			right = 1;
		}
		else{
			njpos = paths[jPath->pathID].leftJunction;
			right = 0;
		}
		// set path to aS
		pos_len[0]++;
		printf("-----> Set new edge: %i to Junction %i (pos: %i uniLength: %i)\n",jPath->pathID,njpos,pos_len[0],pos_len[1]);
		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
		scaffedgenew->ID = jPath->pathID;
		scaffedgenew->depth = pos_len[0];
		scaffedgenew->len = paths[jPath->pathID].len;
		scaffedgenew->targetJunction = njpos;
		scaffedgenew->bridge = NULL;
		scaffedgenew->next = NULL;
		scaffedge->next = scaffedgenew;
		scaffedge = scaffedgenew;
		lastpath = jPath->pathID;
		aS->scaff[aS->num].len += paths[jPath->pathID].len;
		uniPath[pos_len[0]].ID = jPath->pathID;
		if(pos_len[1]<=pos_len[0]) pos_len[1] = pos_len[0]+1;

		pos_len[1] = scaffold_lookForward(right,jPath->pathID,uniPath,pos_len[0],pos_len[1]);

		// search new jPath and next junction
		j_anno = (struct j_anno*)reads[njpos].annotation;
		jPath = j_anno->inEdge;
		tempJpath = jPath;
		while(jPath){
			if(jPath->pathID == lastpath){
				jPath = j_anno->outEdge;
				out = 1;
				break;
			}
			jPath = jPath->next;
		}
		if(!jPath){
			jPath = j_anno->inEdge;
			out = 0;
		}

		scaffold_deleteJunctionEdge(scaffedgenew->ID,reads,circleID);

		if(jpos == njpos) return pos_len;

		char oncurse = 0;
		do{
			int altCirleID = -1;
			int circleNum = 0;
			printf("Pos: %i, Len: %i\n",pos_len[0]+1,pos_len[1]);
			if(pos_len[0]+1<pos_len[1]){
				correctCircle = j_anno->circle;
				while(correctCircle){
					if(correctCircle->pathID == uniPath[pos_len[0]+1].ID){
						altCirleID = correctCircle->ID;
					}
					if(correctCircle->ID == circleID && correctCircle->pathID == uniPath[pos_len[0]+1].ID) oncurse = 1;
					printf("Circle: %i over Path %i -> Next path: %i \n",correctCircle->ID,correctCircle->pathID,uniPath[pos_len[0]+1].ID);
					circleNum++;
					correctCircle = correctCircle->next;
				}
				if(!oncurse && altCirleID!=-1){
					printf("Follow another Circle\n");
					tempJpath = jPath;
					while(tempJpath){
						if(tempJpath->pathID == uniPath[pos_len[0]+1].ID){
							break;
						}
						tempJpath = tempJpath->next;
					}
					if(tempJpath){
						printf(" --> Go Intermediate Circle (CID: %i)\n",altCirleID);
						circleList[circleListLen++] = circleID;
						if(circleListLen == circleListMaxLen){
							circleListMaxLen *= 2;
							circleList = (int*)realloc(circleList,sizeof(int)*circleListMaxLen);
						}
						pos_len = scaffold_gocirclePath(scaffedge,aS,tempJpath,altCirleID,right,njpos,uniPath,pos_len,uniPathMaxLen,reads);
						while(scaffedge->next){
							scaffedge = scaffedge->next;
						}
						circleListLen--;
						printf(" --> Intermediate Circle Finished\n");
					}
					else{
						printf("Circle in Circle, But correct jpath not found!\n");
						exit(1);
					}
				}
				else oncurse = 1;
			}
			else oncurse=1;
			printf("Number of Circles in this Junction (J: %i): %i\n",njpos,circleNum);
		} while(!oncurse);

		correctCircle = j_anno->circle;
		circleID1 = 0;
		circleID2 = 0;
		while(correctCircle){
			printf("Circle: %i over Path %i\n",correctCircle->ID,correctCircle->pathID);
			if(correctCircle->ID == circleID){
				if(circleID1){
					circleID2 = correctCircle->pathID;
					printf("Found 2nd\n");
					break;
				}
				else{
					circleID1 = correctCircle->pathID;
					printf("Found 1st\n");
				}
			}
			correctCircle = correctCircle->next;
		}

//		tempJpath = jPath;
//
//		while(tempJpath){
//			tempJpath = tempJpath->next;
//		}

//		jPath = tempJpath;

		if(out) jPath = j_anno->outEdge;
		else jPath = j_anno->inEdge;
		while(jPath){
			printf("Comp: %i - %i / %i\n",jPath->pathID,circleID1,circleID2);
			if(jPath->pathID == circleID1 || jPath->pathID == circleID2) break;
			jPath = jPath->next;
		}

		// jPath seems to be NULL at some point; Find out why!

	} while(njpos != jpos);
	return pos_len;
}

static inline char scaffold_comparePaths(struct contigCircle* circle1, struct contigCircle* circle2, int junction, struct reads* reads){
	char verbose = 1;

	int end = junction;

	int ID1 = circle1->ID;
	int ID2 = circle2->ID;
	int target = junction;

	struct j_anno* j_anno;
	char found1;
	char found2;

	int pathID = circle1->pathID;

	do{

		j_anno = (struct j_anno*)reads[target].annotation;

		// Next Path of circle 1
		found1 = 0;
		circle1 = j_anno->circle;
		while(circle1){
			if(circle1->ID == ID1){
				if(circle1->pathID != pathID || found1 == -1){
					break;
				}
				else{
					found1 = -1;
				}
			}
			circle1 = circle1->next;
		}

		// Next Path of cirlce 2
		found2 = 0;
		circle2 = j_anno->circle;
		while(circle2){
			if(circle2->ID == ID2){
				if(circle2->pathID != pathID || found2 == -1){
					break;
				}
				else{
					found2 = -1;
				}
			}
			circle2 = circle2->next;
		}

		// Compare
		if(!circle1 || !circle2){
			printf("Circle Not Found\nAbort\n");
			exit(1);
		}

		// If same go further
		if(circle1->pathID == circle2->pathID){
			if(verbose) printf("Same next Path: %i (junction: %i)\n",circle1->pathID,circle1->dest);
			target = circle1->dest;
			pathID = circle1->pathID;
		}
		// If different return 0 directly
		else{
			if(verbose) printf("Different next Paths: %i / %i (junction: %i / %i)\n",circle1->pathID,circle2->pathID,circle1->dest,circle2->dest);
			return 0;
		}

	} while(end!=target);

	return 1;
}

static inline int* scaffold_findcirclePath(struct scaffEdge* scaffedge, struct jPath* jPath, int jpos, char initRight, struct j_anno* j_anno,int nextpathNum, int uniPathMaxLen, struct contigScaff* uniPath, struct scaffold_set* aS, struct reads* reads,int* pos_len){

	struct path* path;
	struct j_anno* j_annoTemp;
	int circNum;
	int foundID = 0;
	int circleID1;
	int circleID2;
	char goon;
	struct contigCircle* circle;
	struct contigCircle* correctCircle = NULL;


	if(nextpathNum == 1){
		// Go definitely the entire circle
		correctCircle = j_anno->circle;
		foundID = correctCircle->ID;
//		scaffold_gocirclePath(aS,correctCircle,foundID,initRight,jpos,uniPath,pos,uniPathLen);
		pos_len = scaffold_gocirclePath(scaffedge, aS,jPath,foundID,initRight,jpos,uniPath,pos_len,uniPathMaxLen,reads);
	}
	else{
		for(int i=0; i < pos_len[1]; i++){
			path = &paths[uniPath[i].ID];
			if(path->circfreq){
				j_annoTemp = (struct j_anno*)reads[path->leftJunction].annotation;
				circle = j_annoTemp->circle;
				circNum = 0;
				while(circle){
					if(circle->pathID == uniPath[i].ID){
						foundID = circle->ID;
						circNum++;
					}
					circle = circle->next;
				}
				if(circNum>1){
					goon = -1;
					printf("Test if all possible circles runing the same paths?\n");
					circle = j_annoTemp->circle;
					while(circle){
						if(circle->pathID == uniPath[i].ID){
							correctCircle = circle->next;
							while(correctCircle){
								if(correctCircle->ID != circle->ID && correctCircle->pathID == uniPath[i].ID){
									printf("Compare: Path %i / %i (circle: %i/ %i)\n",circle->pathID,correctCircle->pathID,circle->ID,correctCircle->ID);
									if(goon == -1){
										goon = scaffold_comparePaths(circle,correctCircle,circle->dest,reads);
									}
									else{
										if(!scaffold_comparePaths(circle,correctCircle,circle->dest,reads)) goon = 0;
									}
									if(!goon) break;
								}
								correctCircle = correctCircle->next;
							}
							if(!goon) break;
						}
						circle = circle->next;
					}
					if(goon) circNum = 1;
				}
				if(circNum == 1){
					printf("Found The correct from path %i circle with ID: %i\n",uniPath[i].ID,foundID);
					circleID1 = 0;
					circleID2 = 0;
					correctCircle = j_anno->circle;
					while(correctCircle){
						printf("Circle: %i over Path %i\n",correctCircle->ID,correctCircle->pathID);
						if(correctCircle->ID == foundID){
							if(circleID1){
								circleID2 = correctCircle->pathID;
								printf("Found 2nd\n");
								break;
							}
							else{
								circleID1 = correctCircle->pathID;
								printf("Found 1st\n");
							}
						}
						correctCircle = correctCircle->next;
					}
					while(jPath){
						if(jPath->pathID == circleID1 || jPath->pathID == circleID2) break;
						jPath = jPath->next;
					}
					if(jPath) pos_len = scaffold_gocirclePath(scaffedge,aS,jPath,foundID,initRight,jpos,uniPath,pos_len,uniPathMaxLen,reads);
					else{
						printf("No Correct JPATH found: Abort\n");
						exit(1);
					}
					break;
				}
			}
		}
	}
	return pos_len;
}

static inline void scaffold_uniqPath(char initRight, char casus, int i, int nextpathNum, int jpos, struct jPath* jPath, struct reads* reads , struct scaffold_set* aS){
	struct scaffEdge* scaffedge = aS->scaff[aS->num].first;
	struct scaffEdge* scaffedgenew;
    struct jPath* tempjPath;
    int difpathNum;
    int lastPath;
    struct j_anno* j_anno;
    struct contigCircle* circle;
    int circleNum;
    int depth = 0;
    int old_depth;
//    char casus; // 0 if outgoing, 1 if incomming
    int* pos_len = (int*)malloc(sizeof(int)*2);
    int lastjpos;

	int uniPathLen = 1;
    static int uniPathMaxLen = 1000;
    static struct contigScaff* uniPath = NULL;
    if(!uniPath){
    	uniPath = (struct contigScaff*)malloc(sizeof(struct contigScaff)*uniPathMaxLen);
    }

    uniPath[0].ID = i;

	uniPathLen = scaffold_lookForward(initRight,i,uniPath,0,uniPathLen);
	if(uniPathLen > uniPathMaxLen - 20){
		uniPathMaxLen *= 2;
		struct contigScaff* temp = (struct contigScaff*)realloc(uniPath,sizeof(struct contigScaff)*uniPathMaxLen);
		if(temp) uniPath = temp;
		else{
			printf("Reallocation of UniPath failed -> ABORT\n");
			exit(1);
		}
	}

	while(1){
//		tempjPath = jPath;
		// Use while loop

    	if(nextpathNum) printf("Number of following paths: %i (TargetJunction: %i)\n",nextpathNum,jpos);
    	else{
    		printf("end of Path reached\n: %i (TargetJunction: %i)\n",nextpathNum,jpos);
    		break;
    	}
    	tempjPath = jPath;
    	difpathNum = 0;
    	while(jPath){
    		printf("\t-> Path %i, (Circle %i)\n",jPath->pathID,paths[jPath->pathID].circfreq);
    		jPath = jPath->next;
    		difpathNum++;
    	}


		j_anno = (struct j_anno*)reads[jpos].annotation;

		jPath = tempjPath;
		if(difpathNum>1 && uniPathLen > depth+1){
			printf("DifPaths: %i\n",difpathNum);
			while(jPath){
				if(jPath->pathID == uniPath[depth+1].ID){
					printf("Correct Path Found: %i\n",jPath->pathID);
					difpathNum = 1;
					break;
				}
				jPath = jPath->next;
			}
		}


    	// Next Path is unique
    	if(difpathNum == 1){
    		while(paths[jPath->pathID].circfreq){
//        		jPath = tempjPath;
    			printf("Start Circle searching for circles (CicleNum: %i, restNum: %i)!\n",paths[jPath->pathID].circfreq,paths[jPath->pathID].freq);
    	    	// Look if there is any circle. Is It possible to assign uniquely one of the circles?
    	    	if(paths[jPath->pathID].circfreq){
    	    		printf("Start Circle searching for circles!\n");
        			j_anno = (struct j_anno*)reads[jpos].annotation;
        			old_depth = depth;
        			pos_len[0]=depth;
        			pos_len[1]=uniPathLen;
        			if(!circleList){
        				circleListMaxLen = 50;
        				circleList = (int*)malloc(sizeof(int)*circleListMaxLen);
        			}
    	    		pos_len = scaffold_findcirclePath(scaffedge,jPath,jpos,initRight,j_anno,nextpathNum,uniPathMaxLen,uniPath,aS,reads,pos_len);
    	    		circleListLen = 0;
    	    		depth = pos_len[0];
    	    		uniPathLen = pos_len[1];
    	    		if(depth != old_depth){
    		    		while(scaffedge->next){
    		    			scaffedge = scaffedge->next;
    		    		}
    		    		printf("Circles paths successfully end\n\n\n");

    		    		printf("New Cirlc???\n");
    	    		}
    	    		else{
    	    			printf("Circle not worked\n");
    	    			break;
    	    		}
    	    		if(!casus && !j_anno->outDegree) break;
    	    		if(casus && !j_anno->inDegree) break;
        			printf("Start Circle searching on path %i for circles (CicleNum: %i)!\n",jPath->pathID,paths[jPath->pathID].circfreq);

    	    	}
    		}

    		if(!casus && !j_anno->outDegree){
    			printf("Junction: %i -> CASUS: %i Out: %i\n",jpos,casus,j_anno->outDegree);
    			break;
    		}
    		if(casus && !j_anno->inDegree){
    			printf("Junction: %i -> CASUS: %i Out: %i\n",jpos,casus,j_anno->inDegree);
    			break;
    		}

//    		jPath = tempjPath;
    		lastPath = jPath->pathID;
    		lastjpos = jpos;

    		if(paths[jPath->pathID].leftJunction == jpos){
    			jpos = paths[jPath->pathID].rightJunction;
    			initRight = 1;
    		}
    		else{
    			jpos = paths[jPath->pathID].leftJunction;
    			initRight = 0;
    		}
    		j_anno = (struct j_anno*)reads[jpos].annotation;

    		depth++;
    		printf("-----> Set new edge: %i to Junction %i (unipos: %i uniLength: %i)\n",jPath->pathID,jpos,depth,uniPathLen);
    		scaffedgenew = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
    		scaffedgenew->ID = jPath->pathID;
    		scaffedgenew->len = paths[jPath->pathID].len;
    		scaffedgenew->bridge = NULL;
    		scaffedgenew->next = NULL;
    		scaffedgenew->depth = depth;
    		scaffedgenew->targetJunction = jpos;
    		scaffedge->next = scaffedgenew;
    		scaffedge = scaffedgenew;
    		aS->scaff[aS->num].len += paths[jPath->pathID].len;
    		uniPath[depth].ID = jPath->pathID;
    		if(uniPathLen<=depth) uniPathLen = depth+1;
    		uniPathLen = scaffold_lookForward(initRight,jPath->pathID,uniPath,depth,uniPathLen);



    		// Find the right edges of the next junction
    		jPath = NULL;
    		nextpathNum = 0;

			if(jpos == lastjpos){
				char same = 0;
				printf("Last Path was a selfloop\n");
				jPath = j_anno->inEdge;
				while(jPath){
					if(jPath->pathID == lastPath){
						same = 1;
						break;
					}
					jPath = jPath->next;
				}
				jPath = j_anno->outEdge;
				while(jPath){
					if(jPath->pathID == lastPath){
						if(same) same = 2;
						else same = 1;
						break;
					}
					jPath = jPath->next;
				}
				if(same == 2){
					printf("Was straight\n");
				}
				else{
					printf("Was Reverse (Go way around)\n");
					casus = !casus;
				}
			}
			else{
				jPath = j_anno->inEdge;
				casus = 0;
				while(jPath){
					if(jPath->pathID == lastPath){
						printf("--> Next path is OUTgoing\n");
						jPath = j_anno->outEdge;
						if(!jPath){
							printf("Next Path is outgoing: But is Empty! Path End\n");
							casus = -1;
							break;
						}
						casus = 0;
						break;
					}
					jPath = jPath->next;
				}
				if(casus == -1){
					scaffold_deleteJunctionEdge(scaffedgenew->ID,reads,-1);
					break;
				}
				if(!jPath){
					printf("--> Next path is INcomming\n");
					jPath = j_anno->inEdge;
					if(!jPath){
						scaffold_deleteJunctionEdge(scaffedgenew->ID,reads,-1);
						break;
					}
					casus = 1;
				}
			}
			scaffold_deleteJunctionEdge(scaffedgenew->ID,reads,-1);
			if(!casus){
				printf("--> Next path is OUTgoing\n");
				nextpathNum = j_anno->outDegree;
				jPath = j_anno->outEdge;
				printf("NexpathNum: %i\n",nextpathNum);

			}
			else{
				printf("--> Next path is INcomming\n");
				nextpathNum = j_anno->inDegree;
				jPath = j_anno->inEdge;
				printf("NexpathNum: %i\n",nextpathNum);
			}
    	}
    	else{
    		// If I reached a circle containing Junction test whether to go one of the circles first
    		if(j_anno->circle){
    			circleNum = 0;
    			circle = j_anno->circle;
    			while(circle){
    				circleNum++;
    				circle = circle->next;
    			}
    			if(circleNum != nextpathNum){
    				printf("Test if one of the circles should be worked first\n");
    			}
    		}
    		printf("Next Path is ambiguous --> End Path\n");
    		break;
    	}

    	// Break if no path is possible
    	if(!jPath){
    		printf("NO FOLLOWING PATH FOUND\n");
    		break;
    	}
	}
	free(pos_len);
//    free(uniPath);
}

struct scaffold_set* scaffold_init3(struct reads* reads){
	printf("CHECKPOINT: Run Scaffolding3\n");

	connectPathStats();

    int i;
    struct scaffold_set* aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
    aS->num = 0;
    aS->nummax = 1000;
    aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*aS->nummax);

//    char verbose = 0;
    char initRight = 0;

//    struct pathEdge* edge; // = (struct pathEdge*)malloc(sizeof(struct pathEdge));
	struct scaffEdge* scaffedge;

    prepPathsFlag();

    struct j_anno* j_annoL;
    struct j_anno* j_annoR;
    struct j_anno* j_anno;
    struct jPath* jPath;

    int jpos;
    int lastPath;
    int nextpathNum;
    char incomming = 0;

    int x = 10;

    while(x){
        for(i=1;i<pathsNum;i++){
        	if(paths[i].circfreq + paths[i].freq == 0) continue;
        	j_annoL = (struct j_anno*)reads[paths[i].leftJunction].annotation;
        	j_annoR = (struct j_anno*)reads[paths[i].rightJunction].annotation;
        	initRight = -1;
        	lastPath = i;
        	j_anno = NULL;
        	// Left Path side is death end
        	if(j_annoL->inDegree + j_annoL->outDegree == 1){
        		printf("\nFound Start point of a touring at path %i (end on left side) -> Pathfreq: %i\n",i,paths[i].circfreq + paths[i].freq);
        		initRight = 1;
        		j_anno = j_annoR;
        		jpos = paths[i].rightJunction;
        	}
        	// Right Path side is death end
        	if(j_annoR->inDegree + j_annoR->outDegree == 1){
        		printf("\nFound Start point of a touring at path %i (end on right side) -> Pathfreq: %i\n",i,paths[i].circfreq + paths[i].freq);
        		initRight = 0;
        		j_anno = j_annoL;
        		jpos = paths[i].leftJunction;
        	}
    		jPath = NULL;
    		nextpathNum = 0;
        	if(initRight >= 0){
        		jPath = j_anno->inEdge;
        		while(jPath){
        			if(jPath->pathID == lastPath){
        				printf("Next Path is OUTGOING\n");
        				jPath = j_anno->outEdge;
        				incomming = 0;
        				if(!jPath){
        					nextpathNum = -1;
        					break;
        				}
        				nextpathNum = j_anno->outDegree;
        				break;
        			}
        			jPath = jPath->next;
        		}
        		if(!jPath && nextpathNum != -1){
        			jPath = j_anno->inEdge;
        			incomming = 1;
    				if(!jPath){
    					nextpathNum = -1;
    				}
    				else{
    	    			printf("Next Path is INCOMING\n");
    	    			nextpathNum = j_anno->inDegree;
    				}
        		}
        		scaffedge = (struct scaffEdge*)malloc(sizeof(struct scaffEdge));
        		scaffedge->ID = i;
        		scaffedge->len = paths[i].len;
        		scaffedge->depth = 0;
        		scaffedge->next = NULL;
        		scaffedge->bridge = NULL;
        		aS->scaff[aS->num].first = scaffedge;
        		aS->scaff[aS->num].len = paths[i].len;
        		aS->scaff[aS->num].ID = aS->num;
        		aS->scaff[aS->num].scaffoldID = aS->num;

        		if(initRight){
        			aS->scaff[aS->num].startJunction = paths[i].leftJunction;
        			aS->scaff[aS->num].endJunction = paths[i].rightJunction;
        			scaffedge->targetJunction = paths[i].rightJunction;
        		}
        		else{
        			aS->scaff[aS->num].startJunction = paths[i].rightJunction;
        			aS->scaff[aS->num].endJunction = paths[i].leftJunction;
        			scaffedge->targetJunction = paths[i].leftJunction;
        		}

        		// delete set junction edges
        		scaffold_deleteJunctionEdge(i,reads,-1);

//        		paths[i].freq--;
    //    		if(nextpathNum == -1) aS->num++;
            	if(nextpathNum != -1 && jPath){
            		scaffold_uniqPath(initRight, incomming, i, nextpathNum, jpos, jPath, reads, aS);
            	}
            	aS->num++;
        	}
        }
        x--;
    }

    printf("Finished Scaffold:\n");
    int startJunction;
    for(int a=0;a<aS->num;a++){
		startJunction = aS->scaff[a].startJunction;
		printf("Scaffold: %i (len: %i bp) Type: %i\n",a,aS->scaff[a].len,aS->scaff[a].type);
		scaffedge = aS->scaff[a].first;
		printf(KRED"%i"KNRM,startJunction);
		while(scaffedge){
			printf(" -> "KGRN"%i"KNRM,scaffedge->ID);
			printf(" -> "KRED"%i"KNRM,scaffedge->targetJunction);
			scaffedge = scaffedge->next;
		}
		printf("\n");
    }

	return aS;
}

void scaffold_printfreqs(struct reads* reads, struct myovlList* G){
	printf("CHECKPOINT: Print Junction Connections\n");

	int i;
	struct j_anno* j_anno;
	struct jPath* jPath;
	int freq,circfreq;

    for(i=1; i <= G->V; i++){
    	if(G->read[i] && G->read[i]->flag == JUNCTION){
    		j_anno = (struct j_anno*)reads[i].annotation;
    		printf("Junction: %i (InDegree: %i / OutDegree: %i)\n",i,j_anno->inDegree,j_anno->outDegree);
    		printf("-> In-Edges:\n");
    		jPath = j_anno->inEdge;
    		while(jPath){
    			freq = paths[jPath->pathID].freq;
    			circfreq = paths[jPath->pathID].circfreq;
    			printf(" --> Path: %i (freq: %i, circ: %i)\n",jPath->pathID,freq,circfreq);
    			jPath = jPath->next;
    		}
    		printf("-> Out-Edges\n");
    		jPath = j_anno->outEdge;
    		while(jPath){
    			freq = paths[jPath->pathID].freq;
    			circfreq = paths[jPath->pathID].circfreq;
    			printf(" --> Path: %i (freq: %i, circ: %i)\n",jPath->pathID,freq,circfreq);
    			jPath = jPath->next;
    		}
    		printf("\n");
    	}
    }
//    exit(1);
}

static inline void deleteJunctionEdge(int pathID, struct reads* reads, int circleID){
	printf("CHECKPOINT: Delete junction edges to paths\n");
	// delete leftside edge
	char del = 0;
	char deledge = 0;
	if(paths[pathID].flag == 1) deledge = 1;
	struct contigCircle* circle;
	struct contigCircle* pcircle;
	int jun = paths[pathID].leftJunction;
	struct j_anno* j_anno = (struct j_anno*)reads[jun].annotation;
	struct jPath* jPath = j_anno->inEdge;
	struct jPath* pjPath = NULL;
	while(jPath){
		if(jPath->pathID == pathID){
			if(deledge){
				if(pjPath) pjPath->next = jPath->next;
				else j_anno->inEdge = jPath->next;
				free(jPath);
			}
			j_anno->inDegree--;
			del = 1;
			break;
		}
		pjPath = jPath;
		jPath = jPath->next;
	}
	if(!del){
		jPath = j_anno->outEdge;
		pjPath = NULL;
		while(jPath){
			if(jPath->pathID == pathID){
				if(deledge){
					if(pjPath) pjPath->next = jPath->next;
					else j_anno->outEdge = jPath->next;
					free(jPath);
				}
				j_anno->outDegree--;
				del = 1;
				break;
			}
			pjPath = jPath;
			jPath = jPath->next;
		}
	}
	if(!del){
		printf("Path of the junctions (%i) not found (InDegree: %i / Outdegree: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);
		exit(1);
	}

	printf("Junction: %i (in: %i, out: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);

	// rigth side
	jun = paths[pathID].rightJunction;
	j_anno = (struct j_anno*)reads[jun].annotation;
	jPath = j_anno->outEdge;
	pjPath = NULL;
	del = 0;
	while(jPath){
		if(jPath->pathID == pathID){
			if(deledge){
				if(pjPath) pjPath->next = jPath->next;
				else j_anno->outEdge = jPath->next;
				free(jPath);
			}
			j_anno->outDegree--;
			del = 1;
			break;
		}
		pjPath = jPath;
		jPath = jPath->next;
	}
	if(!del){
		jPath = j_anno->inEdge;
		pjPath = NULL;
		while(jPath){
			if(jPath->pathID == pathID){
				if(deledge){
					if(pjPath) pjPath->next = jPath->next;
					else j_anno->inEdge = jPath->next;
					free(jPath);
				}
				j_anno->inDegree--;
				del = 1;
				break;
			}
			pjPath = jPath;
			jPath = jPath->next;
		}
	}
	if(circleID>=0){
		circle = j_anno->circle;
		pcircle = NULL;
		while(circle){
			if(circle->ID == circleID && circle->pathID == pathID){
				if(pcircle) pcircle->next  = circle->next;
				else j_anno->circle = circle->next;
				free(circle);
				break;
			}
			pcircle = circle;
			circle = circle->next;
		}
	}
	if(circleID>=0) paths[pathID].flag--;
	printf("Junction: %i (in: %i, out: %i)\n",jun,j_anno->inDegree,j_anno->outDegree);
}

// Check if the next path is outgoing or incomming
static inline char nextpath_isOut(int junction, struct reads* reads, int pathID){

	struct j_anno* j_anno = (struct j_anno*)reads[junction].annotation;
	struct jPath* jPath;

	char out = -1;

	jPath = j_anno->inEdge;
	while(jPath){
		if(jPath->pathID == pathID){
			out = 1;
			break;
		}
		jPath = jPath->next;
	}
	jPath = j_anno->outEdge;
	while(jPath){
		if(jPath->pathID == pathID){
			if(out == 1) out = 2;
			else out = 0;
			break;
		}
		jPath = jPath->next;
	}

	return 0;
}

static int setUniqueNeighbor(struct contigScaff* side, int elem, int currentPath, int currentJ, struct reads* reads){
	char verbose = 1;

	struct j_anno* j_anno = (struct j_anno*)reads[currentJ].annotation;
	struct jPath* jPath;
	struct jPath* jPathcool;

	jPath = j_anno->inEdge;
	char found = 0;
	while(jPath){
		if(jPath->pathID == currentPath){
			found = 1;
			break;
		}
		jPath = jPath->next;
	}
	jPath = j_anno->outEdge;
	while(jPath){
		if(jPath->pathID == currentPath){
			if(found) return elem;
			else found = 2;
			break;
		}
		jPath = jPath->next;
	}
	if(found == 1) jPath = j_anno->outEdge;
	if(found == 2) jPath = j_anno->inEdge;
	int diffnum = 0;
	while(jPath){
		if(jPath && paths[jPath->pathID].flag > 0){
			diffnum++;
			jPathcool = jPath;
		}
		jPath = jPath->next;
	}

	if(diffnum == 1){
		if(verbose) printf("Set Unique Neighbor Path from %i to %i\n",currentPath,jPathcool->pathID);
		side[elem].ID = jPathcool->pathID;
		side[elem].sameside = 1;
		side[elem].bridge = NULL;
		elem++;
	}

	return elem;
}

// Combined for PE and SE data
struct scaffold_set* scaffold_init6(struct scaffold_set* aS, struct reads* reads, char bridging){
	printf("Checkpoint: Scaffold 6\n");
	prepPathsFlag();
    int i,j;
    if(!aS){
    	aS = (struct scaffold_set*)malloc(sizeof(struct scaffold_set));
        aS->num = 0;
        aS->nummax = 1000;
        aS->scaff = (struct scaffold*)malloc(sizeof(struct scaffold)*aS->nummax);
    }

    char verbose = 0;
    char stopit;

    int depth = 0;
    int len = 0;

    int ltemp;
    int rtemp;
    int lastID;
    char found;

    int lpos = 0;
    int rpos = 0;
    int lelem = 0;
    int relem = 0;
    int lmaxelem = 1000;
    int rmaxelem = 1000;
    int siblNum;
    int current;
    struct contigScaff* left = (struct contigScaff*)malloc(sizeof(struct contigScaff)*lmaxelem);
    struct contigScaff* right = (struct contigScaff*)malloc(sizeof(struct contigScaff)*rmaxelem);
    struct pathEdge* edge; // = (struct pathEdge*)malloc(sizeof(struct pathEdge));
    struct pathEdge* sibl;

    for(j=0;j<5;j++){
        for(i=1;i<pathsNum;i++){
        	if(paths[i].flag){
        		if(j || ((paths[i].scaffflag & 32) && (paths[i].scaffflag & 2))){
        			continue;
        		}
        		if(verbose)
        			printf("\t NEW SCAFFOLD --> Startpath: %i (Freq: %i)\n",i,paths[i].flag);
        		paths[i].flag--;
        		// initial left
        		lpos = 0;
        		rpos = 0;
        		lelem = 0;
        		relem = 0;
        		// TODO: Try to connect to the one unique solution to left
        		lelem = setUniqueNeighbor(left,lelem,i,paths[i].leftJunction,reads);
        		// TODO: Then following code, update lelem before
        		edge = paths[i].leftPath;
        		while(edge){
        			siblNum = 0;
        			if(paths[edge->ID].flag){
        				siblNum = 1;
        				sibl = edge;
        			}
        			while(edge->sibl){
        				if(verbose) printf("Sible Found: %i -> %i\n", edge->ID,edge->sibl->ID);
        				edge = edge->sibl;
        				if(paths[edge->ID].flag){
        					if(verbose) printf("is Flagged\n");
        					siblNum++;
        					if(siblNum == 1) sibl = edge;
        					if(siblNum > 1) break;
        				}
        			}
        			if(siblNum != 1) break;
        			else edge = sibl;

        			if(edge->junctionCon<0 || bridging){
        				if(edge->depth == lelem + 1){
                			if(verbose)
                				printf("(%i) Set left %i\n",i,edge->ID);
                			left[lelem].ID = edge->ID;
                			if(edge->targetJunction == paths[edge->ID].leftJunction) left[lelem].sameside = 1;
                			else left[lelem].sameside = 0;
                			if(edge->junctionCon>=0){
                				if(verbose) printf("\t--> is bridge\n");
                				left[lelem].bridge = edge;
                			}
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
        			}
        			else break;
//        			paths[edge->ID].flag--;
        			edge = edge->next;
        		}
        		// initial right
        		// TODO: Try to connect to the one unique solution to right
        		relem = setUniqueNeighbor(right,relem,i,paths[i].rightJunction,reads);
        		// TODO: Then following code, update relem before
        		edge = paths[i].rightPath;
        		while(edge){
        			siblNum = 0;
        			if(paths[edge->ID].flag){
        				siblNum = 1;
        				sibl = edge;
        			}
        			while(edge->sibl){
        				edge = edge->sibl;
        				if(paths[edge->ID].flag){
        					siblNum++;
        					if(siblNum == 1) sibl = edge;
        					if(siblNum > 1) break;
        				}
        			}
        			if(siblNum != 1) break;
        			else edge = sibl;
        			if(edge->junctionCon<0 || bridging){
        				if(edge->depth == relem + 1){
                			if(verbose)
                				printf("(%i) Set right %i\n",i,edge->ID);
                			right[relem].ID = edge->ID;
                			if(edge->targetJunction == paths[edge->ID].rightJunction) right[relem].sameside = 1;
                			else right[relem].sameside = 0;
                			if(edge->junctionCon>=0){
                				if(verbose) printf("\t--> is bridge\n");
                				right[relem].bridge = edge;
                			}
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
        			}
        			else break;
//        			paths[edge->ID].flag--;
        			edge = edge->next;
        		}
        		stopit = 0;
        		while((lpos<lelem || rpos<relem) && !stopit){
        			if(verbose) printf("Something was set, go deeper\n");
        			// left
        			while(lpos<lelem){
        				if(!paths[left[lpos].ID].flag){
        					stopit = 1;

        					break;
        				}
        				current = left[lpos].ID;
        				paths[current].flag--;
        				// left -> left
        				if(verbose) printf("LL\n");
        				if(left[lpos].sameside) edge = paths[current].leftPath;
        				else edge = paths[current].rightPath;
        				while(edge){
                			siblNum = 0;
                			if(paths[edge->ID].flag){
                				siblNum = 1;
                				sibl = edge;
                			}
                			while(edge->sibl){
                				edge = edge->sibl;
                				if(paths[edge->ID].flag){
                					siblNum++;
                					if(siblNum == 1) sibl = edge;
                					if(siblNum > 1) break;
                				}
                			}
                			if(siblNum != 1) break;
                			else edge = sibl;

                			if(edge->junctionCon<0 || bridging){
            					if(edge->depth + lpos == lelem){
            						if(verbose) printf("LL: (%i) (lpos: %i, lelem: %i) Set left left %i (target: %i)\n",i,lpos,lelem,edge->ID,edge->targetJunction);
    //        						if(verbose) printf("(%i) (lpos: %i) Set left left %i (target: %i)\n",i,lpos,edge->ID,edge->targetJunction);
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
    //        		    			paths[edge->ID].flag--;
            					}
                			}
                			else break;

        					edge = edge->next;
        				}
        				// left -> right
        				if(verbose) printf("LR: %i\n",current);
        				if(left[lpos].sameside) edge = paths[current].rightPath;
        				else edge = paths[current].leftPath;
        				ltemp = lpos-1;
        				rtemp = 0;
        				while(edge){
        					if(edge->depth < lpos + relem +2){
        						found = 0;
        						if(ltemp>=0){
        							lastID = left[ltemp].ID;
        							ltemp--;
        						}
        						else if(ltemp == -1){
        							lastID = i;
        							ltemp--;
        						}
        						else{
        							lastID = right[rtemp].ID;
        							rtemp++;
        						}
        						if(verbose) printf("LastId: %i\n",lastID);
        						while(edge){
        							if(edge->ID == lastID){
        								if(verbose) printf("EdgeID: %i\n",edge->ID);
        								found = 1;
        								edge = edge->next;
        								break;
        							}
        							edge = edge->sibl;
        						}
        						if(found) continue;
        						else{
        							if(verbose) printf("Correct Backpath not found\n");
        							break;
        							exit(1);
        						}
        					}
                			siblNum = 0;
                			if(paths[edge->ID].flag){
                				siblNum = 1;
                				sibl = edge;
                			}
                			while(edge->sibl){
                				edge = edge->sibl;
                				if(paths[edge->ID].flag){
                					siblNum++;
                					if(siblNum == 1) sibl = edge;
                					if(siblNum > 1) break;
                				}
                			}
                			if(siblNum != 1) break;
                			else edge = sibl;
                			if(edge->junctionCon<0 || bridging){
            					if(edge->depth == lpos + relem +2){
            						if(verbose) printf("LR: (%i) (rpos: %i, relem: %i) Set right %i (target: %i)\n",i,rpos,relem,edge->ID,edge->targetJunction);
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
    //        		    			paths[edge->ID].flag--;
            					}
                			}
                			else break;
        					edge = edge->next;
        				}
        				lpos++;
        			}
        			// right
        			while(rpos < relem){
        				if(!paths[right[rpos].ID].flag){
        					stopit = 1;
        					break;
        				}
        				current = right[rpos].ID;
        				paths[current].flag--;
        				// right -> right
        				if(verbose) printf("RR\n");
        				if(right[rpos].sameside) edge = paths[current].rightPath;
        				else edge = paths[current].leftPath;
        				while(edge){
                			siblNum = 0;
                			if(paths[edge->ID].flag){
                				siblNum = 1;
                				sibl = edge;
                			}
                			while(edge->sibl){
                				edge = edge->sibl;
                				if(paths[edge->ID].flag){
                					siblNum++;
                					if(siblNum == 1) sibl = edge;
                					if(siblNum > 1) break;
                				}
                			}
                			if(siblNum != 1) break;
                			else edge = sibl;
                			if(edge->junctionCon<0 || bridging){
            					if(edge->depth + rpos == relem){
            						if(verbose) printf("RR: (%i) (rpos: %i, relem: %i) Set right %i (target: %i)\n",i,rpos,relem,edge->ID,edge->targetJunction);
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
    //        		    			paths[edge->ID].flag--;
            					}
                			}
                			else break;
        					edge = edge->next;
        				}
        				// right -> left
        				if(verbose) printf("RL: %i\n",current);
        				if(right[rpos].sameside) edge = paths[current].leftPath;
        				else edge = paths[current].rightPath;
        				ltemp = 0;
        				rtemp = rpos-1;
        				while(edge){
        					if(edge->depth < rpos + lelem +2){
        						found = 0;
        						if(rtemp>=0){
        							lastID = right[rtemp].ID;
        							rtemp--;
        						}
        						else if(rtemp == -1){
        							lastID = i;
        							rtemp--;
        						}
        						else{
        							lastID = left[ltemp].ID;
        							ltemp++;
        						}
        						if(verbose) printf("LastId: %i\n",lastID);
        						while(edge){
        							if(verbose) printf("EdgeID: %i\n",edge->ID);
        							if(edge->ID == lastID){
        								if(verbose) printf("-> Correct\n");
        								found = 1;
        								edge = edge->next;
        								break;
        							}
        							edge = edge->sibl;
        						}
        						if(found) continue;
        						else{
        							printf("Correct Backpath not found\n");
        							break;
        							exit(1);
        						}
        					}
                			siblNum = 0;
                			if(paths[edge->ID].flag){
                				siblNum = 1;
                				sibl = edge;
                			}
                			while(edge->sibl){
                				edge = edge->sibl;
                				if(paths[edge->ID].flag){
                					siblNum++;
                					if(siblNum == 1) sibl = edge;
                					if(siblNum > 1) break;
                				}
                			}
                			if(siblNum != 1) break;
                			else edge = sibl;
                			if(edge->junctionCon<0 || bridging){
            					if(edge->depth == rpos + lelem +2){
            						if(verbose) printf("RL: (%i) (lpos: %i, lelem: %i) Set left left %i (target: %i)\n",i,lpos,lelem,edge->ID,edge->targetJunction);
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
    //        		    			paths[edge->ID].flag--;
            					}
                			}
                			else break;
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

