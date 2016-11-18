/*
 ============================================================================
 Name        : DBSGA.c
 Author      : Kämpfer, Philipp
 Version     :
 Copyright   : GPLv3 (GNU general public licence)
 Description : Scaffolding unit of the DBSGA
 ============================================================================
 */

#ifndef DBGRAPH_SCAFFOLD_H_
#define DBGRAPH_SCAFFOLD_H_

#include "readDB.h"
#include "DBGraph_stringer.h"
#include "FileReader.h"

/** Relative orientation of the reads in the graph to each other */
//#define FR 0 	//  --->   <--- (forward - reverse) -> Default
//#define RF 1	//  <---   ---> (reverse - forward)
//#define FF 2	//  --->   ---> (forward - forward) or
//#define RR 2	// 	<---   <--- (reverse - reverse)
//// FF =^ RR
//
//static const char* peOri[3]={"FR","RF","FF"};

extern struct path* paths;
extern int pathsNum;

struct pc_anno{ 				// Annotation for non junction reads (proper or contained)
	int pathID;					// ID of path, the read belongs to
	int lJunctionDist;			// Distance to the left junction in Path
	int rJunctionDist;			// Distance to the right junction in Path
};

struct j_anno{		 			// Annotation for Junction reads
//	int inDegree;				// number of paths into the junction
//	int outDegree;				// number of paths out of the junction
	char vFlag;					// Visited Flag
//	int** counts;				// Number of read pair support from in to out path
};

struct pathEdge{				// PE - Connection of the path to the following paths on the way in concordance to the read pairings (order of the pathedges should be consistent with the depth of the path)
	int ID;						// ID of the connected Path
	int counter;				// Number of read pairs supporting this path connection
	int depth;					// distance of the connected path to the origin path
	int targetJunction;			// .--->.<---. (readR1 on first path, target would be the last dot (junction))
	int estLen;					// Estimated Length of the path edge, can span several paths, than the sum of length or be a virtual path, than estimated by insert size
	struct pathEdge* next;		// next path connection
	struct pathEdge* sibl;		// Sibling -> alternative path in same depth. Forms a closed ring
	char junctionCon;			// Path Connection Junction (Left or right Junction of the path): 0: L; 1: R
//	char junctionDir;			// Connection Junction Orientation: 0: FF --> -->; 1: FR --> <--; 2: RF <-- -->; 3: RR <-- <--
};

struct path{
	int ID;						// ID of the path -> not necessary, because identical with the index
	int len;					// length of the path/contig
	char pathdir;				// 0 if the path was created over bread->sideflag == 0 from leftJunction; 1 if bread->sideflag was 1 from leftJunction
	int leftJunction;			// left Junction read ID
	struct pathEdge* leftPath;	// adjacency list of possible following paths with counter of read pairs supporting this connection
	int rightJunction;			// right Junction read ID
	struct pathEdge* rightPath;	// adjacency list of possible following paths with counter of read pairs supporting this connection
	int component;				// Component ID, the path is a part of
	char flag;					// nothing at the moment - just there
};

struct scaffEdge{
	int ID;
	int len;
	int depth;
	int targetJunction;
	struct scaffEdge* next;
};

struct scaffold{
	int ID;
	int len;
	char type; 					// 0 -> Scaffold, 1 -> Singleton (unscaffolded contig)
	int startJunction;
	int endJunction;
	struct scaffEdge* first;
};

struct scaffold_set{
	struct scaffold* scaff;
	int num;
	int nummax;
};

struct scaffold_set* contigs_init(struct myovlList* G); // , struct reads* reads

void initScaff(struct myovlList* G, struct reads* reads);

void free_schaffoldSet(struct scaffold_set* aS);

void readTouring(struct myovlList* G, struct readFiles* files, struct reads* reads);

void scaffGraphDot(struct myovlList* G, struct reads* reads, char* dotfile);

#endif /* DBGRAPH_SCAFFOLD_H_ */
