/*
 ============================================================================
 Name        : DBGraph_scaffold.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Scaffolding unit of the DBSGA
 ============================================================================
 */

#ifndef DBGRAPH_SCAFFOLD_H_
#define DBGRAPH_SCAFFOLD_H_

#include "readDB.h"
#include "DBGraph_stringer.h"
#include "FileReader.h"

extern struct path* paths;
extern int pathsNum;
extern uint64_t pathsTotBase;
extern uint64_t pathsTotCov;
extern uint64_t pathsTotLen;

struct contigCircle{
	int ID;						// ID of the circle
	int dest;					// Junction readID of the next node of the circle
	int pathID;					// ID of the following path, as part of the circle
	char vFlag;					// visited Flag - Edges not allowed to visit 2 times
	struct contigCircle* next;	// Next circle, if the node is origin of more than one circles
};

struct pathEdge{				// PE - Connection of the path to the following paths on the way in concordance to the read pairings (order of the pathedges should be consistent with the depth of the path)
	int ID;						// ID of the connected Path
	int counter;				// Number of read pairs supporting this path connection
	int targetcounter;			// Number of read pairs targeting this path directly
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
	uint32_t cov;				// Contig frequency (Means initially the coverage, to estimate the frequency the contig appears in the assembly)
	uint32_t freq;				// Number of allowed repeat usage of this contig path (without circle usage)
	uint32_t circfreq;			// Circle frequency, not used in scaffold path finding
	char pathdir;				// 0 if the path was created over bread->sideflag == 0 from leftJunction; 1 if bread->sideflag was 1 from leftJunction
	int leftJunction;			// left Junction read ID
	struct pathEdge* leftPath;	// adjacency list of possible following paths with counter of read pairs supporting this connection
	int rightJunction;			// right Junction read ID
	struct pathEdge* rightPath;	// adjacency list of possible following paths with counter of read pairs supporting this connection
	int component;				// Component ID, the path is a part of
	char flag;					// nothing at the moment - just there
	char scaffflag;				// first 4 bit left, last 4 bit right (0 - End, 1 - Proper, 2 - Ambiguous)
};

struct pc_anno{ 				// Annotation for non junction reads (proper or contained)
	int pathID;					// ID of path, the read belongs to
	int lJunctionDist;			// Distance to the left junction in Path
	int rJunctionDist;			// Distance to the right junction in Path
};

struct jPath{
	uint32_t pathID;
	uint32_t sJID;
	uint32_t eJID;
	struct jPath* next;
};

struct j_anno{		 			// Annotation for Junction reads
	int inDegree;				// number of paths into the junction
	int outDegree;				// number of paths out of the junction
	char vFlag;					// Visited Flag
	struct contigCircle* circle; // Circles in the Contig Graph, not investigated in the Scaffold touring but substituted when one node is visited
	struct jPath* inEdge;		// Paths out of the Junction
	struct jPath* outEdge;		// Paths into the Junction
//	int** counts;				// Number of read pair support from in to out path
};

struct scaffEdge{
	int ID;
	int len;
	int depth;
	int targetJunction;
	struct pathEdge* bridge;
	struct scaffEdge* next;
};

struct scaffold{
	int ID;
	int len;
	uint32_t estim_Cov;			// Estimated Coverage
	float true_Cov;				// True Coverage
	uint32_t testVar_delete;
	char type; 					// 0 -> Scaffold, 1 -> Singleton (unscaffolded contig)
	int startJunction;
	int endJunction;
	int scaffoldID;
	struct scaffEdge* first;
	int next;
};

struct contigScaff{
	int ID;
	char sameside;				// 1 - if same side as start contig, 0 - otherwise
	struct pathEdge* bridge;	// gets entire path information if the edge is a bridge, NULL otherwise
};

struct scaffold_set{
	struct scaffold* scaff;
	int numMinLen;
	int num;
	int numbridge;
	int nummax;
};

struct scaffold_set* contigs_init(struct myovlList* G); // , struct reads* reads

void initScaff(struct myovlList* G, struct reads* reads);

void free_schaffoldSet(struct scaffold_set* aS);

void readTouring(struct myovlList* G, struct readFiles* files, struct reads* reads);

void scaffGraphDot(struct myovlList* G, struct reads* reads, char* dotfile);

#endif /* DBGRAPH_SCAFFOLD_H_ */
