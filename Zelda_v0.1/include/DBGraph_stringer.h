/*
 ============================================================================
 Name        : DBGraph_stringer.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Heart of Zelda: Converts the dBG to a String Graph by calculating
 	 	 	   transitively irreducible overlaps from the dBG.
 ============================================================================
 */

#ifndef DBGRAPH_STRINGER_H_
#define DBGRAPH_STRINGER_H_

#include "DBGraph.h"

/********  OVERLAP / UNITIG GRAPH OBJECT ********/

#define VERT(v) ( ((v) + 1) >> 1 )
#define MATE(v) ( ((v) % 2) ? (v) + 1 : (v) - 1 )

#define WIDOWED   0   					/* No overlaps with any other read                          */
#define CONTAINED 1   					/* Contained by at least one other read                     */
#define STRAY     2   					/* Not contained, but all overlaps are with contained reads */
#define PROPER    3   					/* None of the other states                                 */
#define JUNCTION  4   					/* After reduction, node has an in or out degree != 1       */

struct node{							// Only used in struct partGraph, equivalent to Kannte
	int ID;								// Node-ID
	int layer;							// Number of edges to go to come to this node (0->1->2->3->4)
	int depth;							// Global position of the node
	int flag;							// Important to know which edges are allowed to go from the given node
	struct edge *head;					// Edges to the parents
	struct edge *tail;					// Edges to the children
};

struct partGraph{						// Used in Stringer and provides a small snapshot downwards a given node
	int V;								// Num of Nodes
	struct node *array;					// Nodes somehow different to Kannten important for the partGraph
};

struct dest{							// Pointer is set to the first proper read outgoing of a junction read, describing the end of a path
	int ID;								// ID the junction read at the end of a unique path
	int len;							// length of the unique path in basepairs
	int flag;							// is 0 if the path was not toured at this moment, 1 otherwise. No second touring
	struct bread* counterbread;			// Link to corresponding destination of the back path
//	struct path* path;
	int pathID;
};

struct bread{							// List of transitively reduced breads overlapping the aread
	char dir;							// Direction relatively to the aread
	int ID;								// ID of the bread
	int overhang;						// Difference between the "end" of tha aread to the "end" of the bread
	char flag;							//// Only for testing
	char sideflag;						// Readside overhang
	struct bread* next;					// Next bread of this particular aread
	struct dest* dest;					// If the bread is the first read of a path to a junction the dest is set, NULL otherwise
};

struct aread{							// Aread if exists, NULL otherwise (aread have had to many error, so it wasn't connected to the graph)
	int length;							// Read length
	char dir;							// Direction of the read relatively to the graph -> critical point (a read could by used more than one time in different directions)
	char flag;							// Read-overlap status (values 0-4 (widowed, contained, stray, proper, junction))
	int stk;							// !!! Only for viz to debug !!! Start node of the read
	int endk;							// !!! Only for viz to debug !!! End node of the read
	struct bread *first;				// First bread overlap if exists, NULL otherwise
};

struct myovlList{						// Struct contains all valid overlap informations (aread array is sorted by aread ID)
	int V;								// Number of reads
	struct aread **read;				// Array of the areads (position of aread in the array is the ID of the reads in the database)
};

// copy of reducer structs
struct readend{
    int   top; 							/* Index of last edge on adjacency sub-array 							*/
    int   mark;     					/* Initialized to 0.  Used by lots of steps as a mark field but always
											restored so that it is the index of its node in the string graph if
                     	 	 	 	 		it is a branch point, and 0 otherwise.   							*/
};

struct overhang{
	int 			nreads;				/*	# of reads supporting this overhang						*/
    int           	target;   			/*  Index of read end at the head of the edge 				*/
    int           	length;   			/*  Length of the overhang to this read end   				*/
    unsigned char	reduced;  			/*  Flag that is set if the edge is transitively reducible 	*/
    short         	divergence;  		/*  Alignment divergence percentage 						*/
};

// Final String Graph

struct string_graph{
	int             nverts;   			// # of reads (only junctions left? -> No, all (+ widowed and contained))
	int             nends;    			// # of read ends = 2*nverts
	int             nedges;   			// # of overhangs = sides[nends].top+1
    struct readend  *side; 		    	// [1..nends] readend records, side[0].top = -1 to simplify boundary
    struct overhang *edge;     			// [0..nedges-1] overhang records
    unsigned char  *status;			// [1..nverts] status bytes
    int             *length;   			// [1..nverts] length of reads
    int             *ID;   				// [1..nverts] Read ID
};

struct overlap_graph{
    int             nverts;   			// # of reads
    int             nends;    			// # of read ends = 2*nverts
    int             nedges;   			// # of overhangs = sides[nends].top+1
    struct readend  *side; 		    	// [1..nends] readend records, side[0].top = -1 to simplify boundary
    struct overhang *edge;     			// [0..nedges-1] overhang records
    unsigned char  *status;			// [1..nverts] status bytes
    int             *length;   			// [1..nverts] length of reads
    int             *offset;   			// [1..nverts] left trim offset
};

// String graph construction

int startEnd2(int i,int *stp, int *endp, int *stk, int *endk, struct ReadNode *rnode,struct partGraph *childReadIds);

int lockKids(int i,struct partGraph *childReadIds);

int case3Lock(int aendk, struct partGraph *childReadIds);

void downstreamTree(int i, struct partGraph *childReadIds);

int stringer2(struct myovlList*);

int stringCase2_2(int i, int astp, int aendp, struct ReadNode *nextrNode,struct ReadNode *rNode,struct partGraph *childReadIds, struct myovlList *ovlGraph);

int stringCase3_2(int i, int astp,int aendp, int aendk, struct ReadNode *nextrNode,struct ReadNode *rNode,struct partGraph *childReadIds, struct myovlList *ovlGraph);

int stringCase4_2(int astp, int aendp, struct ReadNode *rNode,struct partGraph *childReadIds, struct myovlList *ovlGraph);

int stringCase4_3(int astp, int aendp, struct ReadNode *rNode,struct partGraph *childReadIds, struct myovlList *ovlGraph);

int stringCase4_4(int astp, int aendp, struct ReadNode *rNode,struct partGraph *childReadIds,struct myovlList *ovlGraph);

// String graph processing, refinement and output

void countRemainingNodes();

struct myovlList* initOVLgraph(int);

struct string_graph* initStringGraph(struct myovlList*, char*);

void printStringGraph(struct myovlList*, char*);

void printOVLgraph(struct myovlList*, int , char*);

struct string_graph* catOVLgraph(struct myovlList*, char*);

void ovlToString(struct myovlList *O, struct string_graph *S, char* pathAssembly);

void print_string_graph_list(struct string_graph* G, char* label);

void printReducedOVLgraph(struct myovlList *ovlGraph,int dir, char *ovlPath);

void print_overlap_graph_dot(struct string_graph* G, const char* pcFile, const char* title);

char* vtx_name(int v, int width);

// Based on the ideas of Gene
// Look-ups in other direction
int stringer3(struct myovlList *ovlGraph);

void upstreamTree(int i, struct partGraph *childReadIds);

int endStart(int i,int *stp, int *endp, int *stk, int *endk, struct ReadNode *rnode, struct partGraph *childReadIds,int* reallen);

void tag_A_Contained(struct myovlList *ovlGraph);

void freeMyOvlList(struct myovlList* G, struct string_graph* S);

#endif /* DBGRAPH_STRINGER_H_ */


