/*
 * DBgraph.h
 *
 *  Created on: Oct 13, 2014
 *      Author: kaempfpp
 */

#ifndef DBGRAPH_H_
#define DBGRAPH_H_

#include "kmer.h"
#include "uthash.h"
//#include "DBGraph_reduced.h"
//#include "DBGraph_error.h"

#define MINCOMP 100
#define MAX_DEPTH 200

extern struct Reads *allReads;
extern int minovlLen;
extern uint32_t tabindex;

struct Reads{							// Reads with connection to the Kannten; part of the ReadNodes
	struct KannteNode *headkannte;		// Pointer to Kannten the reads are contained
	int ID;								// Read ID and hash-key in the order the reads came in during the file reading
	UT_hash_handle hhb;					// hash handle variable
};

struct ReadNode{						// Start- or endpoint of a read; part and listed by the Kannte
	int pos;							// Pos in the Kannte
	char dir;							// Direction of the read in relation to the graph direction
	char flag;
	struct Reads* read;					// Connection to the hole read information
	struct ReadNode* next;				// The next ReadNode
};

struct edge{							// Edge between the Kannten
	int dest;							// The ID of the Kannte the edge is pointing on
	struct edge* next;					// tThe next Edge of the given Kannte
};

struct KannteNode{						// The link from reads to the corresponding Kannten
	int dest;							// Kannten ID
	int pos;							// Pos of the Read link in the Kannte
	struct ReadNode* ReadNode;			// The corresponding ReadNode of this KannteNode
	struct KannteNode* next;			// The next KannteNode
};

struct readEdge{
	int KannteID;
	struct readEdge* next;
};

struct Kannte{
	int len;							// Number of containing kmers of this Kannte
	int minpos;							// Mintag Pos
	int maxpos;							// Maxtag Pos
	struct edge* head;					// Edges to the parents
	struct edge* tail;					// Edges to the children
	struct ReadNode* headread;			// Pointer to the read start- and end positions
	struct readEdge* readedge;			// Edges to the nodes the reads are connected with
};

struct redGraph{						// The reduced graph
    int V;								// Number of the contained nodes in the entire graph
    char* vFlag;						// Flag to mark visited nodes
    struct Kannte* array;				// the List of the nodes
};

struct tempHash{						// For real time-check if unreduced junction is already in reduced graph;
	int oldID;							// ID of the initial DeBruijn graph; the hash value
	int numIDs;							// ID of the new reduced Graph
	int* newID; 						// all new edges starting in this old node
	UT_hash_handle hhb;					// hash handle variable
};

struct AdjListNode{
    int dest;
    char trans;
    struct AdjListNode* next;
};

struct LinkListNode{
	readID ID;
	struct LinkListNode *next;
};

struct AdjList{
    struct AdjListNode *head;  // pointer to head node of list (node parents)
    struct AdjListNode	*tail;	// pointer to head node of list (node childs)
    struct LinkListNode *Link;
};

struct Graph{
    int V;
    char *vFlag;
    struct AdjList* array;
};

struct seekOvl{
	int oldlen;
	int newlen;
	int *childs;
	int *parents;
	int *cbefpos;
	int *pbefpos;
};

extern struct Graph *graph;
extern struct redGraph* redGraph;

void writePartialDot(int ori, int width);

void writeDot();

void printGraph();

void countKmers();

void printSorrounding(int ori);

struct LinkListNode* newLinkListNode(readID ID);

void addEdge(int src, int dest, char base, char base2);

void addLinks(readID *ID, int index);

void createGraph(int V);

void freeGraph();

void rekFindChild(struct hashTable *u);

void hashToTabDFS();

void travDFS(struct hashTable* s);

void travDownDFS(struct hashTable* s);

void deleteComp(int start);

void hashToTabBFS();

char isChild(int ori, int chi);

char isParent(int ori, int chi);

void travUpBFS(struct hashTable* s);

void travDownBFS(struct hashTable* s);

// Adjacency list bit operation functions

static inline char toBinTrans(char base){
	return (char)(codes[(int)base]);
}

static inline char toCharTrans(char base){
	return (rev_codes[base & TRANS_MASK]);
}

static inline void setVFlag(int s){
	graph->vFlag[s] |= 0x80;
}

static inline void delVFlag(int s){
	graph->vFlag[s] |= 0x7F;
}

static inline int isVFlagged(int s){
	if(__builtin_clz(graph->vFlag[s])) return 0;
	else return 1;
}

void printReads();

char* revRead(char *read);

void revReadSt(char *read, char* revread);

#endif /* DBGRAPH_H_ */
