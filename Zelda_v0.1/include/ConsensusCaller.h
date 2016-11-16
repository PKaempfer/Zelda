/*
 * ConsensusCaller.h
 *
 *  Created on: Feb 26, 2016
 *      Author: kaempfpp
 */


#ifndef CONSENSUS_H_
#define CONSENSUS_H_

#include "DBGraph_stringer.h"
#include "readDB.h"
//#include "../poaV2/poa.h"

#define GAP_PENALTY  -1


// Substitution Matrix 1 		  A, C, G, T, -
//static const int SM1[5][5] = 	{
//								 {1, 0, 0, 0, 0},	// A
//								 {0, 1, 0, 0, 0},	// C
//								 {0, 0, 1, 0, 0},	// G
//								 {0, 0, 0, 1, 0},	// T
//								 {0, 0, 0, 0, 0}	// -
//};
//static const int SM1[5][5] = 	{
//								 {2, 1, 1, 1, 1},	// A
//								 {1, 2, 1, 1, 1},	// C
//								 {1, 1, 2, 1, 1},	// G
//								 {1, 1, 1, 2, 1},	// T
//								 {1, 1, 1, 1, 2}	// -
//};
static const int SM1[5][5] = 	{
								 {1, -1, -1, -1, -1},	// A
								 {-1, 1, -1, -1, -1},	// C
								 {-1, -1, 1, -1, -1},	// G
								 {-1, -1, -1, 1, -1},	// T
								 {-1, -1, -1, -1, -1}	// -
};

struct contig{
	int len;
	int stReadID;
	int endReadID;
	char* seq;
};

struct contigList{
	int num;
	int maxnum;
	struct contig* contig;
};

// POA structs -> See ../poaV2/poa.h
// Make my own data structure

struct LetterSource_S {
	int iseq;						/** index of the sequence, referencing the source_seq[] array*/
	uint32_t ipos;					/** index of the corresponding position in that sequence */
	struct LetterSource_S *next; 	/** next node in the linked list */
};

struct LetterEdge{
	uint32_t dest;
	unsigned char counter;
	struct LetterEdge* next;
};

/** Structure for storing individual Letters*/
struct Letter_T {
	struct LetterEdge* left; 		/** ADJACENT LETTER(S) TO THE LEFT */
	struct LetterEdge* right; 		/** ADJACENT LETTER(S) TO THE RIGHT */
	struct LetterSource_S source;	/** SOURCE SEQ POSITION(S) */
	struct Letter_T* align_ring; 	/** CIRCULAR LIST OF ALIGNED POSITIONS */
	uint16_t counter;				/** Number reads supporting this letter (number of sources)  */
	int* ml;						/** Line of Alignment matrix between the local PO graph and the new sequence*/
	int score;						/** SCORE FOR BALANCING PARTIAL ORDER EFFECTS ON MATRIX NEUTRALITY */
	char letter;					/** THE ACTUAL RESIDUE CODE! */
	char junction;
	char vFlag;
};

/** holder for an LPO sequence, its letters, and associated information */
struct Sequence {/** */
	int length;						/** Backbone/Ref length */
	struct LetterEdge startLetter;	/** First Letter of the contig graph1*/
	uint32_t readleft;				/** */
	uint32_t readright;				/** */
	char* title;					/** */
	char* sequence;					/** Title of the contig*/
	char* name;						/** Name of Contig */
	int nsource_seq;				/** Number of sequences in representing this PO graph */
//	LPOSourceInfo_T *source_seq;	/** */
};

struct POG{
	uint32_t contigNum;
	uint32_t maxNum;
	struct Sequence* contig;
};

struct pairAlign{
	char* refSeq;
	char* readSeq;
	struct Letter_T* current;
	int len;
	int j;
};

extern struct Letter_T* Letters;
extern uint32_t numNodes;
// For ConsensusCaller2 -> Temporary delete when merged
extern uint32_t maxNumNodes;
extern int **alMatrix;
extern int *alMatrix_Best;
extern struct Letter_T** alMatrix_Letter;
extern struct timespec ts_start;
extern struct timespec ts_finish;
extern long sumMatrix;
extern long sumTrace;

void buildBackBone(struct myovlList* ovlgraph, struct string_graph* S, struct reads* reads);

void buildBackBone2(struct myovlList* G, struct reads* reads);

void buildBackBone3(struct myovlList* G, struct reads* reads);

struct contigList* realBackbone(struct myovlList* G, struct reads* reads);

struct contigList* realBackbone2(struct myovlList* G, struct reads* reads);

//Sequence_T* make_poa(struct myovlList* G, struct read* reads);

void poa_part_toDot(char* dotFile, struct Sequence* contig);

void poa_catBackbone2(struct Sequence* contig, struct myovlList *G, struct reads* read, char* seq, int leftID, int rightID);

void testFunct();

// 1 align Function
int poa_align_prepro(struct Sequence* contig, int len, int overhang);

void poa_heuristic_align(struct Sequence* contig, struct reads* read, char* seq, char backbone, int insNum);

void poa_heuristic_align2(struct Sequence* contig, struct reads* read, char* seq, char backbone, int insNum, int overhang);

struct POG* make_poa(struct myovlList* G, struct reads* reads);

struct POG* make_poaScaff(struct myovlList* G, struct reads* reads, char scaffolding);

void poa_toDot(char* dotFile);

void poa_consensus(struct Sequence* contig);

void poa_printContigs(struct POG* pog, char* contigFile);

void free_POG(struct POG* contigs_pog);

#endif
