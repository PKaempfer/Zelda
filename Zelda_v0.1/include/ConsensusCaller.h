/*
 ============================================================================
 Name        : ConsensusCaller.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Consensus Caller utilizes a re-implementation of POA - Algorithms
 	 	 	   (Lee et al.,2002) for multiple sequence alignments to create a
 	 	 	   Layout
 ============================================================================
 */

#ifndef CONSENSUS_H_
#define CONSENSUS_H_

#include "DBGraph_stringer.h"
#include "readDB.h"

#define GAP_PENALTY  -1
#define MATRIX_MAX_BR 100
//#define TIMEM
//#define RUSAGE

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
static const char SM1[5][5] = 	{
								 {1, -1, -1, -1, -1},	// A
								 {-1, 1, -1, -1, -1},	// C
								 {-1, -1, 1, -1, -1},	// G
								 {-1, -1, -1, 1, -1},	// T
								 {-1, -1, -1, -1, -1}	// -
};

static const char* const varType[] = {"snp\0", "mnp\0", "ins\0", "del\0", "complex\0"};

//int16_t data2[maxReadLen*MATRIX_MAX_BR+1][maxReadLen+1];

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

//struct LetterSource_S {
//	int iseq;						/** index of the sequence, referencing the source_seq[] array*/
//	uint32_t ipos;					/** index of the corresponding position in that sequence */
//	struct LetterSource_S *next; 	/** next node in the linked list */
//};

struct LetterEdge{
	uint32_t dest;
	unsigned char counter;
	char vFlag;
	struct LetterEdge* next;
};

/** Structure for storing individual Letters*/
struct Letter_T {
	struct LetterEdge* left; 		/** ADJACENT LETTER(S) TO THE LEFT */
	struct LetterEdge* right; 		/** ADJACENT LETTER(S) TO THE RIGHT */
//	struct LetterSource_S source;	/** SOURCE SEQ POSITION(S) */
	struct Letter_T* align_ring; 	/** CIRCULAR LIST OF ALIGNED POSITIONS */
	uint16_t counter;				/** Number reads supporting this letter (number of sources)  */
	int16_t* ml;					/** Line of Alignment matrix between the local PO graph and the new sequence*/
	int score;						/** SCORE FOR BALANCING PARTIAL ORDER EFFECTS ON MATRIX NEUTRALITY */
	unsigned char letter;			/** THE ACTUAL RESIDUE CODE! */
	char junction;					/** Number of Passing Paths through this node, continue aligning behind when last path reached the node*/
	char vFlag;						/** Flags the Nodes when used in Consensus Calling, Tags the Consensus Path*/
};

struct Variation{
	uint32_t pos;					/** Variant Position relatively to reference sequence */
	char* refSeq;					/** Reference sequence at the given pos */
	char* altSeq;					/** Reference sequence at the given pos */
	uint16_t dp;					/** Totol Read Depth at this position */
	uint16_t ao;					/** Number of reads supporting the alternative */
	uint16_t ro;					/** Number of reads supporting the reference */
	uint16_t len;					/** Length of truly different bases */
	unsigned char type;			/** Type of alternative (snp, mnp, ins, del, complex) */
	struct Variation* next;			/** Next Variation Position*/
};

struct sequenceEdge{
	int nextScaff;					/** ID of the next scaffold over a bridge*/
	int insertLen;					/** Bridge Length insert of Ns */
	char ori;						/** 0 -> same, 1 reverse, 2 reverse complement */
};

struct contigPart{
	uint32_t startpos;
	uint32_t endpos;
	uint32_t startread;
	uint32_t endread;
	struct contigPart* next;
};

/** holder for an LPO sequence, its letters, and associated information */
struct Sequence{
	int length;						/** Backbone/Ref length */
	struct LetterEdge startLetter;	/** First Letter of the contig graph1*/
	uint32_t readleft;				/** */
	uint32_t readright;				/** */
	char* title;					/** */
	char* sequence;					/** Title of the contig*/
	char* name;						/** Name of Contig */
	float avgCov;					/** Average Coverage of Contig*/
	int nsource_seq;				/** Number of sequences in representing this PO graph */
	char vflag;						/** 0 if free to use, 1 go ahead */
	struct Variation* var;			/** Linked List of all Variations of the Contig/Scaffold */
	struct Variation* lastvar;		/** Last Variation of the Contig/Scaffold */
	struct sequenceEdge* seqEdge;	/** Edge to the next part of the scaffold connected by a bridge*/
	struct contigPart* first;		/** First contig of a linked list with all contigs in this scaffold*/
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
extern int gdepth;
extern struct Letter_T* Letters;
extern uint32_t numNodes;
// For ConsensusCaller2 -> Temporary delete when merged
extern uint32_t maxNumNodes;
extern int16_t** alMatrix;
extern int *alMatrix_Best;
extern int alignedReads;
extern int unalignedReads;
extern struct Letter_T** alMatrix_Letter;
extern struct timespec ts_start;
extern struct timespec ts_finish;
extern long sumMatrix;
extern long sumTrace;
extern long alignmentTime;

void balancePaths(struct myovlList* G, struct reads* reads);

void contig_repeatFinder();

void hierholzerTourAll(struct myovlList* G, struct reads* reads);

//static inline int* scaffold_findcirclePath(struct scaffEdge* scaffedge, struct jPath* jPath, int jpos, char initRight, struct j_anno* j_anno,int nextpathNum, int uniPathMaxLen, struct contigScaff* uniPath, struct scaffold_set* aS, struct reads* reads,int* pos_len);

void connectPathStats();

struct scaffold_set* scaffold_init2();

struct scaffold_set* scaffold_init3(struct reads* reads);

struct scaffold_set* scaffold_init4(struct scaffold_set* aS);

struct scaffold_set* scaffold_init6(struct scaffold_set* aS, struct reads* reads, char bridging);

void scaffold_printfreqs(struct reads* reads, struct myovlList* G);

struct scaffold_set* scaffold_stats(struct scaffold_set* aS);

void buildBackBone(struct myovlList* ovlgraph, struct string_graph* S, struct reads* reads);

void buildBackBone2(struct myovlList* G, struct reads* reads);

void buildBackBone3(struct myovlList* G, struct reads* reads);

struct contigList* realBackbone(struct myovlList* G, struct reads* reads);

struct contigList* realBackbone2(struct myovlList* G, struct reads* reads);

//Sequence_T* make_poa(struct myovlList* G, struct read* reads);

void poa_part_toDot(char* dotFile, struct Sequence* contig);

void poa_catBackbone2(struct Sequence* contig, struct myovlList *G, char* seq, int leftID, int rightID); // Parameter 3 , struct reads* read

void testFunct();

// 1 align Function
//int poa_align_prepro(struct Sequence* contig, int len, int overhang);

void poa_heuristic_align(struct Sequence* contig, struct reads* read, char* seq, char backbone, int insNum);

char poa_heuristic_align2(struct Sequence* contig, struct reads* read, unsigned char* seq, char backbone, char heuristic, int overhang, int backoverhang);

struct POG* make_poa(struct myovlList* G, struct reads* reads);

struct POG* make_poaScaff(struct myovlList* G, struct reads* reads, char scaffolding, struct para* para, char heuristic);

void poa_toDot(char* dotFile);

void poa_reportVariant(struct POG* pog, char* vcfFile, char* ref);

void poa_consensus(struct Sequence* contig);

void poa_printContigs(struct POG* pog, char* contigFile);

void poa_deleteVariant(struct POG* contigs_pog);

void free_POG(struct POG* contigs_pog);

#endif
