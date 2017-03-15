/*
 * CC.h
 *
 *  Created on: Mar 9, 2017
 *      Author: lkaempfpp
 */

#ifndef CC_H_
#define CC_H_

#include "ConsensusCaller.h"

struct POGreads{
	int ID;				// ID of the read in struct reads (for sequence and length of the read)
	char ori;			// Orientation in regard to the backbone sequence (defined by the first junction read)
	uint32_t start;		// Index of the Backbone node, where the read alignment starts
	uint32_t end;		// Index of the Backbone node, where the read alignment ends
};

struct POGreadsSet{
	uint32_t number;				// Number of reads constructing the contig
	uint32_t size;					// Maximum Container size of reads - threshold for resize
	struct POGreads* pogreads;		// Read Container
};

struct POGseq{
	int length;						/** Backbone/Ref length */
	struct LetterEdge startLetter;	/** First Letter of the contig graph1*/
	char* sequence;					/** Consensus Sequence */
	char* name;						/** Name of Contig */
	float avgCov;					/** Average Coverage of Contig*/
	int nsource_seq;				/** Number of sequences in representing this PO graph */
	struct Variation* var;			/** Linked List of all Variations of the Contig/Scaffold */
	struct Variation* lastvar;		/** Last Variation of the Contig/Scaffold */
	struct sequenceEdge* seqEdge;	/** Edge to the next part of the scaffold connected by a bridge*/
	char vflag;						/** 0 if free to use, 1 go ahead */
	//	uint32_t readleft;				/** */
	//	uint32_t readright;				/** */
	//	char* title;					/** */
	//	struct contigPart* first;		/** First contig of a linked list with all contigs in this scaffold*/
};

struct POGset{
	uint32_t contigNum;
	uint32_t maxNum;
	struct POGseq* contig;
};

struct POG* OLC(struct myovlList* G, struct reads* reads, char scaffolding, char heuristic, struct para* para);

char POG_align(struct reads* reads, struct POGreadsSet* pogreadsSet, char heuristic, uint32_t contigLen);

void POG_alignConsensus(struct POGseq* contig);

void POG_alignConsensus(struct POGseq* contig);

void POG_showMatrix(int row, int column, char* seq);

//void POG_variantCalling(struct POGseq* contig);

//void POG_avgCov(struct POGseq* contig);

#endif /* CC_H_ */
