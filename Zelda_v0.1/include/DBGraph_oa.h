/*
 ============================================================================
 Name        : DBGraph_oa.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Converts the dBG from (random access) hash table to cache
 	 	 	   coherent adjacency list and calculates edges
 ============================================================================
 */

#ifndef DBGRAPH_OA_H_
#define DBGRAPH_OA_H_

uint32_t getKmer_oa(KmerBitBuffer kmer);

void hashToTabDFS_oa();

void travDFS_oa(uint32_t);

void goUpDFS_2_oa(uint32_t**,uint32_t**, uint32_t*, uint32_t*, uint32_t*, uint32_t*, int);

#endif /* DBGRAPH_OA_H_ */
