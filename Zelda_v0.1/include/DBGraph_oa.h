/*
 * DBGraph_oa.h
 *
 *  Created on: Feb 10, 2016
 *      Author: kaempfpp
 */

#ifndef DBGRAPH_OA_H_
#define DBGRAPH_OA_H_

uint32_t getKmer_oa(KmerBitBuffer kmer);

void hashToTabDFS_oa();

void travDFS_oa(uint32_t);

void goUpDFS_2_oa(uint32_t**,uint32_t**, uint32_t*, uint32_t*, uint32_t*, uint32_t*, int);

#endif /* DBGRAPH_OA_H_ */
