/*
 ============================================================================
 Name        : read_filter.h
 Author      : Kämpfer, Philipp
 Version     : v0.1
 Copyright   : GPLv3 (general public license)
 Description : Prefilter to correct erroneous reads based on coverage drops
 	 	 	   in the initial build hash table.
 ============================================================================
 */

#ifndef READ_FILTER_H_
#define READ_FILTER_H_

#include "readDB.h"
#include "FileReader.h"

struct filter_block{
	int pthr_id;
	int pthr_num;
	long start;
	long end;
	struct reads* reads;
};

void filter_reads(struct reads* reads, const int NUM_READS,pthread_t* threads);

#endif /* READ_FILTER_H_ */
