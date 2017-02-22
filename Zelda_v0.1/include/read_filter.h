/*
 * read_filter.h
 *
 *  Created on: Feb 22, 2017
 *      Author: lkaempfpp
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
