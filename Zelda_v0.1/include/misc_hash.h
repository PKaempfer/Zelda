/*
 * misc_hash.h
 *
 *  Created on: Aug 11, 2016
 *      Author: kaempfpp
 */

#ifndef MISC_HASH_H_
#define MISC_HASH_H_

#define IS_REPRESENTIBLE_IN_D_BITS(D, N)                \
  (((unsigned long) N >= (1UL << (D - 1)) && (unsigned long) N < (1UL << D)) ? D : -1)

#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (31                                     \
                 + IS_REPRESENTIBLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(16, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(17, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(18, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(19, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(20, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(21, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(22, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(23, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(24, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(25, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(26, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(27, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(28, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(29, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(30, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(31, N)    \
                 + IS_REPRESENTIBLE_IN_D_BITS(32, N)    \
                 )                                      \
   )

struct bucket{
	uint32_t key;
	uint32_t value;
	struct bucket* next;
};

struct hashTable32{
	uint32_t elements;		// Number of elements in HashTable
	uint32_t size;			// Current number of bucktes
	uint16_t maxrp;			// max bucket size before resize (number of reprobes)
	uint32_t bitmask;		// Mask for fit the hash value to the bucket size, number of used bits
	struct bucket** bucket;	// First element of the bucket;
};

struct hashTable32* hashTable32_create(int size);

void hashTable32_stats(struct hashTable32* ht);

void hashTable32_clean(struct hashTable32* ht);

void hashTale32_insert(struct hashTable32* ht, uint32_t key, uint32_t value);

uint32_t hashTable32_find(struct hashTable32* ht, uint32_t key);

uint32_t misc_myhash(uint32_t value);

#endif /* MISC_HASH_H_ */
