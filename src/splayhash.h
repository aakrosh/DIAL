#ifndef SPLAYHASH_H
#define SPLAYHASH_H

#include "utilities.h"
#include "common.h"

typedef struct splaybin_st
{
	kmer name;				/*the kmer*/
	int counts;				/*number of instances of the kmer*/
	int16_t score;			/*socre of the kmer*/
	struct splaybin_st* left;
	struct splaybin_st* right;
}splaybin;

typedef struct splayhash_st
{
	int po2size;	/*size of the hash table in power of 2*/
	int size;		/*size of the hash table*/
	int mask;		/*the mask to be applied*/

	int elcount;	/*number of elements in the hash table*/
	splaybin** bins;
}splayhash;

/*allocate a new splayhash*/
splayhash* new_splayhash(const int po2size);

/*look up the hash table and return the bin if it does exist*/
splaybin* lookup_splayhash(splayhash* const hash, 
				           const kmer name);

/*add the following and return the bin corresponding to it*/
splaybin* add_splayhash(splayhash* const hash,  	/*the hashtable*/
					    const kmer name,			/*the string*/
				   		const int counts,	   	 	/*the index (0-based)*/
						const int score);		/*the position(0-based)*/

/* free all the resources held by the splayhash */
void free_splayhash(splayhash** phash);
#endif
