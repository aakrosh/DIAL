/* Routines for reading and writing maf files */

#ifndef MAF_H
#define MAF_H

#include "common.h"
#include "utilities.h"

typedef struct mafblock_st
{
	char* fptr;		/* buffer for reading */
	size_t n;		/* size of the buffer */

	int score;		/* score of the block */
	
	uchar* name1;	/* name of the first sequence */
	uint s1;		/* 0-based start of the alignment */
	uint o1;		/* length of overlap from the first sequence */
	char strand1;	/* which strand of the first sequence*/
	uint len1;		/* total length of the first sequence */
	uchar* aln1;	/* alignment letters from the first sequence */
	uint amax;		/* memory allocated for the first alignment */
	uint alen;		/* length of the alignment*/


	uchar* name2;	/* name of the second sequence */
	uint s2;		/* 0-based start of the alignment from the second seq*/
	uint o2;		/* length of overlap from the second sequence */
	char strand2;	/* which strand of the second sequence*/
	uint len2;		/* length of the second sequence */
	uchar* aln2;	/* alignment from the second sequence */
	uint bmax; 		/* memory allocated for the second alignment */
	uint blen;		/* length of the second alignment */

	FILE* fp;		/* pointer to the file we are reading */
}mafblock;

/* open the maf file corresponding to this name */
mafblock* open_maf_file(const char* const file);

/* open the maf file and return the first block */
mafblock* read_maf_file(const char* const file);

/* get the next maf block from the file */
mafblock* get_next_block(mafblock* const mp);

/*clone a maf block*/
mafblock* copy_mafblock(const mafblock* const mp);

/*just print out this mafblock*/
void print_mafblock(const mafblock* const mp);

/* free all the resources used by this maf block */
void close_maf_file(mafblock* mp);
#endif
