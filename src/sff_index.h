//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sff_index.h
//
//----------

#ifndef sff_index_H				// (prevent multiple inclusion)
#define sff_index_H

// other files

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff

// establish ownership of global variables

#ifdef sff_index_owner
#define global
#else
#define global extern
#endif

// sized data types

#include <stdint.h>
typedef uint8_t  u8;
typedef uint32_t u32;
typedef uint64_t u64;

//----------
//
// data structures and types
//
//----------

// header

typedef struct sffindex
	{
	FILE*	f;					// file pointer for the open index file
	int		swapBytes;			// true => swap byte order on reads
	u32		numFiles;			// the number of files in the file index
	u32		numBuckets;			// the number of hash buckets (not including
								// .. the sentinel bucket)
	u32		numSequences;		// the number of sequences in the files
	u32		nameBytes;			// number of bytes for each sequence's name
	u32		filePosBytes;		// number of bytes for each sequence's file
								// .. position
	u32		fileTableOffset;	// position of the beginning of the file table
	u32		hashTableOffset;	// position of the beginning of the hash table
	u32		seqTableOffset;		// position of the beginning of the sequence
								// .. index table
	} sffindex;

// filename table

typedef struct sffft
	{
	u32		num;				// the number of names in the table
	char*	name[1];			// (variable length array);  name[i] is a
								// .. pointer to the ith name (pointing to
								// .. memory within this same allocated block;
								// .. indexed by 0..num-1
	} sffft;

// hash bucket table

typedef uint32_t f31;			// 31-bit equivalent of fpos_t
#define f31_empty 0x80000000	// 'empty' marker for f31 type (see bucket[]
								// .. below)

typedef struct sffbt
	{
	u32		num;				// the number of names in the table
	f31		bucket[1];			// (variable length array);  bucket[h] is the
								// .. position, within the index file, of the
								// .. first sequence with hash value h; if the
								// .. msbit (f31_empty) is set, this bucket is
								// .. empty and the other 31 bits indicate the
								// .. first sequence in the next non-empty
								// .. bucket;  indexed by 0..num (not limited
								// .. to num-1)
	} sffbt;

// seqeunce info table

typedef struct sffst
	{
	u32		allocSize;			// number of bytes allocated
	u32		size;				// the number of sequence info records allocated
	u32		num;				// the number of sequence info records that are
								// .. actually being used
	u32		recordSize;			// the number of bytes in each record of info[]
	char	info[1];			// (variable length array);  the ith record is
								// in info[i*R..(i+1)*R-1], where R is the
								// recordSize and i is the range 0..num-1 (or
								// 0..size-1)
	} sffst;

// layout of data in the file
// $$$ need more here, plus a diagram of the layout

#define sffIndexMagicBig    0x7095D252
#define sffIndexMagicLittle 0x52D29570

#define sffIndexVersion     0x00000100

//----------
//
// prototypes for routines in sff_index.c
//
//----------

sffft*    sff_index_filename_table      (u32 num, char* filenames[]);
sffbt*    sff_index_bucket_table        (u32 num);
sffst*    sff_index_sequence_table      (u32 size, u32 recordSize);
u32       sff_index_hash                (const char* s);
sffst*    sff_index_add_sequence        (sffst* st, u32 h, char* name,
                                         u8 fileNum, f31 filePos);
u32       sff_index_sequence_record_h   (sffst* st, u32 seqIx);
void      sff_index_sort_sequences      (sffst* st);
f31       sff_index_write_header        (FILE* f, sffft* ft, sffbt* bt, sffst* st);
f31       sff_index_write_filenames     (FILE* f, sffft* ft);
f31       sff_index_write_buckets       (FILE* f, sffbt* bt);
f31       sff_index_write_sequence      (FILE* f, sffst* st, u32 seqIx);
f31       sff_index_write_sequences_pad (FILE* f, sffst* st);
sffindex* sff_index_read_header         (FILE* f);
sffft*    sff_index_read_filenames      (sffindex* header);
int       sff_index_lookup_sequence     (sffindex* header, char* seqName,
                                         u32* h, u8* fileNum, u32* filePos);

#undef global
#endif // sff_index_H
