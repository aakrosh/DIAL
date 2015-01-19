//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sff_index.c
//
//----------
//
// sff_index--
//	Support for indexing genomic sequences across several sff files, using a
//	hash-based techinique.
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <string.h>				// standard C string stuff

#include "utilities.h"

#define  sff_index_owner		// (make this the owner of its globals)
#include "sff_index.h"			// interface to this module

//----------
//
// prototypes for private functions
//
//----------

static int qSortByHashAndName (const void* recA, const void* recB);
static u32 hassock_hash       (const void* key, u32 len);

//----------
//
// sff_index_filename_table--
//	Allocate a filename table containing the specified names.
//
//----------
//
// Arguments:
//	u32		num:			The length of filenames[].
//	char*	filenames[]:	The file names.
//
// Returns:
//	A pointer to the newly allocated table.  NULL indicates allocation failed.
//	The caller is repsonsible for de-allocating this table, with a call to
//	free().
//
//----------

sffft* sff_index_filename_table
   (u32		num,
	char*	filenames[])
	{
	sffft*	ft;
	u32		bytesNeeded;
	char*	scan;
	u32		fIx;

	if (num < 1)
		return NULL;

	// allocate, including enough room for the pointer array and names
	// $$$ should check for overflow

	bytesNeeded = sizeof(sffft);
	for (fIx=0 ; fIx<num ; fIx++)
		bytesNeeded += sizeof(char*) + strlen (filenames[fIx]) + 1;

	ft = (sffft*) malloc (bytesNeeded);
	if (ft == NULL) return NULL;

	// initialize 'header' fields

	ft->num = num;

	// set up the pointer array

	scan = (char*) (&ft->name[num]);
	for (fIx=0 ; fIx<num ; fIx++)
		{
		ft->name[fIx] = scan;
		strcpy (/*to*/ scan, /*from*/ filenames[fIx]);
		scan += strlen (filenames[fIx]) + 1;
		}



	return ft;
	}

//----------
//
// sff_index_bucket_table--
//	Allocate an empty hash bucket table.
//
//----------
//
// Arguments:
//	u32		num:	The number of buckets, *not* including the extra bucket at
//					.. the end.  Zero is a legitimate (though unlikely) value.
//
// Returns:
//	A pointer to the newly allocated table.  NULL indicates allocation failed.
//	The caller is repsonsible for de-allocating this table, with a call to
//	free().
//
//----------
//
// Notes:
//	(1)	The contents of the bucket array are not intialized in any way.
//
//----------

sffbt* sff_index_bucket_table
   (u32		num)
	{
	sffbt*	bt;
	u32		bytesNeeded;

	// allocate, including enough room for the array
	// $$$ should check for overflow

	bytesNeeded = sizeof(sffbt) + (num+1) * sizeof(bt->bucket[0]);

	bt = (sffbt*) malloc (bytesNeeded);
	if (bt == NULL) return NULL;

	// initialize 'header' fields

	bt->num = num;

	return bt;
	}

//----------
//
// sff_index_sequence_table--
//	Allocate an empty sequence info table.
//
//----------
//
// Arguments:
//	u32		size:		The number of entries to pre-allocate for.
//	u32		recordSize:	The number of byes in each record.
//
// Returns:
//	A pointer to the newly allocated table.  NULL indicates allocation failed.
//	The caller is repsonsible for de-allocating this table, with a call to
//	free().
//
//----------

sffst* sff_index_sequence_table
   (u32		size,
	u32		recordSize)
	{
	sffst*	st;
	u32		bytesNeeded;

	// $$$ change these hard-wired sizes
	//     .. currently I have 0..3   hash value
	//     ..                  4..18  file name (with terminating zero)
	//     ..                  19     file number
	//     ..                  20..23 file pos

	if (recordSize != 24)
		return NULL;

	// allocate, including enough room for the info array
	// $$$ should check for overflow

	bytesNeeded = sizeof(sffst) + size*recordSize;

	st = (sffst*) malloc (bytesNeeded);
	if (st == NULL) return NULL;

	// initialize 'header' fields

	st->allocSize  = bytesNeeded;
	st->size       = size;
	st->num        = 0;
	st->recordSize = recordSize;

	return st;
	}

//----------
//
// sff_index_add_sequence--
//	Add a sequence info record to a table.
//
//----------
//
// Arguments:
//	sffst*	st:			The table to add to.
//	u32		h:			The sequence's hash code.
//	char*	name:		The sequence's name.
//	u8		fileNum:	The sequence's file number
//	f31		filePos:	The position of the sequence withing the target file.
//
// Returns:
//	A pointer to the table, which may have been moved in the heap (re-
//	allocated).  NULL indicates re-allocation failed (in which case this
//	routine has de-allocated the incoming table).
//
//----------

#define expansionFactor 1.10   // expand by 10%
#define expansionBump   64000  // ... and by 64K entries

sffst* sff_index_add_sequence
   (sffst*	st,
	u32		h,
	char*	name,
	u8		fileNum,
	f31		filePos)
	{
	u32		size, bytesNeeded;
	sffst*	newSt;
	char*	infoRecord;

	if (strlen (name) + 1 > 15)
		return NULL;

	// if we don't have room for this sequence, re-allocate
	// $$$ should check for overflow

	if (st->num == st->size)
		{
		size = (st->size * expansionFactor) + expansionBump;
		bytesNeeded = sizeof(sffst) + size*st->recordSize;

		newSt = (sffst*) realloc (st, bytesNeeded);
		if (newSt == NULL)
			{ free (st);  return NULL; }

		st = newSt;
		st->allocSize  = bytesNeeded;
		st->size       = size;
		}

	// write the info record into the table
	// $$$ change these hard-wired sizes

	infoRecord = st->info + (st->num * st->recordSize);

	*(u32*) (infoRecord + 0)  = h;
	*(u8*)  (infoRecord + 19) = fileNum;
	*(u32*) (infoRecord + 20) = filePos;
	strncpy (/*to*/ infoRecord+4, /*from*/ name, /*limit*/ 15);

	st->num++;

	return st;
	}

//----------
//
// sff_index_sequence_record_h--
//	Extract the previously-computed hash value of a sequence info record.
//
//----------
//
// Arguments:
//	sffst*	st:		The sequence table.
//	u32		seqIx:	The index of the sequence to write.  This is a number in
//					.. the range 0..st->num-1.
//
// Returns:
//	The hash of that sequence info record's name.
//
//----------

u32 sff_index_sequence_record_h
   (sffst*	st,
	u32		seqIx)
	{
	char*	infoRecord;
	u32		h;

	// $$$ change these hard-wired sizes

	infoRecord = st->info + (seqIx * st->recordSize);

	h = *(u32*) (infoRecord + 0);

	return h;
	}

//----------
//
// sff_index_sort_sequences--
//	Sort a sequence info table.
//
//----------
//
// Arguments:
//	sffst*	st:	The table to sort.
//
// Returns:
//	(nothing)
//
//----------

void sff_index_sort_sequences
   (sffst*	st)
	{
	qsort (st->info, st->num, st->recordSize, qSortByHashAndName);
	}

//----------
// [[-- comparison function for the standard c function qsort --]]
//
// qSortByHashAndName--
//	Compare two sequence infor records, so that qsort will sort by increasing
//	hash value, with ties broken by alphabetically increasing name.
//
//----------
//
// Arguments:
//	const void* recA:	Pointer to one record.
//	const void* recB:	Pointer to record.
//
// Returns:
//	> 0 => recA is greater than recB.
//	= 0 => recA and recB are the same.
//	< 0 => recA is less than recB.
//
//----------

// $$$ change these hard-wired sizes, passing them in as globals
// $$$ need to test that this routine creates the right order

static int qSortByHashAndName
   (const void*	_recA,
	const void*	_recB)
	{
	char*	recA = (char*) _recA;
	char*	recB = (char*) _recB;
	u32		hashA, hashB;
	char*	nameA, *nameB;

	hashA = *(u32*) (recA + 0);
	hashB = *(u32*) (recB + 0);

	if      (hashA < hashB) return -1;
	else if (hashA > hashB) return  1;

	nameA = recA + 4;
	nameB = recB + 4;

	return strcmp (nameA, nameB);
	}

//----------
//
// hassock_hash--
//	Compute a variant of Austin Appleby's MurmurHash2.
//
//----------
//
// Arguments:
//	const void*	key:	The data block to hash.
//	u32			len:	The length of that block.
//
// Returns:
//	A hash of the block.
//
//----------
//
// Notes:
//
//	(1)	As of Apr/2009, information about this hash fucntion can be found at
//		  murmurhash.googlepages.com
//	(2) This implementation is based on an implementation found at
//	      murmurhash.googlepages.com/MurmurHashNeutral2.cpp
//	    It differs in the following ways:
//	      (a) The "seed" is hardwired.
//		  (b) We parse the data block in reverse;  this allows the caller to
//            prepend an additional seed pattern to his buffer, potentially
//	          getting better mixing for the bits in the final incorporated
//	          bytes.
//		  (c) The last three bytes are incorporated in a different order than
//	          they were in MurmurHash2, because the code just works out better
//	          this way.
//
//----------

u32 sff_index_hash (const char* s)
	{ return hassock_hash (s, strlen(s)); }


static u32 hassock_hash
   (const void*	key,
	u32			len)
	{
	const u32	seed = 0x5C3FC4D3;
	const u32	m    = 0x87C10417;
	const int	r    = 24;
	const u8*	data = ((const u8*) key) + len;
	const u8*	stop = ((const u8*) key) + 4;
	u32			h, k;

	h = seed ^ len;
	while (data >= stop)
		{
		k  = *(--data);
		k |= *(--data) << 8;
		k |= *(--data) << 16;
		k |= *(--data) << 24;

		k *= m;
		k ^= k >> r;
		k *= m;

		h *= m;
		h ^= k;

		len  -= 4;
		}

	switch (len)
		{
		case 3: h ^= *(--data) << 16;
		case 2: h ^= *(--data) << 8;
		case 1: h ^= *(--data);
		        h *= m;
		};

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	//printf ("%08X %s\n", h, (char*) key);

	return h;
	}

//----------
//
// sff_index_write_header--
//	Write the header for an index file.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to write to.
//	sffft*	ft:	The filename table.
//	sffbt*	bt:	The bucket table.
//	sffst*	st:	The sequence table.
//
// Returns:
//	The number of bytes written to the file.  0 indicates a failure.
//
//----------
//
//	offset 0x00: 70 95 D2 52        big endian magic number
//									.. (52 D2 95 70 => little endian)
//	offset 0x04: 00 00 01 00        version 1.0
//	offset 0x08: 00 00 00 24        header length (in bytes, including this
//									.. field, but not the magic number or
//									.. version)
//	offset 0x0C: 00 00 00 xx        NB, bytes for each sequence's name
//									.. (including terminating zero)
//	offset 0x10: 00 00 00 xx        PB, bytes for each sequence's file position
//	offset 0x14: 00 00 00 xx        FN, number of files
//	offset 0x18: xx xx xx xx        FO, offset to file table
//	offset 0x1C: xx xx xx xx        HN, number of hash buckets (without extra
//									.. sentinel bucket)
//	offset 0x20: xx xx xx xx        HO, offset to hash table
//	offset 0x24: xx xx xx xx        SN, number of sequences
//	offset 0x28: xx xx xx xx        SO, offset to sequence index table
//
//----------

static u32 pad_for_32 (u32 n);
static u32 pad_for_32 (u32 n) { return (32 - (n % 32)) % 32; }

static char zeros[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };


f31 sff_index_write_header
   (FILE*	f,
	sffft*	ft,
	sffbt*	bt,
	sffst*	st)
	{
	u32		seqRecordSize = 20;			// $$$ de-hardwire these
	u32		nameBytes     = 15;
	u32		filePosBytes  = 4;
	u32		magic, version;
	u32		headerLength, headerPad, headerSize;
	u32		fileTableOffset, fileTableLength, fileTablePad, fileTableSize;
	u32		seqTableOffset, seqTableLength, seqTablePad, seqTableSize;
	u32		hashTableOffset;
	u32		fIx;

	// determine where things will layout in the file

	headerLength    = 0x24;
	headerPad       = pad_for_32(0x08 + headerLength);
	headerSize      = 0x08 + headerLength + headerPad;

	fileTableOffset = headerSize;
	fileTableLength = 0;
	for (fIx=0 ; fIx<ft->num ; fIx++)
		fileTableLength += strlen (ft->name[fIx]) + 1;
	fileTablePad    = pad_for_32(fileTableLength);
	fileTableSize   = fileTableLength + fileTablePad;

	seqTableOffset  = fileTableOffset + fileTableSize;
	seqTableLength  = st->num * seqRecordSize;
	seqTablePad     = pad_for_32(seqTableLength);
	seqTableSize    = seqTableLength + seqTablePad;

	hashTableOffset = seqTableOffset + seqTableSize;

	// write the header
	// $$$ should check for errors

	magic   = sffIndexMagicBig;
	version = sffIndexVersion;

	fwrite (&magic,           4, 1, f);
	fwrite (&version,         4, 1, f);
	fwrite (&headerLength,    4, 1, f);
	fwrite (&nameBytes,       4, 1, f);
	fwrite (&filePosBytes,    4, 1, f);
	fwrite (&ft->num,         4, 1, f);
	fwrite (&fileTableOffset, 4, 1, f);
	fwrite (&bt->num,         4, 1, f);
	fwrite (&hashTableOffset, 4, 1, f);
	fwrite (&st->num,         4, 1, f);
	fwrite (&seqTableOffset,  4, 1, f);

	if (headerPad > 0) fwrite (zeros, headerPad, 1, f);

	return headerSize;
	}

//----------
//
// sff_index_write_filenames--
//	Write the filename table for an index file.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to write to.
//	sffft*	ft:	The filename table.
//
// Returns:
//	The number of bytes written to the file.  0 indicates a failure.
//
//----------

f31 sff_index_write_filenames
   (FILE*	f,
	sffft*	ft)
	{
	u32		nameLength, fileTableLength, fileTablePad;
	u32		fIx;

	fileTableLength = 0;
	for (fIx=0 ; fIx<ft->num ; fIx++)
		{
		nameLength = strlen (ft->name[fIx]);
		fwrite (ft->name[fIx], nameLength + 1, 1, f);
		fileTableLength += nameLength + 1;
		}
	fileTablePad = pad_for_32(fileTableLength);

	if (fileTablePad > 0) fwrite (zeros, fileTablePad, 1, f);

	return fileTableLength + fileTablePad;
	}

//----------
//
// sff_index_write_buckets--
//	Write the hash bucket table for an index file.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to write to.
//	sffbt*	bt:	The bucket table to write.
//
// Returns:
//	The number of bytes written to the file.  0 indicates a failure.
//
//----------

f31 sff_index_write_buckets
   (FILE*	f,
	sffbt*	bt)
	{
	u32		hashTableLength;

	hashTableLength  = (bt->num+1) * sizeof(bt->bucket[0]);
	fwrite (bt->bucket, hashTableLength, 1, f);

	return hashTableLength;
	}

//----------
//
// sff_index_write_sequence--
//	Write one sequence info record for an index file.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write to.
//	sffst*	st:		The sequence table.
//	u32		seqIx:	The index of the sequence to write.  This is a number in
//					.. the range 0..st->num-1.
//
// Returns:
//	The number of bytes written to the file.  0 indicates a failure.
//
//----------

f31 sff_index_write_sequence
   (FILE*	f,
	sffst*	st,
	u32		seqIx)
	{
	char*	infoRecord;
	int		fileNum;
	u32		filePos;
	char*	name;

	// $$$ change these hard-wired sizes

	infoRecord = st->info + (seqIx * st->recordSize);

	fileNum = *(u8*)  (infoRecord + 19);
	filePos = *(u32*) (infoRecord + 20);
	name    =         (infoRecord + 4);

	//fprintf (stderr, "%s %3d %08X\n", name, fileNum, filePos);

	fwrite (name,     15, 1, f);
	fwrite (&fileNum, 1,  1, f);
	fwrite (&filePos, 4,  1, f);

	return 20;	// $$$ change this hard-wired size
	}

//----------
//
// sff_index_write_sequences_pad--
//	Write padding for the sequence info table for an index file.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to write to.
//	sffst*	st:		The sequence table.
//
// Returns:
//	The number of bytes written to the file.  0 indicates a failure.
// $$$ no, 0 might be a legitimate padding value
//
//----------

f31 sff_index_write_sequences_pad
   (FILE*	f,
	sffst*	st)
	{
	u32		seqRecordSize = 20;			// $$$ de-hardwire this
	u32		seqTableLength, seqTablePad;

	seqTableLength = st->num * seqRecordSize;
	seqTablePad    = pad_for_32(seqTableLength);

	if (seqTablePad > 0) fwrite (zeros, seqTablePad, 1, f);

	return seqTablePad;
	}

//----------
//
// sff_index_read_header--
//	Read the header from an index file.
//
//----------
//
// Arguments:
//	FILE*	f:	The file to read from.
//
// Returns:
//	A pointer to the newly allocated header.  NULL indicates read or allocation
//	failed.  The caller is repsonsible for de-allocating this table, with a call
//	to free().
//
//----------

static void swap_32 (u32* _v);
static void swap_32 (u32* _v)
	{
	u32 v = *_v;
	*_v = (( v        & 0x000000FF) << 24)
	    + (((v >>  8) & 0x000000FF) << 16)
	    + (((v >> 16) & 0x000000FF) <<  8)
	    +  ((v >> 24) & 0x000000FF);
	}


sffindex* sff_index_read_header
   (FILE*	f)
	{
	int			swapBytes;
	u32			magic, version, headerLength;
	sffindex*	header;

	// read and validate the magic number, version, and header size

	if(fread (&magic, 4, 1, f) != 1){   
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
	if (magic == sffIndexMagicBig)
		swapBytes = false;
	else if (magic == sffIndexMagicLittle)
		swapBytes = true;
	else
		return NULL;

	if(fread (&version, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&version);
	if (version != sffIndexVersion)
		return NULL;

	if(fread (&headerLength, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }  
    if (swapBytes) swap_32(&headerLength);
	if (headerLength != 0x24)
		return NULL;

	// allocate the header

	header = (sffindex*) malloc (sizeof(sffindex));
	if (header == NULL) return NULL;

	header->f         = f;
	header->swapBytes = swapBytes;

	if(fread (&header->nameBytes, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }  
    if (swapBytes) swap_32(&header->nameBytes);
	if(fread (&header->filePosBytes, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }

    }  
    if (swapBytes) swap_32(&header->filePosBytes);
	if(fread (&header->numFiles, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }  
    if (swapBytes) swap_32(&header->numFiles);
	if(fread (&header->fileTableOffset, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&header->fileTableOffset);
	if(fread (&header->numBuckets, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&header->numBuckets);
	if(fread (&header->hashTableOffset, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&header->hashTableOffset);
	if(fread (&header->numSequences, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&header->numSequences);
	if(fread (&header->seqTableOffset, 4, 1, f) != 1){
        if(feof(f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }
    if (swapBytes) swap_32(&header->seqTableOffset);

	return header;
	}

//----------
//
// sff_index_read_filenames--
//	Read the filename table from an index file.
//
//----------
//
// Arguments:
//	sffindex*	header:	The index to read from.
//
// Returns:
//	A pointer to the newly allocated table.  NULL indicates read or allocation
//	failed.  The caller is repsonsible for de-allocating this table, with a call
//	to free().
//
//----------

sffft* sff_index_read_filenames
   (sffindex*	header)
	{
	u32			bytesNeeded;
	sffft*		ft;
	u32			fIx;
	char*		scan;
	char		ch;

	// read the file name table once, just to get the size

	fseek (header->f, header->fileTableOffset, SEEK_SET);

	bytesNeeded =  sizeof(sffft) + (header->numFiles * sizeof(char*));
	bytesNeeded += header->numFiles; // (one byte for each terminator)
	for (fIx=0 ; fIx<header->numFiles ; fIx++)
		{ while (fgetc (header->f) != 0) bytesNeeded++; }

	// allocate

	ft = (sffft*) malloc (bytesNeeded);
	if (ft == NULL) return NULL;

	// initialize 'header' fields

	ft->num = header->numFiles;

	// re-read the file name table, copying the names and setting up the
	// pointer array

	fseek (header->f, header->fileTableOffset, SEEK_SET);

	scan = (char*) (&ft->name[header->numFiles]);
	for (fIx=0 ; fIx<header->numFiles ; fIx++)
		{
		ft->name[fIx] = scan;
		while ((ch = fgetc (header->f)) != 0) *(scan++) = ch;
		*(scan++) = 0;
		}

	return ft;
	}

//----------
//
// sff_index_lookup_sequence--
//	Lookup one sequence info record from an index file.
//
//----------
//
// Arguments:
//	sffindex*	header:		The index to read from.
//	char*		seqName:	The sequence name to locate.
//	u32*		h:			Place to return the sequence name's hash value.
//							.. This can be NULL if the caller is not interested.
//	u8*			fileNum:	Place to return the number of the file that contains
//							.. the sequence.  This is an index into the index's
//							.. filename table.  If the sequence is not found
//							.. this value is not changed.
//	u32*		filePos:	Place to return the offset within that file of the
//							.. beginning of the sequence.  This is suitable for
//							.. fseek().  If the sequence is not found this value
//							.. is not changed.
//
// Returns:
//	true if the sequence was found;  false if not.
//
//----------

int sff_index_lookup_sequence
   (sffindex*	header,
	char*		seqName,
	u32*		_h,
	u8*			fileNum,
	u32*		_filePos)
	{
	u32			h;
	f31			thisBucket, nextBucket;
	char		seqRecord[20];				// $$$ de-hardwire this
	u32			filePos;

	if (header->nameBytes + 1 + header->filePosBytes != sizeof(seqRecord))
		return false;

	// find the bucket that should contain this sequence

	h = sff_index_hash (seqName) % header->numBuckets;
	if (_h != NULL) *_h = h;

	fseek (header->f, header->hashTableOffset + (4*h), SEEK_SET);
	if(fread (&thisBucket, 4, 1, header->f) != 1){
         if(feof(header->f) != 0){
            fatal("error in reading bytes from the stream, eof reached");
        }else if(ferror(header->f) != 0){
            fatal("error in reading bytes from the stream, error");
        }
       
    }  
    if (header->swapBytes) swap_32(&thisBucket);

	// if this bucket is empty, we don't have that sequence

	if ((thisBucket & f31_empty) != 0)
		return false;

	// get the end of the bucket (the start of the next bucket)

	if(fread (&nextBucket, 4, 1, header->f) != 1){
        if(feof(header->f) != 0){
            fatal("error in reading 4 bytes from the stream, eof reached");
        }else if(ferror(header->f) != 0){
            fatal("error in reading 4 bytes from the stream, error");
        }
    }  
    if (header->swapBytes) swap_32(&nextBucket);
	nextBucket &= ~f31_empty;

	// search the bucket for this sequence name
	// $$$ this could be improved by
	// $$$     ... doign a binary search (instead of this linear scan)
	// $$$     ... reading more than one sequence record at a time

	fseek (header->f, thisBucket, SEEK_SET);

	while (thisBucket < nextBucket)
		{
		if(fread (seqRecord, sizeof(seqRecord), 1, header->f) != 1){
            if(feof(header->f) != 0){
                fatal("error in reading bytes from the stream, eof reached");
            }else if(ferror(header->f) != 0){
                fatal("error in reading bytes from the stream, error");
            }
        }
		thisBucket += sizeof(seqRecord);
		if (strcmp (seqRecord, seqName) < 0)
			continue;

		if (strcmp (seqRecord, seqName) > 0)
			return false;

		*fileNum = *(u8*)  (seqRecord + 15);		// $$$ de-hardwire these
		filePos  = *(u32*) (seqRecord + 16);
		if (header->swapBytes) swap_32(&filePos);
		*_filePos = filePos;
		return true;
		}

	return false;
	}

