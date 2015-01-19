#ifndef COMMON_H
#define COMMON_H

#include "utilities.h"

#define _ -1

#define MAX_READ_SIZE 4096

typedef uint64_t kmer;

static const signed char fasta_encoding[] = {
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,_,_,_,_,_,
_,_,_,_,_,0,_,1,_,_,
_,2,_,_,_,_,_,_,0,_,
_,_,_,_,3,_,_,_,_,_,
_,_,_,_,_,_,_,0,_,1,
_,_,_,2,_,_,_,_,_,_,
0,_,_,_,_,_,3,_,_,_,
_,_,_,_,_,_,_,_,_
};

static const char bit_encoding[] = {'A', 'C', 'G', 'T'};

/*the hash function from http://www.azillionmonkeys.com/qed/hash.html*/
uint32_t superfasthash (const char* data, int len);

/*return the first kmer from this sequence*/
kmer build_index(const uchar* const s, const int klen);

/*update and return the next kmer in the sequence*/
kmer get_next_kmer(const kmer word, 
				   const uchar* const s, 
				   const int klen, 
				   const int i);

/*return the reverse complement of the given kmer of length klen*/
kmer reverse_complement_kmer(const kmer word, const int klen);

void convertNtoA(uchar* const sequence, const int length);

#endif
