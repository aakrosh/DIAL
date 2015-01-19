#include "common.h"

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << UINT32_C(8))\
                      +((const uint8_t *)(d))[0])
#endif

uint32_t superfasthash (const char* data, int len) {
uint32_t hash = len, tmp;
int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

/*return the first kmer from this sequence*/
kmer build_index(const uchar* const s, const int klen)
{
	assert(klen < 32);
	kmer word = 0;
	int i;
	for(i = 0; i < klen; i++){
		word = (word << 2);
		word += fasta_encoding[s[i]];
	}

	return word;
}

/*update and return the next kmer in the sequence*/
kmer get_next_kmer(const kmer word, 
				   const uchar* const s, 
				   const int klen, 
				   const int i)
{
	kmer next = (word << (64 - 2*klen + 2));
	next = (next >> (64 - 2*klen));
	next += fasta_encoding[s[i+klen-1]];
	return next;
}

/*return the reverse complement of the given kmer of length klen*/
kmer reverse_complement_kmer(const kmer word, const int klen)
{
	kmer copy = word;
	kmer rc = 0;

	int i = 0;
	for(i = 0; i < klen; i++){
		rc = (rc << 2);
		rc += 3 - (copy & 3);
		copy = (copy >> 2);
	}

	return rc;
}

void convertNtoA(uchar* const sequence, const int length)
{
	int i;
	for(i = 0; i < length; i++){
		sequence[i] = sequence[i] == 'N' ? 'A' : sequence[i];
	}
}

