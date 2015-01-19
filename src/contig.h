#ifndef contig_H
#define contig_H

#include "utilities.h"
#include "sff.h"
#include "slinklist.h"

typedef enum valid_base_en {A,C,G,T,GAP,N} valid_base;

typedef struct allele_st
{
	struct allele_st* next;
	valid_base   base;
	uint8_t quality;
	sff_read* read;
}allele;

typedef struct information_st
{
	valid_base consensus;
	uint8_t quality;
	allele* alleles;
}information;

typedef struct contig_st
{
	struct contig_st* next;	/*the next contig in the assembly*/
	char* header;			/*header for the contig*/
	unsigned int num_bases;	
	unsigned int num_reads;
	unsigned int num_segments;
	information*  info;		/*information about this site*/	
	sff_read** reads;		/*reads associated with this contig*/
}contig;

typedef struct af_segment_st
{
	struct af_segment_st* next;
	bool complement;
	int offset;
}af_seg;

/*convert the character to  a base*/
valid_base char_2_base(const char a);

/*convert a base to a character*/
char base_2_char(const valid_base a);

/*make a new contig of this size*/
contig* construct_contig(const int num_bases, 
						 const int num_reads, 
						 const int num_segments);

/*create a new allele*/
allele* new_allele(const char base, const char qual, sff_read* const read);

/*free the resources associated with all the contigs*/
void free_contigs(contig** pcontigs);

/*for debugging : print the contig in detail*/
void print_contig(const contig* const contig);
#endif
