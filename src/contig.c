#include "contig.h"

/*convert the character to  a base*/
valid_base char_2_base(const char a)
{
	valid_base base = N;
	switch(a){
		case 'a': base = A; break;
		case 'A': base = A; break;
		case 'c': base = C; break;
		case 'C': base = C; break;
		case 'g': base = G; break;
		case 'G': base = G; break;
		case 't': base = T; break;
		case 'T': base = T; break;
		case '*': base = GAP; break;
		case 'N': base = N; break;
		case 'n': base = N; break;
		default: fatalf("unknown base:%c", a);
	}
	return base;
}

/*convert a base to a character*/
char base_2_char(const valid_base a)
{
	char b;
	switch(a){
		case A: b = 'A'; break;
		case C: b = 'C'; break;
		case G: b = 'G'; break;
		case T: b = 'T'; break;
		case N: b = 'N'; break;	
		default: b = '*'; break;
	}
	return b;
}

/*make a new contig of this size*/
contig* construct_contig(const int num_bases, 
						 const int num_reads, 
						 const int num_segments)
{	
	contig* contig = ckallocz(sizeof(struct contig_st));

	contig->num_bases = num_bases;
	contig->num_reads = num_reads;
	contig->num_segments = num_segments;

	contig->reads = ckallocz(num_reads*sizeof(sff_read**));
	contig->info  = ckallocz(num_bases*sizeof(information));
	return contig;
}

/*create a new allele*/
allele* new_allele(const char base, const char qual, sff_read* const read)
{
	allele* allele = ckalloc(sizeof(struct allele_st));
	allele->base = char_2_base(base);
	allele->quality = qual;
	allele->read = read;
	allele->next = NULL;

	return allele;
}

static void free_contig(contig* contig)
{
	ckfree(contig->header);
	ckfree(contig->reads);

	unsigned int i;   
	information info;	

	for(i = 0; i < contig->num_bases; i++){
		info = contig->info[i];
		slfreelist(&info.alleles);
	}
	ckfree(contig->info);
}

/*free the resources associated with all the contigs*/
void free_contigs(contig** pcontigs)
{
	contig* iter;
	contig* contigs = *pcontigs;
	for(iter = contigs; iter; iter = iter->next){
		free_contig(iter);
	}	
	slfreelist(pcontigs);
	ckfree(*pcontigs);
	*pcontigs = NULL;
}

/*for debugging purposes : print the contig in detail*/
void print_contig(const contig* const contig)
{
	unsigned int i;
	sff_read* read;
	information info;
	allele* iter;

	printf("%s\n", contig->header);
	printf("Number of bases:%d\n", contig->num_bases);
	printf("Number of reads:%d\n", contig->num_reads);
	printf("Number of segments:%d\n", contig->num_segments);
	
	for(i = 0; i < contig->num_reads; i++){
		read = contig->reads[i];
		printf("%s\t%s\n", read->header->name, read->individual);
	}	

	for(i = 0; i < contig->num_bases; i++){
		info = contig->info[i];
		printf("%c(%d)\t", base_2_char(info.consensus), info.quality);
		for(iter = info.alleles; iter; iter = iter->next){
			printf("%c(%d)\t", base_2_char(iter->base), iter->quality);
		}
		printf("\n");
	}
}
