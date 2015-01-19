/*This module reads the following:
 *		A file with the exons and the names of the reads hitting them
 chr1_1497
 E88BQJZ01A813B cedric
 E88BQJZ02FPN4Q cedric
 E86HODY01CHHBZ spirit
 FF4D8XH01DBURR spirit
 FF4D8XH01C2UIA spirit
 FF4D8XH01EZAQD spirit
 FF4D8XH02QU0GF spirit

 chr1_2983
 E88BQJZ02GXSXA cedric
 E86HODY02GRS9A spirit
 E86HODY02HKGQO spirit
 FF4D8XH01DZQX8 spirit
 FF4D8XH01EEVOG spirit

 * 		A directory name with all the sff files in it
 *
 *and outputs the following:
 
 chr1_1497
 >contig00001  length=450   numreads=4
 TTGGCCTTAGGAGAAGcAGATGATTGtGTTTgCCaTGagaGTgATacaTTTTCcCTGGAT
 TTGTCTTCTAGAGaTTTTtcTTgCAGAtCtATCAGGATGAGCATCCAGGCCCcACCcAGA
 CTCCTGGAGCTGgCGGGGCAGAGcTGCTGAGAGACcAagCccTTGtCCATCTCTGCCaTG
 GAGGAGCTGCCCAGGGtGCTCTAtcTCCcACTCTTCaTGGAGGCCTtCcGCAGGAGACAC
 TTCcAGACTCTGACAGTGATGgTGCAGGCCTGGCCcTtCACCcgCCTCCCtCtGGGATCG
 CTGATGAAGACGCTtCATCTGgaGaCcTTAAAaGCATTGCTGGAAGGGCTTCACATGCTG
 CTTACaCAGAaGGATCGCCCCaGGtgAGGTGACCCAGGAAGGCTGGTAGATGGGGCTCAG
 GTGTCcAGGGAAaGAAcgAGCaGGGTCAgg
 38      39      cedric:G=2(34) spirit:A=2(34)
 45      46      cedric:A=2(60) cedric:C=2(40) 

 chr1_2983
 
 * which is the list of contigs and alleles in the exon. It also list 
 * the singletons and the alignments if the proper options are specified.
 * For the most part this should be the same as assemble.c, but this version 
 * is only handy when I want to speed up a lot of stuff and hence, I only 
 * update this version once in a while
 **/

//#include <sys/dir.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <errno.h>
#include <unistd.h>

#include "utilities.h"
#include "sequences.h"
#include "hashtable.h"
#include "sff.h"
#include "contig.h"
#include "sff_index.h"

#define USE "assemble1 sff_directory exons.txt [-singletons] [-alignments] [-onlyassemble] [--newbler=directory_with_newbler] [--index=index_file]"

#define READS(A) ((A[0]) || (A[1]) || (A[2]) || (A[3]) || (A[4]))

/*number of flows per read for the gs20 sequences*/
#define GS20 400

/*number of flows per read for the gsflx sequences*/
#define GSFLX 800

/*a minimum of 4 reads should hit the position for it to be viable for a
 * SNP/allele candidate*/
#define MIN_HITS 4

/*a minimum of 2 reads should corroborate the base for an allele*/
#define MIN_GAPS 2

/*read these many exons in batch */
#define BATCH 1500000

unsigned int DSIZE = 4096;

/*this is the directory where 454 binaries are kept*/
char* rig_directory="/usr/local/rig/bin";

/* is there an index file to reduce the IO used by this program. This
 * (theoretically) should make this pipeline more ameable to clusters */
char* index_file = NULL;

/*command line arguments*/
bool Singletons     = FALSE;
bool Alignments     = FALSE;
bool OnlyAssemble   = FALSE;

/*the generic buffer used*/
char* buffer;

/*the mapping in a contig from contig positions in the ace file to the one in
 * AllContigs.fna*/
int* mapping;

/*the command to be run using "system"*/
char* command;

/*The accession number of the read which is read from multiple files and used as
 * an unique identifier*/	
char accession[4096];

/*The information accompanying accession identifier for the read*/
char gossip[4096];

/*a linked list of read names*/
typedef struct readlist_st
{
	struct readlist_st* next;
	char* name;
	bool isgap;
}readlist;

typedef struct snp_st
{
	struct snp_st* next;
	int ps[4];
}snp;

/*print the singletons for this exon*/
static void print_exon_singletons(hashtable* const hash, 
								  const char* const dir)
{
	/*open the 454ReadStatus file to find which of the reads are singletons. 
	 *Then use sffinfo -s to dump the reads from 1.sff and 2.sff and print 
	 *the ones which were Singletons */
	sprintf(buffer,"%s/exon/454ReadStatus.txt", dir);
	FILE* fp = ckopen(buffer,"r");
	size_t n = 1;
	char* fptr = ckalloc(n+1);		
	sff_read* read = NULL;
	int count = 1;
	uint32_t i = 0;
	unsigned int j = 0;
	char* ptr;

	while(getline(&fptr, &n, fp) != -1 && 
		  sscanf(fptr,"%s %s\n",accession, gossip) == 2){
		fprintf(stderr,"%s", fptr);
		if(strcmp(gossip,"Singleton") == 0){
			/*print this read*/
			printf(">singleton%05d\t%s\n", count++, accession+1);
		   	if((ptr = strchr(accession, '.')) != NULL){
				*ptr = '\0';
			}
			if((ptr = strchr(accession, '_')) != NULL){
				*ptr = '\0';
			}

			read = must_find_hashtable(hash, accession,strlen(accession));
			for(i = 0, j = 0; i < read->header->num_bases; i++, j++){
				buffer[i] = read->data->bases[i];
			}
			format((uchar*)buffer, j);
		}
	}

	fclose(fp);
	ckfree(fptr);
}

/*print this aligned read*/
static void print_aligned_reads(sff_read** const reads, 
						 		const unsigned int num_reads,
								hashtable* const individuals)
{
	unsigned int i, j, k;
	int binindex;
	ancilliary* aln;	
	bin* iter;

	for(binindex = 0; binindex < individuals->size; binindex++){
		iter = individuals->bins[binindex];
		while(iter){
			/*print all the reads which are from the individual iter->name*/
			printf("%s\n", iter->name);
			for(i = 0; i < num_reads; i++){
				aln = reads[i]->information;
				if(strcmp(reads[i]->individual, iter->name) == 0){
					j = aln->begin + aln->offset;		
					for(k = 0; k < j; k++){
						printf(" ");
					}
					for(k = aln->begin; k <= aln->end; k++){
						printf("%c", aln->sequence[k].seq);
					}
					printf("\n");
				}
			}
			
			iter = iter->next;
		}
	}
}

/*return the count of gaps or nongap bases in this structure*/
static inline int counts(const information info, const bool isgap)
{
	allele* iter;
	int count;

	if(isgap){
		for(iter = info.alleles, count = 0; iter; iter = iter->next){
			if(iter->base == GAP){
				count++;
			}
		}
	}else{
		for(iter = info.alleles, count = 0; iter; iter = iter->next){
			if(iter->base != GAP){
				count++;
			}
		}
	}
	return count;
}

/* return a list of names of reads at this position of the contig. */
static readlist* select_reads(const information info)
{
	readlist* list = NULL;
	readlist* node;

	allele* iter = info.alleles;

	for(;iter; iter = iter->next){
		node = ckalloc(sizeof(readlist));
		node->isgap = iter->base == GAP ? TRUE : FALSE;
		node->name = iter->read->header->name;
		sladdhead(&list, node);
	}
	slreverse(&list);

	return list;
}

/*what is the consensus base in this position in the contig. Newbler has just
 * one base in a position, and that means that as soon as we encounter the first
 * non gap, non N base, we can return*/
static valid_base call(const information info)
{
	allele* iter;
	for(iter = info.alleles; iter; iter = iter->next){
		if(iter->base < GAP){
			return iter->base;
		}
	}
	return GAP;
}

/* return the index of the last base of this stretch of homopolymers which
 * starts at start (0-based)*/
static int move_ptr(const contig* const contig, const int start)
{
	valid_base consensus = call(contig->info[start]);
	uint index = start + 1;
	while(index < contig->num_bases && call(contig->info[index]) == consensus){
		index++;
	}

	return index - 1;
}

/* given a position on  a contig and a list of reads, return TRUE if all the
 * reads in the list exist AND all the reads in the list that have a non-GAP
 * base,  have a gap in that location. */
static bool check_location(const information info, readlist* reads)
{
	/* at least one of the reads in the list gaps, should have a non-GAP base */
	readlist* gap = reads;
	
	allele* iter = info.alleles;
	bool mod = FALSE;
	
	for(iter = info.alleles; iter; iter = iter->next){
		if(gap && 0 == strcmp(iter->read->header->name, gap->name)){
			if(gap->isgap == TRUE && iter->base != GAP){
				mod = TRUE;
				break;
			}
			gap = gap->next;
		}
	}
	if(mod == FALSE){
		return FALSE;
	}

	/* all the reads in the nongap should have a GAP in this location */
	gap = reads;
	for(iter = info.alleles; iter; iter = iter->next){
		if(gap && 0 == strcmp(iter->read->header->name, gap->name)){
			if(gap->isgap == FALSE && iter->base != GAP){
				return FALSE;
			}
			if(gap->isgap == TRUE && iter->base != GAP){
				gap->isgap = FALSE;
			}
			gap = gap->next;
		}
	}

	return TRUE;
}

/* is this a location of a candidate snp? */
static bool check_locations(contig* const contig,
							const int start,
							int* const ploc,
							int* const pmax)
{
	/*bunch of initialisation*/
	memset(ploc, -1, 4*sizeof(int));
	*pmax = -1;
	int arrayptr = 0;

	/*This would be the last base of a a homopolymer stretch. Move to the last
	 * base of the next homopolymer stretch. If it is a stretch of 1, then there
	 * might be more stretches that need to be examined. Also it cannot be a
	 * stretch of a base that has already occured in this examination. */
	char check[5] = {0,0,0,0,0};
	valid_base base = N;

	valid_base consensus = call(contig->info[start]);
	check[consensus] = 1;

	if(consensus == N || consensus == GAP){
		return FALSE;
	}

	ploc[arrayptr++] = start;
	*pmax = *pmax + 1;

	int end = start + 1;
	while(contig->info[end].consensus == GAP &&
	      call(contig->info[end]) == consensus){
		end++;
	}

	/* is this a true candidate SNP? */
	bool result = TRUE;

	/* create a list of reads that have a gap at this location */
	readlist* reads = select_reads(contig->info[start]);

	int index = end - 1;	
	int prev  = end - 1;
	while(1){
		prev = index;

		if((index + 1) >= (int)contig->num_bases){
			break;
		}

		/* move to the next index */
		index = move_ptr(contig, index + 1);
		base  = call(contig->info[index]);
		if((index - prev) > 1 || check[base] == 1){
			break;
		}
		if(contig->info[index].consensus == N || base == N || base == GAP){
			continue;
		}
		check[base] = 1;

		if(counts(contig->info[index], TRUE)     >= MIN_GAPS &&
		   slcount(contig->info[index].alleles)  >= MIN_HITS){
			if(check_location(contig->info[index], reads)){
				if(counts(contig->info[index], FALSE) >= MIN_GAPS){
					ploc[arrayptr++] = index;
					*pmax = *pmax + 1;
				}
			}else{
				result = FALSE;
				break;
			}
		}else{
			result = FALSE;
			break;
		}
	}

	if(result == TRUE && check[base] != 1){
		/*examine the last stretch*/
		if(counts(contig->info[index], TRUE)     >= MIN_GAPS &&
		   slcount(contig->info[index].alleles) >= MIN_HITS){
			if(check_location(contig->info[index], reads)){
				if(counts(contig->info[index], FALSE) >= MIN_GAPS){
					ploc[arrayptr++] = index;
					*pmax = *pmax + 1;
				}
			}else{
				result = FALSE;
			}
		}
	}

	if(arrayptr > 1){
		result = TRUE;
		/*lets update the consensus sequence */
		contig->info[start].consensus = call(contig->info[start]);
		contig->info[start].quality   = 0;
		int i;
		for(i = 1; i <= *pmax; i++){
			contig->info[ploc[i]].consensus = GAP;
			contig->info[ploc[i]].quality   = 0;
		}
	}else{
		result = FALSE;
	}

	slfreelist(&reads);

	return result; 
}

/* print the quality values for this individual with this allele */
static inline void print_quality_values(const contig* const contig,
										const snp* const candidate,
										const char* const name,
										const int max,
										const valid_base base)
{
	int i, loc, count = 0;
	allele* iter;
	for(i = 0; i < 4; i++){
		if((loc = candidate->ps[i]) != -1){
			for(iter = contig->info[loc].alleles; iter; iter = iter->next){
				if(0 == strcmp(name,iter->read->individual) && 
				   iter->base == base){
				   printf("%d", iter->quality);
				   count++;
				   if(count < max){
				   		printf(",");
				   }
				}
			}
		}
	}
}

/* find and publish the alleles for this individual */
static inline void publish_profile(const contig* const contig,
								   const snp* const candidate,
								   const char* const name)
{
	int i, a=0, c=0, g=0, t=0;
	int loc;
	allele* iter;

	/* identify the alleles for this individual */
	for(i = 0; i < 4; i++){
		if((loc = candidate->ps[i]) != -1){
			for(iter = contig->info[loc].alleles; iter; iter = iter->next){
				if(0 == strcmp(name,iter->read->individual)){
					switch(iter->base){
						case A : a++; break;
						case C : c++; break;
						case G : g++; break;
						case T : t++; break;
						default  : break;
					}
				}
			}
		}
	}

	if(a > 0){
		printf("%s:A=%d(", name, a);
		print_quality_values(contig, candidate, name, a, A);
		printf(") ");
	}
	if(c > 0){
		printf("%s:C=%d(", name, c);
		print_quality_values(contig, candidate, name, c, C);
		printf(") ");
	}
	if(g > 0){
		printf("%s:G=%d(", name, g);
		print_quality_values(contig, candidate, name, g, G);
		printf(") ");
	}
	if(t > 0){
		printf("%s:T=%d(", name, t);
		print_quality_values(contig, candidate, name, t, T);
		printf(") ");
	}
}

/*print the position of this candidate allele*/
static void print_position(const contig* const contig,
						   const snp* const candidate,
						   hashtable* const individuals)
{
	int i;
	bin* iter;

	/* print the site of the candidate snp */
	int site = candidate->ps[0];
	printf("%d\t%d\t", mapping[site], mapping[site] + 1);

	for(i = 0; i < individuals->size; i++){
		iter = individuals->bins[i];
		while(iter){
			/* find and publish the alleles for this individual */
			publish_profile(contig, candidate, iter->name);
			
			iter = iter->next;
		}
	}
	printf("\n");
}

/* find a position where you have at least two reads which agree with a GAP.
 * Remember the consensus base at this position. Now look in the next
 * homopolymer stretch. Check to see if there is a position in there which has
 * some other base in the same reads where we had gaps earlier. Also the
 * consensus base should be a GAP here if the earlier consensus base was not a
 * GAP and vice versa*/
static void analyze_contig(contig* const contig, 
						   hashtable* const individuals,
						   const char* const contigname,
						   const char* const contignumreads)
{
	uint i;
	int j;

	snp* possible_snp;	/* candidate snp */ 
	snp* list = NULL;	/* list of snps for this contig */

	/* indexes which could be candidate snp locations */
	int* locations = ckalloc(4*sizeof(int));
	int max;

	/* coverage and number of gaps at this location */
	int cov, gaps;

	for(i = 0; i < (contig->num_bases - 1); i++){
		/*if any location has MIN_GAPS gaps or more, and more than MIN_HITS
		 * reads at that location, then we should check it out*/
		if((gaps = counts(contig->info[i], TRUE))    >= MIN_GAPS &&
		   (cov = slcount(contig->info[i].alleles)) >= MIN_HITS &&
		   (cov - gaps) >= MIN_GAPS){

			/*check the stretches of A,C,G,T,N to check for SNPs*/
			if(check_locations(contig, i, locations, &max)){
				possible_snp = ckalloc(sizeof(snp));
				for(j = 0; j < 4; j++){
					possible_snp->ps[j] = locations[j];
				}
				sladdhead(&list, possible_snp);
				i = locations[max];
			}
		}
	}
	
	/*lets create the mapping between positions in the ace file, fna file*/
	j = 1;
	mapping[0] = 0;
	if(contig->info[0].consensus == GAP){
		j--;
	}
	for(i = 1; i < contig->num_bases; i++){
		mapping[i] = contig->info[i].consensus == GAP ? mapping[i-1] : j++;
	}
	
	/*lets print the consensus for this contig*/
	for(i = 0, j = 0; i < contig->num_bases; i++){
		if(contig->info[i].consensus != GAP){
			if(contig->info[i].quality <= 20){
				buffer[j] = base_2_char(contig->info[i].consensus)+'a'-'A';
			}else{
				buffer[j] = base_2_char(contig->info[i].consensus);
			}
			j++;
		}
	}
	printf(">%s  length=%d  %s\n", contigname, j, contignumreads);
	format((uchar*)buffer,j);

	/*lets print the snps*/
	slreverse(&list);

	for(possible_snp = list; possible_snp; possible_snp = possible_snp->next){
		/* at least two options */
		if(possible_snp->ps[1] != -1){
			print_position(contig, possible_snp, individuals);
		}
	}

	slfreelist(&list);
	ckfree(locations);
}

/*print the contig sequence and all the candidate alleles in this contig*/
static void print_contig_alleles(char** const pfptr, 
							     size_t* const pn, 
							     FILE* const fp, 
							     hashtable* const hash,
								 hashtable* individuals,
								 const char* const contigname,
								 const char* const contignumreads)
{
	/*properties of a contig*/
	contig* contig;
	unsigned int num_bases, num_reads, num_segments;
	unsigned int reads_handled = 0;
	
	/*list of af segments for a contig*/
	af_seg* seg_list = NULL;
	af_seg* af_seg = NULL;

	/*variables to loop*/
	unsigned int index = 0, i;
	int j;
	char* ptr;
	ssize_t num_read;

	bin* hel;
	information info;
	position position;
	allele* allele;
	
	/*properties of a read*/
	sff_read* read = NULL;
	int offset;		
	char complement;
	unsigned int begin, end;

	while(getline(pfptr, pn, fp) != -1){
		/*lets get the information about the number of contigs*/
		if(strncmp(*pfptr, "AS", 2) == 0){
			continue;
		}

		/*is this a new contig*/
		if(strncmp(*pfptr, "CO", 2) == 0 &&
		   sscanf(*pfptr,"%*s %s %d %d %d ", 
		   		  gossip, &num_bases, &num_reads, &num_segments) == 4){
			if(num_bases > DSIZE){
				buffer = ckrealloc(buffer,  MAX(2*DSIZE,num_bases));
				mapping = ckrealloc(mapping,MAX(2*DSIZE,num_bases)*sizeof(int));
				DSIZE = MAX(2*DSIZE,num_bases);
			}
			contig = construct_contig(num_bases, num_reads, num_segments);
			contig->header = copy_string(gossip);
			
			/*initialize the af_segments for this contig*/
			if(seg_list != NULL){
				slfreelist(&seg_list);
			}

			/*read in the consensus sequence*/
			index = 0;
			while((num_read = getline(pfptr, pn, fp)) != -1 && num_read != 1){
				for(i = 0; i < (strlen(*pfptr)-1); i++){
					contig->info[index].consensus = char_2_base(*(*pfptr+i));
					index++;
				}
			}
		}

		/*the quality values for the consensus base of the contig*/
		if(strncmp(*pfptr,"BQ", 2) == 0){
			index = 0;
			while((num_read = getline(pfptr, pn, fp)) != -1 && num_read != 1){
				ptr = strtok(*pfptr, " ");
				while(ptr != NULL && *ptr != '\n'){
					while(contig->info[index].consensus == GAP){
						contig->info[index].quality = 0;
						index++;
					}
					contig->info[index].quality = atoi(ptr);
					index++;
					ptr = strtok(NULL, " \n");
				}
			}
		}

		/* we can get the offset in the contig for every read */
		if(strncmp(*pfptr,"AF", 2) == 0 &&
		   sscanf(*pfptr,"%*s %s %c %d\n",accession,&complement,&offset)==3){
		   	if((ptr = strchr(accession, '.')) != NULL){
				*ptr = '\0';
			}
			if((ptr = strchr(accession, '_')) != NULL){
				*ptr = '\0';
			}

			af_seg = ckallocz(sizeof(struct af_segment_st));
			af_seg->offset = offset-1;
			af_seg->complement = complement == 'C' ? TRUE : FALSE;
			sladdhead(&seg_list, af_seg);
			
			hel = lookup_hashtable(hash, accession, strlen(accession));
			assert(hel != NULL);
			read = hel->val;
			read->information->header = hel->name;
		}

		/*dont have any use for the BS segment*/
		if(strncmp(*pfptr,"BS", 2) == 0){
			slreverse(&seg_list);
			af_seg = seg_list;
			while((num_read = getline(pfptr, pn, fp)) != -1 && 
				   num_read != 1){
			}
		}

		/*lets fill in the remaning info from the RD section. The offset is for
		 *the start of the read in this section, the first base of the read 
		 *used in the alignment will not be known till we read the QA section.*/
		if(strncmp(*pfptr,"RD", 2) == 0 && 
		   sscanf(*pfptr,"%*s %s %d ", accession, &num_bases) == 2){
		   	if((ptr = strchr(accession, '.')) != NULL){
				*ptr = '\0';
			}
			if((ptr = strchr(accession, '_')) != NULL){
				*ptr = '\0';
			}

		   	read = must_find_hashtable(hash, accession,strlen(accession));
			read->information->num_bases = num_bases;
			if(read->information->sequence != NULL){
				ckfree(read->information->sequence);
			}
			read->information->sequence = 
								ckallocz(num_bases*sizeof(struct position_st));	
			
			index = 0;
			while((num_read = getline(pfptr, pn, fp)) != -1 && num_read != 1){
				for(i = 0; i < (strlen(*pfptr)-1); i++){
					read->information->sequence[index++].seq = *(*pfptr+i);
				}
			}
			assert(index == read->information->num_bases);
			ptr = (char*)read->data->quality_scores;
			if(read->information->complement){
				j = read->header->num_bases-1;
				for(i = 0; i < read->information->num_bases; i++){
					if(read->information->sequence[i].seq == '*'){
						read->information->sequence[i].qual = 0;
					}else{
						read->information->sequence[i].qual = ptr[j];
						j--;
					}
				}
				assert(j==-1);
			}else{
				j = 0;
				for(i = 0; i < read->information->num_bases; i++){
					if(read->information->sequence[i].seq == '*'){
						read->information->sequence[i].qual = 0;
					}else{
						read->information->sequence[i].qual = ptr[j];
						j++;
					}
				}
				assert(j == (int)read->header->num_bases);
			}	
			read->information->offset = af_seg->offset;
			read->information->complement = af_seg->complement;
			af_seg = af_seg->next;		

			contig->reads[reads_handled] = read;
			reads_handled++;
		}
	
		/*which of the bases in the read are being used in the alignment. Use
		 *them to create a profile of the whole contig*/	
		if(strncmp(*pfptr,"QA", 2) == 0 &&
		   sscanf(*pfptr, "%*s %d %d ", &begin, &end) == 2){
			/*0 based*/
			begin--;
			end--;
			
			/*these are the bases from the read used in the alignment*/
			read->information->begin = begin;
			read->information->end = end;
			index = begin + read->information->offset;
			for(i = begin; i <= end && index < contig->num_bases; i++){
				info = contig->info[index];
				position = read->information->sequence[i];
				allele = new_allele(position.seq, position.qual, read);
				sladdhead(&contig->info[index].alleles, allele);
				index++;
			}
		}

		/*dont have any use for the DS segment which is one line*/
		if(strncmp(*pfptr,"DS", 2) == 0){ 
			if(reads_handled == contig->num_reads){
				contig->next = NULL;
				/*analyze the contigs for snps and alleles*/
				analyze_contig(contig, individuals, contigname, contignumreads);
				
				if(Alignments){
					print_aligned_reads(contig->reads, 
										contig->num_reads, 
										individuals);
				}
#if DEBUG	
				print_contig(contig);
#endif
				free_contigs(&contig);
				if(seg_list){
					slfreelist(&seg_list);
				}
				break; 
			}
			continue;
		}
	}
}

/*print the details including the contig sequences and SNPs*/
static void print_exon_details(hashtable* const hash,
							   const char* const dir,
							   hashtable* const individuals)
{
	sprintf(buffer,"%s/exon/454AllContigs.fna", dir);
	sequence* sp = read_fasta_sequence(buffer);
	
	sprintf(buffer,"%s/exon/454Contigs.ace",dir);
	FILE* fp = ckopen(buffer,"r");		/*ACE file*/
	size_t n = 1;						
	char* fptr = ckalloc(n+1);			

	char name[1024];
	char numreads[1024];


	while(sp && strcmp((char*)sp->header,"") != 0){
		if(OnlyAssemble == TRUE){	
			printf(">%s\n", sp->header);
			format(sp->sequence, sp->slen);
		}else{
			if(sscanf((char*)sp->header,"%s %*s %s\n", name, numreads) != 2){
				fatalf("error in reading the name of the contig");
			}
			/*print the contig sequence candidate alleles*/
			print_contig_alleles(&fptr,&n,fp,hash,individuals,name,numreads);
		}
		
		if(!get_next_sequence(sp)){
			break;
		}
	}
	fclose(fp);
	close_fasta_sequence(sp);
	ckfree(fptr);
}

/*print the singletons for this exon*/
static void clone_header(common_header* const copy, 
					     const common_header* const header)
{
	copy->index_offset = 0;
	copy->index_length = 0;
	copy->num_reads = 0;

	copy->header_length = header->header_length;
	copy->key_length = header->key_length;
	copy->num_flows_per_read = header->num_flows_per_read;
	copy->flow_chars = copy_string(header->flow_chars);
	copy->key_sequence = copy_string(header->key_sequence);
}

/* go and fill in the common header section for the two technologies. Try to
 * iterate through the directory with the sff files and find the common header
 * section */
static void fill_headers(const char* const sff_dir,
						 common_header* const gs20,
						 common_header* const gsflx)
{
	DIR* dp;
	struct dirent* ep;
	FILE* fp;

	common_header* placeholder = ckallocz(sizeof(struct common_header_st));

	bool gs20_flag = TRUE;
	bool gsflx_flag = TRUE;

	if((dp = opendir(sff_dir)) == NULL){
		fatalf("error in listing the contents of the directory %s", sff_dir);
	}
	if (dp != NULL){
		while ((ep = readdir (dp)) != NULL){
			if(strncmp(ep->d_name,".",1)==0 || strncmp(ep->d_name,"..",2)==0){
				continue;
			}
			sprintf(buffer,"%s/%s", sff_dir, ep->d_name);
	
			fp = ckopen(buffer, "r");
			read_common_header(fp, placeholder);
			fclose(fp);
			
			if(gs20_flag && placeholder->num_flows_per_read == GS20){
				clone_header(gs20, placeholder);
				gs20_flag = FALSE;
			}else if(!gs20_flag && placeholder->num_flows_per_read == GS20){
				assert(gs20->header_length == placeholder->header_length);
				assert(gs20->key_length == placeholder->key_length);
//				assert(strcmp(gs20->flow_chars, placeholder->flow_chars) == 0);
				assert(strcmp(gs20->key_sequence,placeholder->key_sequence)==0);
			}
			
			if(gsflx_flag && placeholder->num_flows_per_read == GSFLX){
				clone_header(gsflx, placeholder);
				gsflx_flag = FALSE;
			}else if(!gsflx && placeholder->num_flows_per_read == GSFLX){
				assert(gsflx->header_length == placeholder->header_length);
				assert(gsflx->key_length == placeholder->key_length);
				assert(gsflx->num_flows_per_read == 
					   placeholder->num_flows_per_read);
//				assert(strcmp(gsflx->flow_chars,placeholder->flow_chars) == 0);
				assert(strcmp(gsflx->key_sequence,placeholder->key_sequence)==0);
			}
		
			ckfree(placeholder->flow_chars);
			ckfree(placeholder->key_sequence);
		
			if(!gs20 && !gsflx){	
				break;
			}
		}
		closedir (dp);
	}

	ckfree(placeholder);
}

/* free all resources associated with this read */
static void clear_reads(bin* bin)
{
	sff_read* read = bin->val;
	sff_read_free(&read);
	ckfree(read);
}

static void build_assemblies(const char* const exons_file,
							const char* const dir,
							common_header* gs20_header,
							common_header* gsflx_header)
{
	FILE* fp = ckopen(exons_file, "r");

	/*getline stuff*/
	size_t n = 1;
	ssize_t num_read = 0;
	char* fptr = ckalloc(n+1);

	/*a small hash to store the individuals*/
	hashtable* individuals = new_hashtable(3);	
	bin* hel;

	/*use this arguments for the assembly*/
	char assemble[512];
	sprintf(assemble,"runAssembly -ar -o %s/exon ", dir);

	/*the names of the gs20 sff file and the gsflx sff file*/
	char gs20[512], gsflx[512];
	sprintf(gs20,"%s/1.sff", dir);
	sprintf(gsflx,"%s/2.sff", dir);

	/*the name of the contigs file*/
	char contigsfile[512];
	sprintf(contigsfile,"%s/exon/454AllContigs.fna", dir);

	sff* gs20_sff = ckallocz(sizeof(sff));	/*exon with the gs20 reads*/
	sff* gsflx_sff = ckallocz(sizeof(sff));	/*exon with the flx reads*/
	sff_read* read = NULL;					/*single read*/
	char* header = NULL;					/*header for an exon*/
	struct stat sb;							/*to check if the file is empty*/
	
	gs20_sff->common_header = *gs20_header;
	gsflx_sff->common_header = *gsflx_header;	

	/* prepare to read the index file*/
	FILE* ip = ckopen(index_file, "r");
	if(ip == NULL){
		fatalf("error in reading the index from %s", index_file);
	}

	sffindex* rheader;
	rheader = sff_index_read_header (ip);
	if(rheader == NULL){
		fatalf("error in reading the header from the index %s", index_file);
	}

	sffft* ft;
    ft = sff_index_read_filenames (rheader);
	if(ft == NULL){
		fatalf("error in reading the file names from %s", index_file);
	}	
	
	FILE* sp;
	uint8_t fileNum;
	uint32_t filePos;		
	int flows_per_read;
	common_header* ch;

	hashtable* hash = NULL;

	while(getline(&fptr, &n, fp) != -1){
		if(strlen(fptr) > 1){
			header = copy_string(fptr);	
			read = NULL;
			gs20_sff->reads = NULL;
			gsflx_sff->reads = NULL;			
			hash = new_hashtable(4);
			
			/* read the reads for this exon */
			while((num_read = getline(&fptr, &n, fp)) != -1 && num_read != 1){
				if(sscanf(fptr,"%s %s\n", accession, gossip) != 2){ 
					fatalf("error in reading the read information:%s", fptr); 
				}

				if(sff_index_lookup_sequence (rheader, 
											  accession,
											  NULL,
											  &fileNum,
											  &filePos) != 1){
					fatalf("error in looking up the read %s", accession);
				}

				sp = ckopen(ft->name[fileNum], "r");
		
				ch = ckallocz(sizeof(struct common_header_st));
				read_common_header (sp, ch);
				flows_per_read = ch->num_flows_per_read;
				ckfree(ch->flow_chars);
				ckfree(ch->key_sequence);
				ckfree(ch);
		
				if(fseek(sp, filePos, SEEK_SET) != 0){
					fatalf("error in setting the file %s to %d", 
						    ft->name[fileNum], filePos);
				}

				read              = ckallocz (sizeof(sff_read));
				read->header      = ckallocz (sizeof(read_header));
				read->data        = ckallocz (sizeof(read_data));
				read->information = ckallocz (sizeof(ancilliary));
				read->header->num_flows_per_read = flows_per_read;
				read_sff_rinfo (&read, sp, flows_per_read);

				fclose(sp);
				read->next = NULL;
			
				add_hashtable(hash, accession, strlen(accession), read);

				/*which individual does this read belong to?*/
				if((hel = lookup_hashtable(individuals,
										   gossip,
										   strlen(gossip))) == NULL){
					hel = add_hashtable(individuals, 
										gossip, 
										strlen(gossip),
										NULL);
				}
				read->individual = hel->name;
					
				if(read->header->num_flows_per_read == GS20 && 
				   read->information->is_used == FALSE){
						sladdhead(&gs20_sff->reads,read);
						read->information->is_used = TRUE;
				}else if(read->header->num_flows_per_read == GSFLX && 
				         read->information->is_used == FALSE){
						sladdhead(&gsflx_sff->reads,read);
						read->information->is_used = TRUE;
				}
			}
			
			/*write the sff file for these reads*/				
			if(gs20_sff->reads != NULL && gsflx_sff->reads != NULL){
				sff_write_file(gs20, gs20_sff->common_header, gs20_sff->reads);
				sff_write_file(gsflx,gsflx_sff->common_header,gsflx_sff->reads);
				sprintf(command, "%s/%s %s %s > /dev/null",
								rig_directory,assemble, gs20, gsflx);
			}else if(gs20_sff->reads != NULL && gsflx_sff->reads == NULL){
				sff_write_file(gs20, gs20_sff->common_header, gs20_sff->reads);
				sprintf(command,
						"%s/%s %s > /dev/null",rig_directory,assemble, gs20);
			}else if(gs20_sff->reads == NULL && gsflx_sff->reads != NULL){
				sff_write_file(gsflx,gsflx_sff->common_header,gsflx_sff->reads);
				sprintf(command,
						"%s/%s %s > /dev/null",rig_directory,assemble, gsflx);
			}
			run_command(command);

			printf("%s", header);
			ckfree(header);

			/*print the exon details including the contigs and the SNPs*/
			if(0 == stat(contigsfile ,&sb) && sb.st_size != 0){
				print_exon_details(hash, dir, individuals);
				/*print the singletons if required*/
				if(Singletons){
					print_exon_singletons(hash, dir);
				}
			}

			printf("\n");

			/*lets clean up the mess*/
			sprintf(command,"/bin/rm -rf %s %s %s/exon", gs20, gsflx, dir);
			run_command(command);
			clear_used_flags(gs20_sff->reads);
			clear_used_flags(gsflx_sff->reads);
			func_hashtable(hash, clear_reads);
			free_hashtable(&hash);
		}
	}

	ckfree(ft);
	ckfree(rheader);

	fclose(fp);
	fclose(ip);

	free_hashtable_completely(&individuals);
	ckfree(fptr);	
}
							 

static void assemble1(const char* const sff_dir,  
					  const char* const exons_file)
{
	hashtable* hash = NULL;
	FILE* fp = ckopen(exons_file, "r");

	/* create a directory where this instance will run all its data. This will
	 * be really helpful when we decide to run multiple instances of assemble
	 * for the same dataset.*/
	char template[10] = "JOBXXXXXX";
	char* dir;
	if((dir = mkdtemp(template)) == NULL){
		fatalf("error in creating temp directory called %s : %s", 
				template, strerror(errno));
	}

	common_header* gs20 = ckallocz(sizeof(common_header));
	common_header* gsflx = ckallocz(sizeof(common_header));

	fill_headers(sff_dir, gs20, gsflx);
	fprintf(stderr,"Filled the headers for gs20 and gsflx\n");

	build_assemblies(exons_file, dir, gs20, gsflx);

	ckfree(gs20->flow_chars);	
	ckfree(gs20->key_sequence);	
	ckfree(gsflx->flow_chars);	
	ckfree(gsflx->key_sequence);	
	ckfree(gs20);
	ckfree(gsflx);

	free_hashtable(&hash);
	fclose(fp);

	/* remove the temporary directory */
	sprintf(command, "rm -rf %s", dir);
	run_command(command);
}

int main(int argc, char** argv)
{
	argv0="assemble1";
	/*read the command line arguments*/
	while(argc > 3){
		argc--;
		if(strncmp(argv[argc],"--newbler=", 10) == 0){
			rig_directory = argv[argc] + 10;
		}else if(strncmp(argv[argc],"-singletons",11) == 0){
			Singletons = TRUE;
		}else if(strncmp(argv[argc],"-alignments",11) == 0){
			Alignments = TRUE;
		}else if(strncmp(argv[argc],"-onlyassemble", 13) == 0){
			OnlyAssemble = TRUE;	
		}else if(strncmp(argv[argc],"--index=", 8) == 0){
			index_file = argv[argc] + 8;
		}else{
			fprintf(stderr,"Unknown option %s\n", argv[argc]);
			fatal(USE);
		}
	}

	if(index_file == NULL){
		fatalf("please input the name of the index file");
	}

	if(Alignments && OnlyAssemble){
		fatalf("Please select either of -alignments and -onlyassemble.");
	}
	
	if(argc != 3){
		fatal(USE);
	}

	/*allocate memory for the buffer*/
	buffer = ckalloc(DSIZE*sizeof(char));
	mapping = ckalloc(DSIZE*sizeof(int));
	command = ckalloc(DSIZE*sizeof(char));

	assemble1(argv[1], argv[2]);

	ckfree(buffer);
	ckfree(mapping);
	ckfree(command);
	return EXIT_SUCCESS;
}
