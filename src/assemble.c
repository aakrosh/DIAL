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

#define USE "assemble sff_directory exons.txt [-singletons] [-alignments] [-onlyassemble] [--newbler=directory_with_newbler]"

#define READS(A) ((A[0]) || (A[1]) || (A[2]) || (A[3]) || (A[4]))

/*number of flows per read for the gs20 sequences*/
#define GS20 400

/*number of flows per read for the gsflx sequences*/
#define GSFLX 800

/* number of flows per read for the lr800 sequences */
#define LR800 1400

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

/*command line arguments*/
bool Singletons     = FALSE;
bool Alignments     = FALSE;
bool OnlyAssemble   = FALSE;
bool Debug_flag     = FALSE;

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
	char check[6] = {0,0,0,0,0,0};
	valid_base base = N;

	valid_base consensus = call(contig->info[start]);
	check[consensus] = 1;

	if(consensus == N || consensus == GAP){
		return FALSE;
	}

	ploc[arrayptr++] = start;
	*pmax = *pmax + 1;

	uint end = start + 1;
	while(end < contig->num_bases && 
          contig->info[end].consensus == GAP &&
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

/*iterate through the different exon listings and build assemblies for them. Use
 * that information to call SNPs from them.*/
static void build_assemblies_for_exons(FILE* const fp,
									   const char* const dir,
									   const int toread, 
						        	   hashtable* const hash,
									   common_header* gs20_header,
									   common_header* gsflx_header,
                                       common_header* lr800_header)
{
	/*getline stuff*/
	size_t n = 1;
	ssize_t num_read = 0;
	char* fptr = ckalloc(n+1);

	int num_reads = 0; 				/*number of reads read*/

	/*a small hash to store the individuals*/
	hashtable* individuals = new_hashtable(3);	
	bin* hel;

	/*use this arguments for the assembly*/
	char assemble[512];
	sprintf(assemble,"runAssembly -ar -o %s/exon ", dir);

	/*the names of the gs20 sff file and the gsflx sff file*/
	char gs20[512], gsflx[512], lr800[512];
	sprintf(gs20,"%s/1.sff", dir);
	sprintf(gsflx,"%s/2.sff", dir);
    sprintf(lr800,"%s/3.sff", dir);

	/*the name of the contigs file*/
	char contigsfile[512];
	sprintf(contigsfile,"%s/exon/454AllContigs.fna", dir);

	sff* gs20_sff = ckallocz(sizeof(sff));	/*exon with the gs20 reads*/
	sff* gsflx_sff = ckallocz(sizeof(sff));	/*exon with the flx reads*/
    sff* lr800_sff = ckallocz(sizeof(sff)); /*exon with the lr800 reads */
	sff_read* read = NULL;					/*single read*/
	char* header = NULL;					/*header for an exon*/
	struct stat sb;							/*to check if the file is empty*/
	
	gs20_sff->common_header = *gs20_header;
	gsflx_sff->common_header = *gsflx_header;	
	lr800_sff->common_header = *lr800_header;

	while(getline(&fptr, &n, fp) != -1){
		if(strlen(fptr) > 1){
			header = copy_string(fptr);	
			read = NULL;
			gs20_sff->reads = NULL;
			gsflx_sff->reads = NULL;			
            lr800_sff->reads = NULL;

			/* read the reads for this exon */
			while((num_read = getline(&fptr, &n, fp)) != -1 && num_read != 1){
				if(sscanf(fptr,"%s %s\n", accession, gossip) != 2){ 
					fatalf("error in reading the read information:%s", fptr); 
				}

				read = must_find_hashtable(hash,accession,strlen(accession));
				read->next = NULL;
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
				}else if(read->header->num_flows_per_read == LR800 &&
                         read->information->is_used == FALSE){
						sladdhead(&lr800_sff->reads,read);
						read->information->is_used = TRUE;
                }
				num_reads++;
			}
			
			/*write the sff file for these reads*/				
			sff_write_file(gs20, gs20_sff->common_header, gs20_sff->reads);
			sff_write_file(gsflx, gsflx_sff->common_header, gsflx_sff->reads);
            sff_write_file(lr800, lr800_sff->common_header, lr800_sff->reads);
#if DEBUG
			check_integrity(gs20, gs20_sff->reads);
			check_integrity(gsflx,gsflx_sff->reads);
			check_integrity(lr800,lr800_sff->reads);
#endif
            int ptr = sprintf(command, "%s/%s ", rig_directory,assemble);
            if(gs20_sff->reads != NULL){
                ptr += sprintf(command + ptr, "%s ", gs20);
            }
            if(gsflx_sff->reads != NULL){
                ptr += sprintf(command + ptr, "%s ", gsflx);
            }
            if(lr800_sff->reads != NULL){
                ptr += sprintf(command + ptr, "%s ", lr800);
            }
            sprintf(command + ptr, "> /dev/null");

            /* replaced by the code above 
			if(gs20_sff->reads != NULL && gsflx_sff->reads != NULL){
				sprintf(command, "%s/%s %s %s > /dev/null",
								rig_directory,assemble, gs20, gsflx);
			}else if(gs20_sff->reads != NULL && gsflx_sff->reads == NULL){
				sprintf(command,
						"%s/%s %s > /dev/null",rig_directory,assemble, gs20);
			}else if(gs20_sff->reads == NULL && gsflx_sff->reads != NULL){
				sprintf(command,
						"%s/%s %s > /dev/null",rig_directory,assemble, gsflx);
			}
            */
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
			sprintf(command,"/bin/rm -rf %s %s %s %s/exon", gs20, gsflx, lr800, dir);
			if(!Debug_flag) run_command(command);
			clear_used_flags(gs20_sff->reads);
			clear_used_flags(gsflx_sff->reads);
			clear_used_flags(lr800_sff->reads);
		}
		
		if(num_reads >= toread){
			break;
		}
	}

	free_hashtable_completely(&individuals);
	ckfree(fptr);	
}

/*read the next BATCH names of the reads in the exons_file and add them 
 *to the hash. Return the number of names read */
static int read_names(FILE* const fp , hashtable* const hash)
{
	size_t n = 1;
	char* fptr = ckalloc(n+1);

	int read = 0;

	while(getline(&fptr, &n, fp) != -1){
		while(getline(&fptr, &n, fp) != -1 &&
		      strlen(fptr) != 1            &&
			  sscanf(fptr,"%s ", buffer) == 1){
			if(lookup_hashtable(hash, buffer,strlen(buffer)) == NULL){
				add_hashtable(hash, buffer, strlen(buffer),NULL);
			}
			read++;
		}
		if(read > BATCH){
			break;
		}
	}
	
	ckfree(fptr);
	fprintf(stderr,"Read %d accession IDs\n", read);
	return read;
}

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

/*go and read the sff files in the argv[1] directory and if the name of a
 * read is in the hash, then we should hash that sff_read as the value of
 * that read in the hash. We would need two different common_headerS because
 * we have two different flow types for the reads. */
static void read_sff_files(const char* const sff_dir, 	
						   hashtable* const hash,
						   common_header* const gs20, 
						   common_header* const gsflx,
                           common_header* const lr800)
{
	DIR* dp;
	struct dirent* ep;
	common_header* sequencer = NULL;

	common_header* placeholder = ckallocz(sizeof(common_header));

	bool gs20_flag = TRUE;
	bool gsflx_flag = TRUE;
    bool lr800_flag = TRUE;

	if((dp = opendir(sff_dir)) == NULL){
		fatalf("error in listing the contents of the directory %s", sff_dir);
	}
	if (dp != NULL){
		while ((ep = readdir (dp)) != NULL){
			/*read the contents of this sff file and if the name of the read
			 * exists in the hash, then keep the sff_read structure for that
			 * read.Also make sure that one common_header for each of the gs20
			 * and flx runs is stored.*/
			if(strncmp(ep->d_name,".",1)==0 || strncmp(ep->d_name,"..",2)==0){
				continue;
			}
			fprintf(stderr,"Reading %s/%s\n", sff_dir,ep->d_name);
			sprintf(buffer,"%s/%s", sff_dir, ep->d_name);
			sequencer = sff_read_file(buffer, hash);

			clone_header(placeholder, sequencer);

			if(gs20_flag && placeholder->num_flows_per_read == GS20){
				clone_header(gs20, sequencer);
				gs20_flag = FALSE;
			}else if(!gs20_flag && placeholder->num_flows_per_read == GS20){
				assert(gs20->header_length == sequencer->header_length);
				assert(gs20->key_length == sequencer->key_length);
				assert(strcmp(gs20->flow_chars, sequencer->flow_chars) == 0);
				assert(strcmp(gs20->key_sequence, sequencer->key_sequence)==0);
			}
			
			if(gsflx_flag && placeholder->num_flows_per_read == GSFLX){
				clone_header(gsflx, sequencer);
				gsflx_flag = FALSE;
			}else if(!gsflx && placeholder->num_flows_per_read == GSFLX){
				assert(gsflx->header_length == sequencer->header_length);
				assert(gsflx->key_length == sequencer->key_length);
				assert(gsflx->num_flows_per_read==sequencer->num_flows_per_read);
				assert(strcmp(gsflx->flow_chars, sequencer->flow_chars) == 0);
				assert(strcmp(gsflx->key_sequence, sequencer->key_sequence)==0);
			}

            if(lr800_flag && placeholder->num_flows_per_read == LR800){
				clone_header(lr800, sequencer);
                lr800_flag = FALSE;
            }else if(!lr800_flag && placeholder->num_flows_per_read == LR800){
                assert(lr800->header_length == sequencer->header_length);
                assert(lr800->key_length == sequencer->key_length);
                assert(lr800->num_flows_per_read==sequencer->num_flows_per_read);
                assert(strcmp(lr800->flow_chars, sequencer->flow_chars) == 0);
                assert(strcmp(lr800->key_sequence, sequencer->key_sequence)==0);
            }

            if((gs20_flag == TRUE) && (gsflx_flag == TRUE) && (lr800_flag == TRUE)){
                fprintf(stderr, "This SFF file has %d flows per read, whereas DIAL only supports %d, %d and %d flows per read. Please run\n\tsfffile -20|-flx|-xlr %s\n to convert the flowgrams into GS/FLX/XLR cycles and then add them to the project.\n", placeholder->num_flows_per_read, GS20, GSFLX, LR800, ep->d_name);
                exit(EXIT_FAILURE);
            }
			
			ckfree(sequencer->flow_chars);
			ckfree(sequencer->key_sequence);
			ckfree(sequencer);
		}
		closedir (dp);
	}

	ckfree(placeholder->flow_chars);
	ckfree(placeholder->key_sequence);
	ckfree(placeholder);
}

/* free all resources associated with this read */
static void clear_reads(bin* bin)
{
	sff_read* read = bin->val;
	sff_read_free(&read);
	ckfree(read);
}

static void assemble(const char* const sff_dir,  
					 const char* const exons_file)
{
	hashtable* hash = NULL;
	FILE* fp = ckopen(exons_file, "r");

	long fpos1, fpos2 = 0;	/*the offset in the exons file for this iteration*/
	int num_read;			/*number of exons reads in this iteration*/

	/* create a directory where this instance will run all its data. This will
	 * be really helpful when we decide to run multiple instances of assemble
	 * for the same dataset.*/
	char template[10] = "JOBXXXXXX";
	char* dir;
	if((dir = mkdtemp(template)) == NULL){
		fatalf("error in creating temp directory called %s : %s", 
				template, strerror(errno));
	}

	while(1){
		hash = new_hashtable(20);
		
		fpos1 = fpos2;
		
		/*read the next BATCH names of the reads in the exons_file and 
		 *include them in the hash*/
		if((num_read = read_names(fp, hash)) == 0){
			break;
		}

		fpos2 = ftell(fp);
	
		/*go and read the sff files in the argv[1] directory and if the name 
		 *of a read is in the hash, then we should hash that sff_read as the 
		 *value of that read in the hash. We would need two different 
		 *common_headerS because we have two different flow types for the 
		 *reads. */
		common_header* gs20 = ckallocz(sizeof(common_header));
		common_header* gsflx = ckallocz(sizeof(common_header));
        common_header* lr800 = ckallocz(sizeof(common_header));
		read_sff_files(sff_dir, hash, gs20, gsflx, lr800);

		/* reposition the stream to the old place */
		if(fseek(fp, fpos1, SEEK_SET) != 0){
			fatalf("error in repositioning the file");
		}

		/* iterate through the exons files building assemblies for each of 
		 * them and then finding SNPs*/
		build_assemblies_for_exons(fp, dir, num_read, hash, gs20, gsflx, lr800);
		
		
		/* free all these resources for the next round */
		ckfree(gs20->flow_chars);	
		ckfree(gs20->key_sequence);	
		ckfree(gsflx->flow_chars);	
		ckfree(gsflx->key_sequence);	
        ckfree(lr800->flow_chars);
        ckfree(lr800->key_sequence);
		ckfree(gs20);
		ckfree(gsflx);
        ckfree(lr800);

		func_hashtable(hash, clear_reads);
		free_hashtable(&hash);
	}

	free_hashtable(&hash);
	fclose(fp);

	/* remove the temporary directory */
	sprintf(command, "rm -rf %s", dir);
    if(!Debug_flag) run_command(command);
}

int main(int argc, char** argv)
{
	argv0="assemble";
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
        }else if(strncmp(argv[argc],"-debug", 6) == 0){
            Debug_flag = TRUE;
		}else{
			fprintf(stderr,"Unknown option %s\n", argv[argc]);
			fatal(USE);
		}
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

	assemble(argv[1], argv[2]);

	ckfree(buffer);
	ckfree(mapping);
    ckfree(command);
	return EXIT_SUCCESS;
}
