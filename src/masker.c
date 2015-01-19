/* This module takes as input a fasta file of reads :
 * 		masker 454.fa
 * The output of this module is a file of masked reads */

#include <ctype.h>
#include <math.h>
#include "utilities.h"
#include "common.h"
#include "sequences.h"
#include "splayhash.h"

#define USE "masker reads.fa"

/* Number of kmers in a window */
int WindowSpan = 5;

/* percentiles */
#define THRESHOLD 0.995
#define EXTEND 0.99
#define LOW 0.90
#define HIGH 0.998

typedef struct frequecy_st
{
	kmer kmer;
	int count;
}frequency;

int muls[MAX_READ_SIZE];

char buffer[128];

#if DEBUG
/* convert the kmer into a string. The output is written in "buffer" */
static void stringify(const kmer name, const int klen)
{
	kmer kmer = name;
	
	int i = klen;
	while(i >= 0){
		buffer[--i] = bit_encoding[kmer & 3];
		kmer = kmer >> 2;
	}
}

static void print_kmers_tree(splaybin* bin, const int klen)
{
	stringify(bin->name, klen);
	printf("%s\t%d\n", buffer, bin->counts);
	kmer rc = reverse_complement_kmer(bin->name, klen);
	stringify(rc, klen);
	if(rc != bin->name){
		printf("%s\t%d\n", buffer, bin->counts);
	}
	
	if(bin->left != NULL){
		print_kmers_tree(bin->left, klen);
	}
	if(bin->right != NULL){
		print_kmers_tree(bin->right, klen);
	}
}

static void print_kmers(splayhash* const hash, const int klen)
{
	int i;
	for(i = 0; i < hash->size; i++){
		if(hash->bins[i] != NULL){
			print_kmers_tree(hash->bins[i], klen);
		}
	}
}
#endif

/* find the size of the kmers to be used */
static int find_kmer_size(const char* const readsfile)
{
	sequence* sp;
	if((sp = read_fasta_sequence(readsfile)) == NULL){
		fatalf("error in reading the sequence from %s", readsfile);
	}

	long int total = 0;
	while(sp){
		total += sp->slen;

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);

	return ceil(log(5*total)/log(4));
}

static void process_sequence(const sequence* const sp, 
							 splayhash* const hash, 
							 const int klen)
{
	uint i;
	kmer word, antiword, selection;
	splaybin* bin;
	const uchar* s = sp->sequence;

	/*the first kmer in the string*/
	word = build_index(s, klen);
	antiword = reverse_complement_kmer(word, klen);
	selection = word < antiword ? word : antiword;
	int count = word == antiword ? 2 : 1;
	if((bin = lookup_splayhash(hash, selection)) == NULL){
		add_splayhash(hash, selection, count, 0);
	}else{
		bin->counts += count;
	}

	/*process the remaining kmers*/
	for(i = 1; i <= (sp->slen - klen); i++){
		word = get_next_kmer(word, s, klen, i);
		antiword = reverse_complement_kmer(word, klen);
		selection = word < antiword ? word : antiword;
		count = word == antiword ? 2 : 1;
		if((bin = lookup_splayhash(hash, selection)) == NULL){
			add_splayhash(hash, selection, count, 0);
		}else{
			bin->counts += count;
		}
	}
}

/* calculate the frequencies of the kmers. Also tell us how many kmers have at
 * least once occurence in the reads */
static splayhash* find_kmer_frequencies(const char* const readsfile, 
								  const int klen, 
								  int* const unique)
{
	splayhash* hash = new_splayhash(16);
	*unique = 0;

	sequence* sp;
	if((sp = read_fasta_sequence(readsfile)) == NULL){
		fatalf("error in reading the sequences from %s", readsfile);
	}

	while(sp){
		if(sp->slen >= (uint)klen){
			process_sequence(sp, hash, klen);
		}

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);

	*unique = hash->elcount;

	return hash;
}

static void add_frequencies(splaybin* bin, frequency** const cum, int* a)
{
	frequency* f = ckalloc(sizeof(frequency));
	f->kmer = bin->name;
	f->count = bin->counts;
	cum[*a] = f;
	*a = *a + 1;

	if(bin->left != NULL){
		add_frequencies(bin->left, cum, a);
	}
	if(bin->right != NULL){
		add_frequencies(bin->right, cum, a);
	}
}

static void fill_frequencies(splayhash* const freq, 
							 frequency** const cumulative)
{
	int i, j;
	for(i = 0, j = 0; i < freq->size; i++){
		if(freq->bins[i] != NULL){
			add_frequencies(freq->bins[i], cumulative, &j);
		}
	}
}

static int compar(const void* const a, const void* const b)
{
	return (*(frequency**)a)->count - (*(frequency**)b)->count;
}

static void calculate_cutoffs(splayhash* const freq,
							  const int unique,
							  int* const threshold,
							  int* const extend,
							  int* const low,
							  int* const high)
{
	/* sort the kmer based on their counts */
	frequency** cumulative = ckalloc(unique*sizeof(frequency*));
	fill_frequencies(freq, cumulative);
	qsort(cumulative, unique, sizeof(frequency*), compar); 

	/* construct the cdf */
	int i, total = 0;
	for(i = 0; i < unique; i++){
		total += cumulative[i]->count;
	}

	/* find the thresholds */
	int thresholdpc = THRESHOLD * unique;
	int extendpc    = EXTEND  * unique;
	int lowpc       = LOW  * unique;
	int highpc      = HIGH * unique;

	*low = cumulative[lowpc]->count;
	*extend = cumulative[extendpc]->count;
	*threshold = cumulative[thresholdpc]->count;
	*high = cumulative[highpc]->count;

	for(i = 1; i < unique; i++){
		ckfree(cumulative[i]);
	}
	ckfree(cumulative);
}

static void assign_scores_tree(splaybin* bin, const int low, const int high)
{
	if(bin->counts > high){
		bin->score = high;
	}else if(bin->counts <= low){
		bin->score = ceil(low/2);
	}else{
		bin->score = bin->counts;
	}

	if(bin->left != NULL){
		assign_scores_tree(bin->left, low, high);
	}
	if(bin->right != NULL){
		assign_scores_tree(bin->right, low, high);
	}
}

static void assign_scores(splayhash* const hash,
						  const int low,
						  const int high)
{
	int i;
	for(i = 0; i < hash->size; i++){
		if(hash->bins[i] != NULL){
			assign_scores_tree(hash->bins[i], low,  high);
		}
	}
}

/* mask the 0-based half open interval on the sequences */
static int mask_chars(sequence* const sp, const int start, const int end)
{
	int i, j = 0;
	for(i = start ; i < end; i++){
		if(sp->sequence[i] < 'a'){
			j++;
			sp->sequence[i] = sp->sequence[i] + ('a' - 'A');
		}
	}

	return j;
}

static int  mask_sequence(sequence* const sp,
						  const int klen,
						  splayhash* const hash,
						  const int extend, 
						  const int threshold)
{
	int i = 0, score = 0, j = 0, k;
	kmer word, antiword, selection;
	splaybin* bin;

	int scores[MAX_READ_SIZE];
	int winscores[MAX_READ_SIZE];
	uchar* s = sp->sequence;
	
	/*the first kmer in the string*/
	word = build_index(s, klen);
	antiword = reverse_complement_kmer(word, klen);
	selection = word < antiword ? word : antiword;
	if((bin = lookup_splayhash(hash, selection)) == NULL){
		fatalf("kmer %"PRIu64" should be present\n", selection);
	}else{
		score += bin->score;
		scores[0] = bin->score;
	}

	/*process the remaining kmers*/
	for(i = 1; i <= (int)(sp->slen - klen); i++){
		/* find the score for this window */
		if(i >= WindowSpan){
			winscores[j++] = score / WindowSpan;
			score -= scores[i - WindowSpan];
		}

		word = get_next_kmer(word, s, klen, i);
		antiword = reverse_complement_kmer(word, klen);
		selection = word < antiword ? word : antiword;
		if((bin = lookup_splayhash(hash, selection)) == NULL){
			fatalf("kmer %"PRIu64" should be present\n", selection);
		}else{
			score += bin->score;
			scores[i] = bin->score;
		}
	}

	/* which windows should I mask? */
	int last = 0;
	bool maskable;
	int num_masked = 0;

	for(i = 0; i < j; i++){
		if(winscores[i] > threshold){
			/* mask this window */
			num_masked += mask_chars(sp, i, i+klen+WindowSpan-1);

			/*check every window from index 'last' to 'i'*/
			maskable = TRUE;
			for(k = last; k < i; k++){
				if(winscores[k] <= extend){
					maskable = FALSE;
					break;
				}
			}
			if(maskable == TRUE){
				num_masked += mask_chars(sp, last, i);
			}
			last = i + 1;
		}
	}

	print_fasta_sequence(sp);
	return num_masked;
}

static void add_to_distributions(splaybin* bin)
{
	muls[bin->counts]++;
	if(bin->left != NULL){
		add_to_distributions(bin->left);
	}
	if(bin->right != NULL){
		add_to_distributions(bin->right);
	}
}

static void create_distribution(splayhash* const hash)
{
	int i;
	for(i = 0; i < hash->size; i++){
		if(hash->bins[i] != 0){
			add_to_distributions(hash->bins[i]);
		}
	}
}

static bool check_read(const sequence* const sp, 
					   const int klen,
					   const int extend)
{	
	bool result = FALSE;
	splayhash* hash = new_splayhash(8);
	
	/*find all the kmers in this read */
	if(sp->slen >= (uint)klen){
		process_sequence(sp, hash, klen);
	}else{
		return TRUE;
	}

	/* create the distribution */
	memset(muls, 0, MAX_READ_SIZE*sizeof(int));
	create_distribution(hash);

	int i;
	for(i = extend; i < MAX_READ_SIZE; i++){
		if(muls[i] != 0){
			result = TRUE;
			break;
		}
	}

	free_splayhash(&hash);
	return result;
}

static void mask_sequences(const char* const readsfile, 
						   splayhash* const hash,
						   const int klen,
						   const int extend,
						   const int threshold)
{
	sequence* sp;
	if((sp = read_fasta_sequence(readsfile)) == NULL){
		fatalf("error in reading the sequences from %s", readsfile);
	}

	int num_masked = 0, total = 0;
	bool mask_completely;
	while(sp){
		total += sp->slen;
		
		/*lets check the distribution of kmers in this read. If any of its 
		 * kmers occur more than "extend" times throughout the read, 
		 * then we just mask the read */
		if((mask_completely = check_read(sp, klen,extend)) == TRUE){
			num_masked += mask_chars(sp, 0, sp->slen);
			print_fasta_sequence(sp);
		}else{
			num_masked += mask_sequence(sp, klen, hash, extend, threshold);
		}

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);
#if DEBUG
	fprintf(stderr,"%3.2f%% read bases masked\n", (num_masked*100.0)/total);
#endif
}

static void masker(const char* const readsfile)
{
	/* find the size of the kmers to be used */
	int klen = find_kmer_size(readsfile);
#if DEBUG
	fprintf(stderr,"Masking and Filtering details\n");
	fprintf(stderr,"-----------------------------\n");
	fprintf(stderr,"Size of kmers to be used : %d\n", klen);
#endif

	/* calculate the frequencies of the kmers */
	int unique;
	splayhash* freq = find_kmer_frequencies(readsfile, klen, &unique);
#if DEBUG
	fprintf(stderr,"Calculated the frequency of the kmers\n");
	fprintf(stderr,"%d unique kmers used in the reads\n", 2*unique);
	print_kmers(freq, klen);
#endif

	/* calculate the various cutoof values */
	int threshold, extend, low, high;
	calculate_cutoffs(freq, unique, &threshold, &extend, &low, &high);
#if DEBUG
	fprintf(stderr,"Low count: %d\t", low);
	fprintf(stderr,"Extend count: %d\t", extend);
	fprintf(stderr,"Threshold count: %d\t", threshold);
	fprintf(stderr,"High count: %d\n\n", high);
#endif

	/* assign scores for each of the kmers */
	assign_scores(freq, low, high);

	/* mask the sequences */
	mask_sequences(readsfile, freq, klen, extend, threshold);
}

int main(int argc, char** argv)
{
	argv0="masker";

	while(argc > 2){
		argc--;
		if(strncmp(argv[argc], "--window-span=", 14) == 0){
			WindowSpan = atoi(argv[argc] + 14);
		}else{
			fatalf("unknown argument:%s\n%s",argv[argc],USE);
		}
	}

	if(argc != 2){
		fatal(USE);
	}

	allocate_resources();

	masker(argv[1]);

	return EXIT_SUCCESS;
}
