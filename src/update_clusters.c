/* Module to update the clusters for the half runs */
#include <math.h>

#include "utilities.h"
#include "hashtable.h"
#include "maf.h"
#include "cluster.h"
#include "mfblocks.h"

#define USE "update_clusters new.clt old.clt < fake.maf"

/* minimum length of the required overlap */
uint MinLength = 100;
/* number of bases added in the project till now */
long NumBases = 0;
/* expected size of the genome */
long ExpGenome = 0;
/* if this is a transcript dataset, then we should not set an upper limit on the
 * number of hits */
bool TranscriptDataset = FALSE;

static hashtable* read_clusters(const char* const file)
{
	hashtable* hash = new_hashtable(16);

	cluster* cp;
	cluster* copy;

	if((cp = read_cluster_file(file)) == NULL){
		fatalf("error in reading the clusters from %s", file);
	}

	while(cp){
		copy = copy_cluster(cp);
		assert(copy->num_alloced == copy->num_members);
		add_hashtable(hash, copy->name, strlen(copy->name), copy);

		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);

	return hash;
}

/* simple function to return the factorial of an int */
static long unsigned int factorial(int n) {
    return n>=1 ? n * factorial(n-1) : 1;
}

/* simple function to calculate the maxhits that should be allowed */
static int calculate_maxhits(const float coverage)
{
    // if the coverage is too high, I will not be able to calculate the
    // factorial correctly without a big number library, but this should be a
    // decent approximation for our purposes
    if(coverage > 19){
        return 2 * coverage + 1;
    }

	float total = 0;
	float probability;
	int maxhits = 0;
	while(total < 0.9999){
		probability = pow(coverage, maxhits)*exp(-coverage)/factorial(maxhits);
		total += probability;
		maxhits++;
	}
	return maxhits+1;
}

static void process_blocks(hashtable* const hash,
						   cluster* const cp,
						   mfblocks** const pblocks,
						   const uint maxhits,
						   const char* const currentname,
						   const char* const oldname)
{
	mfblocks* iter = *pblocks;		/* iterate on the mfblocks */

	mafblock* block;
	cluster* tp;
	uint i, j;
	diffs diff;

	while(iter){
		block = iter->block;
		if(block->o1 > MinLength && block->o2 > MinLength){
			assert(strcmp(cp->name, (char*)block->name2) == 0);
			if(cp->num_members <= (maxhits + 1)){
			   	addmember(cp, (char*)block->name1, currentname);
				for(i = 0, j = block->s2; i < block->o2; i++, j++){
					if(block->aln1[i] != block->aln2[i] && 
					   block->aln1[i] < 'a' && 
					   block->aln2[i] < 'a' && 
                       block->aln1[i] != 'N' &&
                       block->aln2[i] != 'N'){
						diff.position = j;
						diff.reference = block->aln2[i];
						diff.target = block->aln1[i];
						adddifference(cp, diff);
					}
				}
			}

			/* add it to the other cluster */
			tp = must_find_hashtable(hash, 
									 (char*)block->name1, 
									 strlen((char*)block->name1));
			if(tp->num_members <= (maxhits + 1)){
				addmember(tp, (char*)block->name2, oldname);
				for(i = 0, j = block->s1; i < block->o1; i++, j++){
					if(block->aln1[i] != block->aln2[i] && 
					   block->aln1[i] < 'a' && 
					   block->aln2[i] < 'a' &&
                       block->aln1[i] != 'N' && 
                       block->aln2[i] != 'N') {
						diff.position = j;
						diff.reference = block->aln1[i];
						diff.target = block->aln2[i];
						adddifference(tp, diff);
					}
				}
			}
		}

		iter = iter->next;
	}
}

/* process the maf file to find new overlaps */
static void process_maf(hashtable* const currentrun,
						const char* const oldrun,
						const int maxhits,
						const char* const currentname,
						const char* const oldname)
{
	cluster* cp;
	if((cp = read_cluster_file(oldrun)) == NULL){
		fatalf("unable to read the clusters from %s", oldrun);
	}

	mafblock* mp;
	if((mp = read_maf_file("/dev/stdin")) == NULL){
		fprintf(stderr, "no putative overlaps found\n");
	}

	mfblocks* blocks;
	mfblocks* mfnode;

	while(cp){
		blocks = NULL;

		/* is the current maf block for this read? */
		if(mp && 
		   strncmp(cp->name, (char*)mp->name2, strlen((char*)mp->name2)) == 0){
			/* fetch all the remaining maf blocks for this read */
			do{
				mfnode = new_mfnode(mp);
				sladdhead(&blocks, mfnode);		
			}while(mp && 
			 get_next_block(mp) && 
			 strncmp(cp->name,(char*)mp->name2, strlen((char*)mp->name2)) == 0);
		}
		
		if(blocks != NULL){
			/* process all the maf blocks for alignments for this cluster */
			process_blocks(currentrun,cp,&blocks,maxhits,currentname,oldname);

			/* free  all the memory */
			free_mfblocks(&blocks);
		}

		print_cluster(cp, oldname, maxhits, stdout);

		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);
}

/* print the clusters in the hashtable in the same order as that in the 
 * seqfile. */
static void print_updated_cluster(hashtable* const hash, 
								  const int maxhits,
								  const char* const seqfile,
								  const char* const name)
{
	sequence* sp;
	if((sp = read_fasta_sequence(seqfile)) == NULL){	
		fatalf("error in reading the sequences from %s", seqfile);
	}

	cluster* cp;
	char* key;
	while(sp){
		key = (char*)sp->header;
		cp = must_find_hashtable(hash, key, strlen(key));
		print_cluster(cp, name, maxhits, stdout);

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);
}

static void update_clusters(const char* const currentrun,
						const char* const currentreads,
						const char* const currentname,
						const char* const oldrun,
						const char* const oldreads __attribute__((unused)),
						const char* const oldname)
{
	/* use the number of bases and expected genome size to calculate the 
	 * maximum number of members that should be allowed in a cluster */
	float coverage = (NumBases*1.0)/ExpGenome;
	int maxhits = TranscriptDataset == TRUE 
				? INT_MAX - 100 : calculate_maxhits(coverage); 
#if DEBUG
	fprintf(stderr,"Maximum number of hits in a cluster:%d\n", maxhits);
#endif

	/* hash the clusters for the current run */
	hashtable* hash = read_clusters(currentrun);

	/* process maf output from lastz */
	process_maf(hash, oldrun, maxhits, currentname, oldname);

	/* print the clusters for the current run */
	printf("----------------------------------------------------\n");
	print_updated_cluster(hash, maxhits, currentreads, currentname);
}

int main(int argc, char** argv)
{
	while(argc > 7){
		argc--;
		if(strncmp(argv[argc], "--minlen=", 9) == 0){
			MinLength = atoi(argv[argc] + 9);
		}else if(strncmp(argv[argc],"--numbases=", 11) == 0){
			NumBases = strtol(argv[argc] + 11, NULL, 10);
		}else if(strncmp(argv[argc],"--expsize=", 10) == 0){
			ExpGenome = strtol(argv[argc] + 10, NULL, 10);
		}else if(strncmp(argv[argc], "-transcript",11) == 0){
			TranscriptDataset = TRUE;
		}else{
			fatalf("Unknown argument:%s\n%s", argv[argc], USE);
		}
	}

	if(ExpGenome == 0 || NumBases == 0){
		fatalf("We need to know the number of bases in the reads && expected size of the target genome.");
	}

	update_clusters(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	return EXIT_SUCCESS;
}
