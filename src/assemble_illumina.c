/* Simple module which takes in two arguments:
 * 	reads_directory : the directory which contains the reads and quality files
 *  clusters.txt : the file which has all the clusters and the components
 *  listed
 *The output is a simple file which has the fastq sequences for the reads*/

#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "utilities.h"
#include "sequences.h"
#include "hashtable.h"

#define USE "assemble_illumina reads_directory clusters.txt"

#define DSIZE 4096

/* the buffer to fill in the shell commands */
char* command = NULL;

/* the place where the Velvet binaries are kept */
char* rig_directory = "/usr/local/rig/bin";

/* command line arguments*/
bool Alignments = FALSE;
bool OnlyAssemble = FALSE;

const unsigned char rcs[] = 
  "                                                                "
  " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
  "                                                                "
  "                                                                ";

/* a single fastq sequence, consists of the bases and the quality values*/
typedef struct fsequence_st
{
	char* sequence;
	char* quality;
	uint start;			/* 0-based start position on the read used in aligment*/
	uint end;			/* 0-based end position on the read used in alignment */
	int offset;			/* 0-based offset of the read on the contig */
	char* person;		/* name of the individual, the read belongs to */
}fsequence;

typedef struct allele_st
{
	struct allele_st* next;
	char base;
	uint8_t quality;
	fsequence* read;
}allele;

typedef struct information_st
{
	char consensus;
	allele* alleles;
}information;

typedef struct contig_st
{
	struct contig_st* next;
	
	/*details of the reads in this contig */
	int num_reads;
	fsequence** reads;
	
	/* details of the bases and quality of the bases in the contig*/
	int num_bases;
	information* info;
}contig;

/* struct to temporary hold the details of the reads */
typedef struct tle_st
{
	struct tle_st* next;
	uint id;
	int offset;
	uint start;
	uint end;
}tle;

/* read the names of the reads that are members of clusters, in this file*/
static hashtable* read_names(const char* const filename)
{
	size_t n = 1;
	char* fptr = ckalloc(n+1);
	FILE* fp = ckopen(filename, "r");

	hashtable* hash = new_hashtable(10);

	int read = 0;
	char name[1024];

	while(getline(&fptr, &n, fp) != -1){
		while(getline(&fptr, &n, fp) != -1 &&
		      strlen(fptr) != 1            &&
			  sscanf(fptr,"%s ", name) == 1){
			if(lookup_hashtable(hash, name,strlen(name)) == NULL){
				add_hashtable(hash, name, strlen(name),NULL);
			}
			read++;
		}
	}
	
	ckfree(fptr);
	fclose(fp);
	fprintf(stderr,"Read %d accession IDs\n", read);
	return hash;
}

/* only consider this file if the substring ".fa" exists in the name of the 
 * file*/
int selector(const struct dirent* dp)
{
	if(strstr(dp->d_name, ".fa") == NULL){
		return 0;
	}

	return 1;
}

static void read_files(const char* const directory,
					   const char* const filename,
					   hashtable* const names)
{
	char* prefix = strdup(filename);

	char* ptr;
	if((ptr = strstr(prefix, ".fa")) == NULL){
		fatalf("expect the name of the sequence files to end with a \"fa\"");
	}
	*ptr = 0;

	/*form the name of the sequence file and read it to select the required
	 * reads from the hash*/
	char sequencefile[1024];
	sprintf(sequencefile, "%s/%s.fa", directory, prefix);

	sequence* sp;
	if((sp = read_fasta_sequence(sequencefile)) == NULL){
		fatalf("error in reading the sequences from %s", sequencefile);
	}

	fsequence* fseq;
	bin* bin;
	while(sp){
		if((bin = lookup_hashtable(names,
						          (char*)sp->header,
						          strlen((char*)sp->header))) != NULL){
			fseq = ckallocz(sizeof(fsequence));
			fseq->sequence = copy_string((char*)sp->sequence);
			fseq->quality = ckallocz(strlen(fseq->sequence)*sizeof(char));
			bin->val = fseq;
		}
		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);

	/* lets try reading the quality values from the files */
	sprintf(sequencefile, "%s/%s.qual", directory, prefix);

	size_t n = 1;
	char* fptr = ckalloc(n+1);
	FILE* fp = ckopen(sequencefile, "r");

	if(getline(&fptr, &n, fp) == -1){
		fatalf("error in reading from the file: %s", sequencefile);
	}

	int index;
	while(1){
		if(fptr[0] == '>'){
			if((bin = lookup_hashtable(names, 
									   fptr + 1, 
									   strlen(fptr+1)-1)) != NULL){
				fseq = bin->val;
				index = 0;

				while(getline(&fptr, &n, fp) != -1 && fptr[0] != '>'){
					ptr = strtok(fptr, " ");
					while(ptr != NULL){
						fseq->quality[index++] = atoi(ptr);
						ptr = strtok(NULL, " \n");
					}
				}

                if(fptr[0] == '>') continue;
			}
		}

		if(feof(fp) || ferror(fp) || getline(&fptr, &n, fp) == -1){
			break;
		}
	}

	fclose(fp);
	ckfree(fptr);
	ckfree(prefix);
}

/* iterate through the files and fill in the sequences and quality values */
static void  read_details(const char* const directory,
			   		      hashtable* const names)
{
	DIR *dp;
	struct dirent *ep;
			       
	dp = opendir (directory);
	if (dp != NULL){
		while((ep = readdir (dp)) != NULL){
			if(selector(ep) != 0){
				read_files(directory, ep->d_name, names);
				fprintf(stderr,"Processed %s\n", ep->d_name);
			}
		}
		closedir (dp);
	}else{
		perror ("Couldn't open the directory");
	}
}

/* print the aligments from the assembly */
static void print_alignments(const contig* const contigs, 
							 hashtable* const names)
{
	const contig* c;
	fsequence* r;
	information* info;
	
	int i, j, k, counter;
	uint start;
	bin* iter;
	for(c = contigs; c; c = c->next){
		printf(">Contig\n");
		printf("Num Bases: %d\n", c->num_bases);
		printf("Num Reads: %d\n", c->num_reads);

		info = c->info;
		for(i = 0; i < c->num_bases; i++){
			printf("%c", info[i].consensus);
		}
		printf("\n");

		for(i = 0; i < names->size; i++){
			iter = names->bins[i];
			while(iter){
				for(counter = 0; counter < c->num_reads; counter++){
					r = c->reads[counter];
					if(strncmp(r->person, iter->name, strlen(r->person)) == 0){
						j = r->offset < 0 ? 0 : r->offset;
						for(k = 0; k < j; k++){
							printf(" ");
						}

						if(r->end > r->start){
							start = r->offset < 0 ?
									r->start - r->offset : r->start;
							while(start < r->end){
								printf("%c", r->sequence[start++]);
							}
							printf("\n");
						}else{
							start = r->offset < 0 ? 
									r->end - r->offset : r->end;
							while(start < r->start){
								k = strlen(r->sequence)-start-1;
								printf("%c", tolower(rcs[(int)r->sequence[k]]));
								start++;
							}
							printf("\n");
						}
					}
				}

				iter = iter->next;
			}
		}
	}
}

//static int compare_tle(const void* const elem1, const void* const elem2)
//{
//	tle* a = *(tle**)elem1;
//	tle* b = *(tle**)elem2;
//
//	return a->id - b->id;
//}

/* read the afg output file and return the contig with all the details */
static contig* read_afg_file(FILE* const fp,
						     hashtable* const details)
{
	size_t n = 1;
	char* fptr = ckalloc(n+1);

	hashtable* read_sequences = new_hashtable(8);
	char* identifier = ckalloc(32*sizeof(char));

	contig* contigs = NULL;	/* the return value, set of assembled contigs */
	contig* result  = NULL;	/* a single contig */

	char consensus[4096]; /* the consensus sequence */
	uint iid; 			  /* idnetifier for a read  */
	int i, j, k, pos;
	char* read_sequence;
	fsequence* fseq;
	char* ptr;
	
	information* info;
	allele* variant;
	
	/* just read the details of all the reads used in this assembly */
	while(getline(&fptr, &n, fp) != -1){
		if(strncmp(fptr, "{RED", 4) == 0){
			getline(&fptr, &n, fp);
			assert(strncmp(fptr, "iid:", 4) == 0);
			if(sscanf(fptr, "iid:%d\n", &iid) != 1){
				fatalf("error in reading the iid of the read:%s", fptr);
			}
			sprintf(identifier, "%d", iid);

			/*if this hash more than one assignment, then I should just ignore *
			 * these contigs*/
			if(lookup_hashtable(read_sequences, 
								identifier, 
								strlen(identifier)) != NULL){
				free_hashtable_completely(&read_sequences);
				ckfree(identifier);
				ckfree(fptr);
				return NULL;
			}

			getline(&fptr, &n, fp);
			assert(strncmp(fptr, "eid:", 4) == 0);
			getline(&fptr, &n, fp);
			assert(strncmp(fptr, "seq:", 4) == 0);

			i = 0;
			while(getline(&fptr, &n, fp) != -1 && fptr[0] != '.'){
				memcpy(consensus + i, fptr, strlen(fptr) - 1);
				i += (strlen(fptr) - 1);
			}
			consensus[i] = 0;

			add_hashtable(read_sequences, 
						  identifier, 
						  strlen(identifier),
						  copy_string(consensus));
		}
	}

	/* rewind the file and then find the details of the contigs */
	if(fseek(fp, 0L, SEEK_SET) != 0){
		fatalf("error in resetting stream %s", strerror(errno));
	}

	/* read the contigs */
	while(getline(&fptr, &n, fp) != -1){
		if(strncmp(fptr, "{CTG", 4) == 0){
			if(result != NULL){
				sladdhead(&contigs, result);
			}

			result = ckallocz(sizeof(contig));	
			getline(&fptr, &n, fp);
			assert(strncmp(fptr,"iid:", 4) == 0);
			getline(&fptr, &n, fp);
			assert(strncmp(fptr,"eid:", 4) == 0);
			getline(&fptr, &n, fp);
			assert(strncmp(fptr,"seq:", 4) == 0);
			
			i = 0;
			while(getline(&fptr, &n, fp) != -1 && fptr[0] != '.'){
				memcpy(consensus + i, fptr, strlen(fptr) - 1);
				i += (strlen(fptr) - 1);
			}
			consensus[i] = 0;
			result->num_bases = i;
			result->num_reads = 0;

			/*fill in the consensus*/
			result->info = ckallocz(i*sizeof(information));
			for(i = 0; i < result->num_bases; i++){
				result->info[i].consensus = consensus[i];
			}

			int open = 1;
			while(getline(&fptr, &n, fp) != -1 && open > 0){
				if(fptr[0] == '{') open++;
				if(fptr[0] == '}') open--;

				if(strncmp(fptr,"{TLE", 4) == 0){
					getline(&fptr, &n, fp);
					assert(strncmp(fptr, "src:", 4) == 0);
					read_sequence = must_find_hashtable(read_sequences, 
														fptr + 4, 
														strlen(fptr + 4) - 1);
					fseq = must_find_hashtable(details, 
											   read_sequence,
											   strlen(read_sequence));
					
					getline(&fptr, &n, fp);
					assert(strncmp(fptr, "off:", 4) == 0);
					fseq->offset = atoi(fptr + 4);
	
					getline(&fptr, &n, fp);
					assert(strncmp(fptr, "clr:", 4) == 0);
					if((ptr = strchr(fptr, ',')) == NULL){
						fatalf("error in format of the clear range");
					}
					*ptr = 0;
					fseq->start = atoi(fptr + 4);
					fseq->end = atoi(ptr + 1);

					result->num_reads++;	
					result->reads = ckrealloc(result->reads,
									 result->num_reads*sizeof(fsequence*));
					result->reads[result->num_reads - 1] = fseq;

					/* add all the read information to the contig */
					j = fseq->offset < 0 ? 0 : fseq->offset; //on the contig
					if(fseq->end > fseq->start){
						k = fseq->offset < 0 ? 
							fseq->start - fseq->offset : fseq->start;
				
						while(k < (int)fseq->end && j < result->num_bases){
							variant = ckalloc(sizeof(allele));		
							variant->base = fseq->sequence[k];
							variant->quality = fseq->quality[k];
							variant->read = fseq;
						
							info = result->info + j;
							sladdhead(&info->alleles, variant); 
							j++; k++;
						}
					}else{
						k = fseq->offset < 0 ? 
							fseq->end - fseq->offset : fseq->end;
						while(k < (int)fseq->start && j < result->num_bases){
							variant = ckalloc(sizeof(allele));		
							pos = strlen(fseq->sequence) - k - 1;
							variant->base = rcs[(int)fseq->sequence[pos]];
							variant->quality = fseq->quality[pos];
							variant->read = fseq;
						
							info = result->info + j;
							sladdhead(&info->alleles, variant); 
							j++; k++;
						}
					}
				}
			}
		}
	}
	
	if(result != NULL){
		sladdhead(&contigs, result);
	}

	free_hashtable_completely(&read_sequences);
	ckfree(identifier);
	ckfree(fptr);	

	return contigs;
}


///* read the afg output file and return the contig with all the details */
//static contig* read_afg_file(FILE* const fp,
//						     hashtable* const details)
//{
//	size_t n = 1;
//	char* fptr = ckalloc(n+1);
//
//	contig* contigs = NULL;
//	contig* result = NULL;
//
//	int read_counter;	/*which read in the contig are we handling*/
//	int num_reads;		/*how many reads does this contig have?*/
//
//	char consensus[4096];
//	char* ptr;
//	int i, j, k, pos;
//
//	tle* rds = NULL;
//	tle* rd = NULL;
//
//	fsequence* fseq;
//	uint iid;
//	information* info;
//	allele* variant;
//
//	while(getline(&fptr, &n, fp) != -1){
//		/* read the contig consensus sequence */
//		if(strncmp(fptr,"{CTG", 4) == 0){
//			if(result != NULL){
//				slfreelist(&rds);
//				sladdhead(&contigs, result);
//			}
//
//			result = ckallocz(sizeof(contig));
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr,"iid:", 4) == 0);
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr,"eid:", 4) == 0);
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr,"seq:", 4) == 0);
//			
//			i = 0;
//			while(getline(&fptr, &n, fp) != -1 && fptr[0] != '.'){
//				memcpy(consensus + i, fptr, strlen(fptr) - 1);
//				i += (strlen(fptr) - 1);
//			}
//			consensus[i] = 0;
//			result->num_bases = i;
//
//			result->num_reads = 0;
//			rds = NULL;
//			num_reads = 0;
//
//			/*fill in the consensus*/
//			result->info = ckallocz(i*sizeof(information));
//			for(i = 0; i < result->num_bases; i++){
//				result->info[i].consensus = consensus[i];
//			}
//		}
//
//		/* read the offsets and clear ranges for the reads */
//		if(strncmp(fptr, "{TLE", 4) == 0){
//			i = 0;  //used as a counter for lines
//			num_reads++;
//			rd = ckalloc(sizeof(tle));
//			while(i++ < 3){
//				if(getline(&fptr, &n, fp) == -1){
//					fatalf("error in reading the tle section");
//				}
//				if(strncmp(fptr, "src:", 4) == 0){
//					rd->id = atoi(fptr + 4);
//				}else if(strncmp(fptr, "off:", 4) == 0){
//					rd->offset = atoi(fptr + 4);
//				}else if(strncmp(fptr, "clr:", 4) == 0){
//					if((ptr = strchr(fptr, ',')) == NULL){
//						fatalf("error in format of the clear range");
//					}
//					*ptr = 0;
//					rd->start = atoi(fptr + 4);
//					rd->end = atoi(ptr + 1);
//				}else{
//					fatalf("error in determining the tle section of read");
//				}
//			}
//			sladdhead(&rds, rd);
//		}
//
//		/* read the sequence of the reads */
//		if(strncmp(fptr, "{RED", 4) == 0){
//			if(result && result->num_reads == 0){
//				/* if any read has a duplicate location assigned to it, then
//				 * its just better to ignore this contig*/
//				if(rds != NULL){
//					slsort(&rds, compare_tle);
//					iid = 0;
//					for(rd = rds; rd; rd = rd->next){
//						if(rd->id == iid){
//							slfreelist(&rds);
//							rds = NULL;
//							num_reads = 0;
//							break;
//						}else{
//							iid = rd->id;
//						}
//					}
//				}
//				result->num_reads = num_reads;
//				result->reads = NULL;
//				if(num_reads > 0){
//					result->reads = ckallocz(num_reads*sizeof(fsequence*));
//				}
//				read_counter = 0;
//			}
//
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr, "iid:", 4) == 0);
//			if(sscanf(fptr, "iid:%d\n", &iid) != 1){
//				fatalf("error in reading the iid of the read:%s", fptr);
//			}
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr, "eid:", 4) == 0);
//			getline(&fptr, &n, fp);
//			assert(strncmp(fptr, "seq:", 4) == 0);
//
//			i = 0;
//			while(getline(&fptr, &n, fp) != -1 && fptr[0] != '.'){
//				memcpy(consensus + i, fptr, strlen(fptr) - 1);
//				i += (strlen(fptr) - 1);
//			}
//			consensus[i] = 0;
//
//			fseq = must_find_hashtable(details, consensus, strlen(consensus));
//			if(rds != NULL){
//				/* fill in the details from rds */	
//				for(rd = rds; rd; rd = rd->next){
//					if(rd->id == iid){
//						break;
//					}
//				}
//				if(rd != NULL){
//					fseq->offset = rd->offset;
//					fseq->start = rd->start;
//					fseq->end = rd->end;
//				}else{
//					continue;
//				}
//			}
//			if(result != NULL && result->num_reads > 0){
//				result->reads[read_counter] = fseq;
//				
//				/* add this read information to the contig */
//				j = fseq->offset < 0 ? 0 : fseq->offset; //on the contig
//				if(fseq->end > fseq->start){
//					k = fseq->offset < 0 ? 
//						fseq->start - fseq->offset : fseq->start;
//				
//					while(k < (int)fseq->end && j < result->num_bases){
//						variant = ckalloc(sizeof(allele));		
//						variant->base = fseq->sequence[k];
//						variant->quality = fseq->quality[k];
//						variant->read = fseq;
//					
//						info = result->info + j;
//						sladdhead(&info->alleles, variant); 
//						j++; k++;
//					}
//				}else{
//					k = fseq->offset < 0 ? 
//						fseq->end - fseq->offset : fseq->end;
//					while(k < (int)fseq->start && j < result->num_bases){
//						variant = ckalloc(sizeof(allele));		
//						pos = strlen(fseq->sequence) - k - 1;
//						variant->base = rcs[(int)fseq->sequence[pos]];
//						variant->quality = fseq->quality[pos];
//						variant->read = fseq;
//					
//						info = result->info + j;
//						sladdhead(&info->alleles, variant); 
//						j++; k++;
//					}
//				}
//			}
//			read_counter++;			
//		}
//	}
//	if(result != NULL){
//		slfreelist(&rds);
//		sladdhead(&contigs, result);
//	}
//
//	ckfree(fptr);
//
//	return contigs;
//}

static bool is_interesting(const char* const counts)
{
	bool interesting = FALSE;

	int i;
	for(i = 0; i < 5; i++){
		if(counts[i] > 1){
			if(interesting == TRUE){
				return TRUE;
			}else{
				interesting =  TRUE;
			}
		}
	}

	return FALSE;
}

static void print_alleles(allele* const alleles, 
						  const char* const name, 
						  const char base)
{
    int count = 0;
    allele* iter;
    for(iter = alleles; iter; iter = iter->next){   
        if(iter->base == base &&
           strcmp(iter->read->person, name) == 0){
            count++;
        }
    }

	if(count < 2){
		return;
	}

	printf("%s:%c=%d(", name, base, count);

	int num_processed = 0;
	for(iter = alleles; iter; iter = iter->next){
		if(iter->base == base && 
		   strcmp(iter->read->person, name) == 0){
		   	printf("%d", iter->quality);
		   	num_processed++;
		   	if(num_processed < count){
		   		printf(",");
		   	}else{
				printf(")\t");
				return;
			}
		}
	}
}

/* remember that this is illumina reads. homopolymer stretches are not as big a
 * problem as for 454. Look for positions where:
 * a) there are at least 4 reads.
 * b) there are 2 reads for each allele
 */
static void analyze_contig(contig* const contig, 
						   hashtable* const names)
{
	int i, j;
	information info;
	allele* iter;

	char* counts = ckalloc(5*sizeof(char));
	bin* bin;

	for(i = 0; i < contig->num_bases; i++){
		info = contig->info[i];
		memset(counts, 0, 5*sizeof(char));
		
		for(iter = info.alleles; iter; iter = iter->next){
			switch(iter->base){
				case 'A': counts[0] = counts[0] + 1; break;
				case 'C': counts[1] = counts[1] + 1; break;
				case 'G': counts[2] = counts[2] + 1; break;
				case 'T': counts[3] = counts[3] + 1; break;
				default : counts[4] = counts[4] + 1; break;
			}
		}

		if(is_interesting(counts) && counts[4] < 2){
			printf("%d\t%d\t", i, i+1);
			for(j = 0; j < names->size; j++){
				bin = names->bins[j];

				while(bin){
					print_alleles(info.alleles, bin->name, 'A');
					print_alleles(info.alleles, bin->name, 'C');
					print_alleles(info.alleles, bin->name, 'G');
					print_alleles(info.alleles, bin->name, 'T');
					bin = bin->next;
				}
			}
			printf("\n");
		}
	}

	ckfree(counts);
}


static void free_contig(contig* c)
{
	ckfree(c->reads);
	
	int i;   
	information info;	

	for(i = 0; i < c->num_bases; i++){
		info = c->info[i];
		slfreelist(&info.alleles);
	}
	ckfree(c->info);
}

static void free_contigs(contig** pcontigs)
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

static void print_contig_alleles(const char* const directory,
								 hashtable* const details,
								 hashtable* const names)
{
	static int clusterindex = 0;

	/* this is the file with the output alignments */
	char assemblyfile[1024];
	sprintf(assemblyfile, "%s/scratch/velvet_asm.afg", directory);

	/* read the contig from the assemblyfile */
	struct stat sb;
	if(stat(assemblyfile, &sb) != 0 || sb.st_size == 0){
		return;
	}

	FILE* fp = ckopen(assemblyfile, "r");
	contig* contigs = read_afg_file(fp, details);

	contig* iter;
	information* info;
	int i, j;

	printf(">Cluster_%d\n", clusterindex++);	
	j = 1;
	for(iter = contigs; iter; iter = iter->next){
		printf(">contig%05d\tlength=%d\tnumreads=%d\n", 
			   j++, iter->num_bases, iter->num_reads);
		info = iter->info;
		for(i = 0; i < iter->num_bases; i++){
			printf("%c", info[i].consensus);
		}
		printf("\n");
	}

	if(OnlyAssemble == TRUE){
		goto clean;
	}else{
		for(iter = contigs; iter; iter = iter->next){
			analyze_contig(iter, names);
		}
		if(Alignments){
			print_alignments(contigs, names);
		}
	}
	printf("\n");

clean:
	fclose(fp);
	free_contigs(&contigs);
}

/* build the assemblies for each cluster and process the output to find the
 * candidate SNPs*/
void build_cluster_assemblies(const char* const clusters_file, 
							  hashtable* const reads)
{
	char template[10] = "JOBXXXXXX";
	char* dir;
	if((dir = mkdtemp(template)) == NULL){
		fatalf("error in creating temp directory called %s : %s", 
				template, strerror(errno));
	}

	char readsfile[1024];
	sprintf(readsfile,"%s/reads.fq", dir);

	size_t n = 1;
	char* fptr = ckalloc(n+1);
	FILE* fp = ckopen(clusters_file, "r");
	FILE* cp = NULL;

	if(getline(&fptr, &n, fp) == -1){
		fatalf("error in reading the cluster from %s", clusters_file);
	}

	/* this hashtable will keep the names of all the different individuals */
	hashtable* names = new_hashtable(3);

	/* this hashtable stores the reads for the current cluster */
	hashtable* mapping = NULL;

	char header[1024];
	char name[1024];
	fsequence* fseq;
	uint i;
	bin* bin;
	while(1){
		if(fptr[0] == '>'){
			cp = ckopen(readsfile, "w");
			mapping = new_hashtable(4);

			while(getline(&fptr,&n,fp) != -1 &&
				  sscanf(fptr, "%s %s\n", header, name) == 2){
				fseq = must_find_hashtable(reads, header, strlen(header));
				fprintf(cp, "@%s\n", header);
				fprintf(cp, "%s\n", fseq->sequence);
				fprintf(cp, "+\n");
				for(i = 0; i < strlen(fseq->sequence); i++){
					fprintf(cp, "%c", fseq->quality[i] + 33);
				}
				fprintf(cp, "\n");

				/*if this individual has not been added to the list, do that*/
				if((bin = lookup_hashtable(names, name, strlen(name)))== NULL){
					bin = add_hashtable(names, name, strlen(name), NULL);
				}
				fseq->person = bin->name;

				/* since velvet removes the names of the reads, and replaces
				 * them with identifiers, we need to construct a mapping
				 * between the bases of a sequence, and the (quals, person) for
				 * that read.*/
				for(i = 0; i < strlen(fseq->sequence); i++){
					if(fseq->sequence[i] == 'N'){
						fseq->sequence[i] = 'A';
					}
				}
				add_hashtable(mapping,
							  fseq->sequence,
							  strlen(fseq->sequence),
							  fseq); 
			}
			fclose(cp);

			sprintf(command, "rm -rf %s/scratch", dir);
			run_command(command);

			/*run velvet on the file with the reads*/
			sprintf(command, 
					"%s/velveth %s/scratch 21 -fastq %s/reads.fq > /dev/null", 
					 rig_directory,dir, dir);
			run_command(command);

			sprintf(command,
			  "%s/velvetg %s/scratch -read_trkg yes -amos_file yes > /dev/null",
	  	      rig_directory,dir);
			run_command(command);

			/*lets display the contigs and parse the alignments*/
			print_contig_alleles(dir, mapping, names);

			/* free all the stuff from the hashtable for this cluster */
			free_hashtable(&mapping);
		}

		if(feof(fp) || ferror(fp) || getline(&fptr, &n, fp) == -1){
			break;
		}
	}

	fclose(fp);
	ckfree(fptr);
	free_hashtable(&names);

	sprintf(command, "rm -rf %s", dir);
	run_command(command);
}

static void free_hash(bin* bin)
{
	ckfree(bin->name);
	fsequence* fseq = bin->val;
	if(fseq != NULL){
		ckfree(fseq->sequence);
		ckfree(fseq->quality);
		ckfree(fseq);
	}
	ckfree((char*)bin);
}

static void assemble_illumina(const char* const reads_directory, 
						      const char* const clusters_file)
{
	/* read the names of the reads to be fetched from the clusters file */
	hashtable* names = read_names(clusters_file);

	/* fill in the hashtable with all the read values */
	read_details(reads_directory, names);

	/* process the clusters one by one, printing out the reads and processing
	 * them with Velvet. The output of velvet would then be processed to 
	 * find the candidate SNPs*/
	build_cluster_assemblies(clusters_file, names);

	/* free all the resources */
	func_hashtable(names, free_hash);
	ckfree(names->bins);
	ckfree(names);
}

int main(int argc, char** argv)
{
	argv0="assemble_illumina";
	/*read the command line arguments*/
	while(argc > 3){
		argc--;
		if(strncmp(argv[argc],"-alignments",11) == 0){
			Alignments = TRUE;	
		}else if(strncmp(argv[argc], "--velvet=", 9) == 0){
			rig_directory = argv[argc] + 9;
		}else if(strncmp(argv[argc],"-onlyassemble", 13) == 0){
			OnlyAssemble = TRUE;	
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

	command = ckalloc(DSIZE*sizeof(char));

	assemble_illumina(argv[1], argv[2]);

	ckfree(command);
	return EXIT_SUCCESS;
}
