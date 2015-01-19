/* This module takes an input an masked fasta file of reads and an length
 * threshold. If the unmasked region of a read is below that length threshold
 * then throw the read away. Another input is the quality file for the same set
 * of reads. The output is the quality scores of the unmaksed bases in the same
 * order */

#include <ctype.h>
#include "common.h"
#include "utilities.h"
#include "sequences.h"
#include "hashtable.h"

#define USE "update reads.fa reads.qual --minlen=<int> --prefix=454"

uint MinLength = 100;
char* Prefix = "454";
char* Repeatsfile = NULL;

/* Examine the read. find the first unmasked region of the read. Mark the
 * longest unmaked region of the read. If it exceeds the length threshold then
 * print the read */
static bool process_read(sequence* const sp, 
						 const int* const quals, 
						 FILE* const sf, 
						 FILE* const qf,
						 FILE* const rf)
{
	if(sp->slen < MinLength){
		return FALSE;
	}

	uint i = 0, j = 0, k = 0, len = 0;

	while(i < sp->slen){
		/* find an masked segment */
		j = i;
		while(j < sp->slen && sp->sequence[j] < 'a'){
			j++;
		}
		
		k = j;
		len = 0;
		while(j < sp->slen && sp->sequence[j] >= 'a'){
			j++;
			len++;
		}

		/* unmask this segment if it is less than 40 bp long */
		if(len < 40){
			for(i = k; i < j; i++){
				sp->sequence[i] = toupper(sp->sequence[i]);
			}
		}

		i = j;
	}

	/* is the length of the unmasked segment less than MinLength ? If so, do not
	 * print this read */
	for(i = 0, j = 0; i < sp->slen; i++){
		if(sp->sequence[i] < 'a'){
			j++;
		}
	}
	if(j < MinLength){
		if(rf){
			fprintf(rf, ">%s\n", sp->header);
			for(i = 0; i < sp->slen; i++){
				fprintf(rf, "%c", toupper(sp->sequence[i]));
			}	
			fprintf(rf, "\n");
		}
		return FALSE;
	}

	fprintf(sf, ">%s\n", sp->header);
    for(i = 0; i < sp->slen; i++){
		fprintf(sf,"%c", sp->sequence[i]);
    }
    fprintf(sf,"\n");
	fprintf(qf, ">%s\n", sp->header);
	for(i = 0; i < sp->slen; i++){
		fprintf(qf,"%d ", quals[i]+1);
	}
	fprintf(qf,"\n");

	return TRUE;
}

/* fill in the quality values for the next read */
static int get_quality_values(const char* const name,
							  char** const fptr, 
							  size_t* n, 
							  FILE* const fp, 
							  int** const quals)
{
	assert((*fptr)[0] == '>');
	int* q = *quals;
	
	/*get to the correct read*/
	while(strncmp(name, *fptr + 1, strlen(name)) != 0){
		if(getline(fptr, n, fp) == -1){
			fatalf("error in reading quality values for %s", name);
		}
	}
	
	char* ptr;
	int i = 0, j = 0;
	while(getline(fptr, n, fp) != -1 && (*fptr)[0] != '>'){
		ptr = strtok(*fptr, " \n");
		j++;
		while(ptr != NULL){
			q[i++] = atoi(ptr);
			ptr = strtok(NULL, " \n");
		}
	}
	return i;
}

static void update(const char* const seqfile, const char* const qualfile)
{
	/* the sequences */
	sequence* sp;
    if((sp = read_fasta_sequence(seqfile)) == NULL){
        fatalf("error in reading the sequences from %s", seqfile);
    }
	
	/* the quality file */
	int* quals = ckalloc(MAX_READ_SIZE*sizeof(int));
	size_t n = 1;
	char* fptr = ckalloc(n+1);
	FILE* fp = ckopen(qualfile,"r");
	while(getline(&fptr, &n, fp) != -1 && fptr[0] != '>');

	/* open the two output files */
	char buffer[128];
	sprintf(buffer,"%s.fa", Prefix);
	FILE* sf = ckopen(buffer,"w");
	sprintf(buffer,"%s.qual", Prefix);
	FILE* qf = ckopen(buffer,"w");
	
	/*do we write the repeats into a special file*/
	FILE* rf = NULL;
	if(Repeatsfile != NULL){
		rf = ckopen(Repeatsfile, "a");
	}
	
	int readsprinted = 0;
	int total = 0;
	uint values;

	while(sp){
		total++;
		if(sp->slen > MAX_READ_SIZE){
			fatalf("I did not expect such huge reads: length %d", sp->slen);
		}
		
		/* get the corresponding quality values from the qualfile */
		values = get_quality_values((char*)sp->header, &fptr, &n, fp, &quals);
		assert(values == sp->slen);

        /* process this read */
	    if(process_read(sp, quals, sf, qf, rf) == TRUE){
			readsprinted++;
		}

        if(!get_next_sequence(sp)){
            break;
        }
    }
    close_fasta_sequence(sp);

	/* print some stats */
	fprintf(stderr, "%d/%d (%3.2f%%) reads with unmasked length less than %d bp\n", total-readsprinted,total,(total-readsprinted)*100.0/total, MinLength);
	
	fclose(sf);
	fclose(qf);
	if(rf){
		fclose(rf);
	}
	fclose(fp);
	ckfree(quals);
}

int main(int argc, char** argv)
{
	argv0="update";

    while(argc > 3){
        argc--;
        if(strncmp(argv[argc],"--minlen=", 9) == 0){
            MinLength = atoi(argv[argc]+9);
		}else if(strncmp(argv[argc],"--prefix=",9) == 0){
			Prefix = argv[argc] + 9;
		}else if(strncmp(argv[argc],"--repeats=", 10) == 0){
			Repeatsfile = argv[argc] + 10;
        }else{
            fprintf(stderr,"unknown option:%s\n", argv[argc]);
            fatal(USE);
        }
    }

    if(argc != 3){
        fatal(USE);
    }
	if(strlen(Prefix) > 120){
		fatalf("Please select a --prefix values less than 120 characters long");
	}
    
	update(argv[1], argv[2]);

    return EXIT_SUCCESS;
}
