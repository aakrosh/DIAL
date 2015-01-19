/* Simple module which goes through all the clusters in question and find the
 * ones with at least one difference block */

#include "hashtable.h"
#include "cluster.h"

hashtable* individuals = NULL;

/*parse the comma separated list of names of the individuals and fill them 
 *in the hashtable "individuals" */
static void parse_names(char* const names)
{
	individuals = new_hashtable(3);

	char* ptr = strtok(names, ", ");
	while(ptr != NULL){
		if(lookup_hashtable(individuals,ptr,strlen(ptr)) == NULL){
			add_hashtable(individuals,ptr,strlen(ptr),NULL);
		}		

		ptr = strtok(NULL,", ");
	}
}

static void initialise(bin* bin)
{
	bin->val = 0 + NULL;
}

/*increase the count for this name in the hashtable "individuals"*/
static inline void update_count(char* const name)
{
	bin* bin = lookup_hashtable(individuals, name, strlen(name));
	assert(bin != NULL);

	bin->val += 1;
}

/*iterate through the names in the hashtable "individuals" and return TRUE if
 *at least 2 individuals have 2 reads each. */
static inline bool check_twonames()
{
	int i;
	bin* iter;
	bin* next;
	
	bool first = FALSE;

	for(i = 0; i < individuals->size; i++){
		iter = individuals->bins[i];
		while(iter){
			next = iter->next;
	
			if(iter->val >= (2 + NULL)){
				if(first == FALSE){
					first = TRUE;
				}else{
					return TRUE;
				}
			}		

			iter = next;
		}
	}

	return FALSE;
}

static bool goodblock(const cluster* const cp)
{
	/* how many substitutions do we have? This is not a good window cluster if
	 * we have more than one substitution difference in a window of 40 bp */
	uint i, j, k = -40;
	bool good = TRUE;
	for(i = 0, j = 0; i < cp->num_diffs; i++){
		if(cp->differences[i].reference != '-' &&
		   cp->differences[i].target != '-'    && 
		   cp->differences[i].reference < 'a'  && 
		   cp->differences[i].target < 'a' && 
		   cp->differences[i].reference != 'N' &&
		   cp->differences[i].target != 'N'){
			j++;
			if((cp->differences[i].position - k) < 40){
				good = FALSE;
			}else{
				k = cp->differences[i].position;
			}
		}
	}

	/*how many individuals do we have?*/
	int num_people = individuals->elcount;

	/* if it is a single person, then the conditions vary from that if there 
	 * are more than one */
	if(num_people == 1){
		if(j >= 1 && j < 3 && cp->num_members >= 4 && good == TRUE){
			for(i = 0; i < cp->num_diffs; i++){
				if(cp->differences[i].reference != '-' &&
				   cp->differences[i].target != '-' &&
				   cp->differences[i].reads > 1){
					return TRUE;
				}
			}
		}
	}else{
		char name[1024];

		/* initialise the hash counts for all the individuals to 0 */
		func_hashtable(individuals, initialise);

		for(i = 0; i < cp->num_members; i++){
			if(sscanf(cp->members[i],"%*s %s", name) != 1){
				fatalf("error in reading the member %s", cp->members[i]);
			}
			update_count(name);
		}

		/* do we have at least two reads each for at least two individuals */
		bool twonames = check_twonames();
		if(twonames && j >= 1 && j < 3 && cp->num_members >= 4 && good == TRUE){
			for(i = 0; i < cp->num_diffs; i++){
				if(cp->differences[i].reference != '-' &&
				   cp->differences[i].target != '-' &&
				   cp->differences[i].reads > 1){
					return TRUE;
				}
			}
		}
	}

	return FALSE;
}

/* read this file and print out the "good" clusters */
static void do_file(const char* const file)
{
	cluster* cp;
	if((cp = read_cluster_file(file)) == NULL){
		fatalf("error in reading the clusters from %s", file);
	}

	uint i;
	while(cp){
		if(cp->num_diffs > 0 && goodblock(cp) == TRUE){
			printf(">%s\n", cp->name);
			printf("%s\n", cp->sequence);
			for(i = 0; i < cp->num_members; i++){
				printf("%s\n", cp->members[i]);
			}
			printf("\n");
		}

		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);
}

int main(int argc, char** argv)
{
	if(strncmp(argv[argc-1], "--names=", 8) == 0){
		argc--;
		parse_names(argv[argc] + 8);
	}

	while(argc > 1){
		argc--;
		do_file(argv[argc]);
	}
	
	if(individuals != NULL){
		free_hashtable(&individuals);
	}

	return EXIT_SUCCESS;
}
