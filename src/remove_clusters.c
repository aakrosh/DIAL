/* Simple module, which reads in the lastz output and if the cluster aligns 
 * with any of the reads with an overlap of more than 100 base pairs, then 
 * it doesnt print the cluster */

#include "mfblocks.h"
#include "cluster.h"

#define USE "remove_clusters clusters.txt < lastz.out"

uint MinOlap = 100;

static void remove_clusters(const char* const clusterfile)
{
	cluster* cp;
	if((cp = read_cluster_file(clusterfile)) == NULL){
		fatalf("unable to read the clusters from %s", clusterfile);
	}

	mafblock* mp;
	if((mp = read_maf_file("/dev/stdin")) == NULL){
		fatalf("error in reading the maf file");
	}

	mfblocks* blocks;
	mfblocks* mfnode;
	uint i, removed = 0, total = 0;

	while(cp){
		blocks = NULL;
		total++;

		/* is the current maf block for this read? */
		if(strncmp(cp->name, (char*)mp->name2, strlen((char*)mp->name2)) == 0){
			/* fetch all the remaining maf blocks for this read */
			do{
				if(mp->o1 > MinOlap && mp->o2 > MinOlap){
					mfnode = new_mfnode(mp);
					sladdhead(&blocks, mfnode);		
				}
			}while(mp && 
			 get_next_block(mp) && 
			 strncmp(cp->name,(char*)mp->name2, strlen((char*)mp->name2)) == 0);
		}
		
		if(blocks != NULL){
			removed++;
			free_mfblocks(&blocks);
		}else{
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
	fprintf(stderr,"%d/%d clusters removed due to similarity to repeats\n",
					removed, total);
}

int main(int argc, char** argv)
{
	while(argc > 2){
		argc--;
		if(strncmp(argv[argc],"--minlen=", 9) == 0){
			MinOlap = atoi(argv[argc] + 9);
		}
	}	

	if(argc != 2){
		fatal(USE);
	}

	remove_clusters(argv[1]);

	return EXIT_SUCCESS;
}
