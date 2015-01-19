/* A module which takes several files of clusters as input. The first argument
 * is the file of clusters that have passed the "filter_clusters" module. The
 * remaining arguments are files of clusters that this module fishes into to
 * find and "grow" the cluster. Then they are output to the screen */

#include "utilities.h"
#include "hashtable.h"
#include "graph.h"
#include "cluster.h"

#define USE "fish_clusters clusters.txt cluster1.txt ..."

/* read the cluster file as a graph */
static graph* read_cluster_graph(const char* const clusterfile)
{
	cluster* cp;
	if((cp = read_cluster_file(clusterfile)) == NULL){
		fatalf("error in reading the first clutser from %s", clusterfile);
	}

	uint i, singles = 0;

	graph* network = new_graph();			/* the graph of reads,individuals */
	hashtable* names = new_hashtable(3);	/* hashtable to store the names of
											 * ... the different individuals */

	char name[1024];
	char read[1024];

	while(cp){
		/* add the first member. NOTE that this member owns the sequence  */
		if(sscanf(cp->members[0],"%s %s", read, name) != 2){
			fatalf("error in reading info about cluster:%s", cp->name);
		}

		/* is the name of the individual in the hash */
		if(lookup_hashtable(names, name, strlen(name)) == NULL){
			add_hashtable(names, name, strlen(name), NULL);
		}

		if(cp->sequence[0] != '?'){
			/* add this node to the graph */
			add_node(network, cp->members[0], NULL);

			for(i = 1; i < cp->num_members; i++){
				if(sscanf(cp->members[i],"%s %s", read, name) != 2){
					fatalf("error in reading info about cluster:%s", cp->name);
				}

				/* is the name of the individual in the hash */
				if(lookup_hashtable(names, name, strlen(name)) == NULL){
					add_hashtable(names, name, strlen(name), NULL);
				}

				/* add this node to the graph */
				add_node(network, cp->members[i], NULL);

				/* make the edges for the graph */
				make_edges(network, cp->members[0], cp->members[i], 0);
			}
		
			if(cp->num_members == 1){
				singles++;
			}
		}
		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);
	fprintf(stderr,"Read all the clusters from %s\n", clusterfile);
	return network;
}

#if 0
static void do_file(graph* const network, const char* const clusterfile)
{
	cluster* cp;
	if((cp = read_cluster_file(clusterfile)) == NULL){
		fatalf("error in reading the first clutser from %s", clusterfile);
	}

	uint i;
	bool to_add;
	bin* bin;

	while(cp){
		to_add = FALSE;
		for(i = 0; i < cp->num_members; i++){
			if((bin = lookup_hashtable(network->nodes, 
									   cp->members[i],
  								   	   strlen(cp->members[i]))) != NULL){
				to_add = TRUE;
				break;
			}
		}	

		if(to_add == TRUE){
			for(i = 0; i < cp->num_members; i++){
				add_node(network, cp->members[i], NULL);
				make_edges(network, bin->name, cp->members[i], 0);
			}
		}
	
		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);
	fprintf(stderr,"Read all the clusters from %s\n", clusterfile);
}
#endif

static void do_file(graph* const network, const char* const clusterfile)
{
	cluster* cp;
	if((cp = read_cluster_file(clusterfile)) == NULL){
		fatalf("error in reading the first clutser from %s", clusterfile);
	}

	bin* bin;

	while(cp){
		if((bin = lookup_hashtable(network->nodes,
								   cp->members[0],
								   strlen(cp->members[0]))) != NULL){
			/*if this cluster is marked as "probable repeat" then we should 
			 *just ignore this cluster*/
			if(cp->sequence[0] == '?'){
				remove_all_edges(network, bin->val);
			}
		}
	
		if(!read_next_cluster(cp)){
			break;
		}
	}
	close_cluster_file(cp);
	fprintf(stderr,"Read all the clusters from %s\n", clusterfile);
}

static void print_clusters(graph* const network)
{
	/* find all the connected components in the graph */
	int components = find_connected_components(network);
	fprintf(stderr,"Found %d connected components in the graph\n", components);
	slsort(&network->nodelist, component_sort);
	
	node* iter;
	node* start;
	node* iter2;
	int count;
	
	char buffer[1024];
	char* ptr;

	for(iter = network->nodelist; iter; iter = iter->next){
		start = iter;
		count = 1;

		while(start && start->next && 
			  start->next->component == iter->component){
			count++;
			start = start->next;
		}

		if(count >= 4){
			sprintf(buffer, "%s", iter->name);
			if((ptr = strchr(buffer, '\t')) == NULL){
				fatalf("no tab in the cluster name %s", buffer);
			}
			*ptr = 0;
			
			printf(">%s\n", buffer);
			printf("GAACTAAAAGCAATAAACCTAAACAGAGGTGCTTCATTCTGCAGGAAGCCTGGGGACTGTCCTTTCTTTGTTCAA\n");
			for(iter2 = iter; iter2 != start->next; iter2 = iter2->next){
				printf("%s\n", iter2->name);
			}					
			printf("\n");
		}
		iter = start;
	}
}

int main(int argc, char** argv)
{
	argv0="fish_clusters";

	allocate_resources();

	/* read the clusters into a graph*/
	graph* network = read_cluster_graph(argv[1]);

	int i;
	for(i = 2; i < argc; i++){
		do_file(network, argv[i]);
	}
	fprintf(stderr,"%d nodes in the graph\n", slcount(network->nodelist));

	print_clusters(network);

	return EXIT_SUCCESS;
}


