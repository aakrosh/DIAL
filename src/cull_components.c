/* This module takes in a clusters.txt file and generates an input file which
 * can be used by the 'assemble' module. */

#include "utilities.h"
#include "hashtable.h"
#include "graph.h"
#include "cluster.h"

#define USE "cull_components clusters.txt [-strict]"

bool Strict = FALSE;

bool Short = FALSE;

typedef struct clique_st
{
	struct clique_st* next;
	int size;
	node** pnodes;
}clique;

static bool mergeable(const clique* const c1, const clique* const c2)
{
	int i, j;

	for(i = 0; i < c1->size; i++){
		for(j = 0; j < c2->size; j++){
			if(c1->pnodes[i] != c2->pnodes[j] && 
			   !edge_exists(c1->pnodes[i], c2->pnodes[j])){
				return FALSE;
			}
		}
	}

	return TRUE;
}

static clique* find_minimal(node* const start,
							node* const end)
{
	node* iter;
	clique* list = NULL;
	clique* part;
	for(iter = start; iter != end->next; iter = iter->next){
		part = ckalloc(sizeof(clique));
		part->size = 1;
		part->pnodes = ckalloc(sizeof(node*));
		part->pnodes[0] = iter;
		sladdhead(&list, part);
	}
	return list;
}

static void merge(clique* const p, clique* const q)
{
	int i, j;
	bool same;
	for(i = 0; i < q->size; i++){
		same = FALSE;
		for(j = 0; j < p->size; j++){
			if(p->pnodes[j] == q->pnodes[i]){
				same = TRUE;
				break;
			}
		}
		if(!same){
			p->pnodes = ckrealloc(p->pnodes, (p->size+1) * sizeof(node*));
			p->size++;
			p->pnodes[p->size - 1] = q->pnodes[i];
		}
	}
}

/* the component starts is inclusinve of start and end nodes. Find all the 
 * cliques in this component.*/
static clique* find_cliques(node* const start, 
						    node* const end)
{
	clique* list = find_minimal(start, end);

	clique* iter1 = list->next;
	clique*  prev = list;
	clique* iter2;
	clique* iter3;

	while(iter1){
		iter2 = list; 
		while(iter2 != iter1){
			if(mergeable(iter1, iter2)){
				merge(iter2, iter1);
				iter3 = slremove(&list, iter1);
				ckfree(iter3->pnodes);
				ckfree(iter3);
				iter1 = prev;
				break;
			}else{
				iter2 = iter2->next;
			}
		}
		prev = iter1;
		iter1 = iter1->next;
	}

	return list;
}

/* if the number of nodes in a component is >= the threshold, then find the
 * cliques in that component */
static void process_components(graph* const network, const int threshold)
{
	slsort(&network->nodelist, component_sort);
	
	node* iter;
	node* start;
	int count;

	clique* cliques;
	clique* citer;

	int i, id = 0;

	int num_components = 0;

	for(iter = network->nodelist; iter; iter = iter->next){
		start = iter;
		count = 1;

		while(start && start->next && 
			  start->next->component == iter->component){
			count++;
			start = start->next;
		}

		if(count >= threshold){	
			num_components++;
			cliques = find_cliques(iter, start);
			for(citer = cliques; citer; citer = citer->next){
				if(citer->size > threshold){
					printf(">Cluster_%d\n", id++);
					for(i = 0; i < citer->size; i++){
						printf("%s\n", citer->pnodes[i]->name);
					}
					printf("\n");
				}
				ckfree(citer->pnodes);
			}
			slfreelist(&cliques);
		}

		iter = start;
	}
	fprintf(stderr,"Number of components which pass the threshold: %d\n", 
			num_components);
}

static void clique_components(const char* const clusterfile)
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

	/* find all the connected components in the graph */
	int components = find_connected_components(network);
	fprintf(stderr,"Found %d connected components in the graph\n", components);

	/* find all the clique with k=4 and above */
	process_components(network, 4);

	free_graph(&network);
	free_hashtable(&names);
}

static void simple_prints(const char* const clusterfile)
{
	cluster* cp;
	if((cp = read_cluster_file(clusterfile)) == NULL){
		fatalf("error in reading the first clutser from %s", clusterfile);
	}
	
	hashtable* names = new_hashtable(16);
	bool print;
	uint i, j = 0;

	while(cp){
		print = TRUE;
		for(i = 0; i < cp->num_members; i++){
			if(lookup_hashtable(names, 
								cp->members[i], 
								strlen(cp->members[i])) != NULL){
				print = FALSE;
				break;
			}else{
				add_hashtable(names, 
							  cp->members[i], 
							  strlen(cp->members[i]), 
							  NULL);
			}
		}
		
		if(print == TRUE && cp->num_members >= 4){
			printf(">Cluster_%d\n", j++);
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
	fprintf(stderr,"Read all the clusters from %s\n", clusterfile);
	free_hashtable(&names);
}

/* if the number of nodes in a component is >= the minthreshold, and <=
 * maxthreshold, then print the component out to stdout */
static void print_components(graph* const network, const int min, const int max)
{
	slsort(&network->nodelist, component_sort);
	
	node* iter;
	node* start;
	node* iter2;
	int count;
	
	int id = 0;

	int num_components = 0;

	for(iter = network->nodelist; iter; iter = iter->next){
		start = iter;
		count = 1;

		while(start && start->next && 
			  start->next->component == iter->component){
			count++;
			start = start->next;
		}

		if(count >= min && count <= max){	
			num_components++;
			printf(">Cluster_%d\n", id++);
				for(iter2 = iter; iter2 != start->next; iter2 = iter2->next){
					printf("%s\n", iter2->name);
				}					
			printf("\n");
		}
		iter = start;
	}
	fprintf(stderr,"Number of components which pass the threshold: %d\n", 
			num_components);

}

/* Try to find the connected components for all the clusters in the input file.
 * Then only choose the clusters that have the size less than the prescribed.
 * Right now I am going to hard code that size....*/
static void connect_components(const char* const clusterfile)
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

	/* find all the connected components in the graph */
	int components = find_connected_components(network);
	fprintf(stderr,"Found %d connected components in the graph\n", components);

	/* find all the clique with k=4 and above */
	print_components(network, 4, 6);

	free_graph(&network);
	free_hashtable(&names);
}

static void cull_components(const char* const clusterfile)
{
	void (*func)(const char* const) = NULL;

	if(Short == TRUE){
		func = connect_components;
	}else if(Strict == TRUE){
		func = clique_components;
	}else{
		func = simple_prints;
	}

	func(clusterfile);
}

int main(int argc, char** argv)
{
	argv0="cull_components";
	while(argc > 2){
		argc--;
		if(strncmp(argv[argc],"-strict", 7) == 0){
			Strict = TRUE;
		}else if(strncmp(argv[argc], "-short", 6) == 0){
			Short = TRUE;
		}else{
			fatalf("unknown option : %s\n%s", argv[argc], USE);
		}
	}

	if(argc != 2){
		fatal(USE);
	}

	allocate_resources();

	cull_components(argv[1]);

	return EXIT_SUCCESS;
}
