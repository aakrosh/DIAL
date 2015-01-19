#include "graph.h"

/* return a new graph */
graph* new_graph()
{
	graph* network = ckallocz(sizeof(graph));

	network->nodes = new_hashtable(16);
	return network;
}

/* add a new node to the graph */
node* add_node(graph* const graph, 
	   	       const char* const name,
			   const uchar* const seq)
{
	node* node = ckalloc(sizeof(struct node_st));
	node->condemned = FALSE;
	node->visited = FALSE;
	node->component = -1;
	node->num_alloced = 0;
	node->num_edges = 0;
	node->edges = NULL;
	node->sequence = NULL;
	if(seq != NULL){
		node->sequence = (uchar*)copy_string((char*)seq);
	}

	bin* bin;
	if((bin = lookup_hashtable(graph->nodes, name, strlen(name))) == NULL){
		bin = add_hashtable(graph->nodes, name, strlen(name), node);
		node->name = bin->name;
		sladdhead(&graph->nodelist, node);
	}

	node = bin->val;
	return node;
}

/* add the name of the nodes to the graph. Return the number of nodes added */
int add_nodes(graph* const graph, const char* const seqfile)
{
	sequence* sp;
	int num_nodes = 0;
	node* node;

	if((sp = read_fasta_sequence(seqfile)) == NULL){
		fatalf("error in reading the sequences from %s", seqfile);
	}

	char* ptr;
	while(sp){
		ptr = strchr((char*)sp->header, ' ');
		if(ptr != NULL){
			*ptr = 0;
		}
		node = add_node(graph, (char*)sp->header, sp->sequence);
		num_nodes++;

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);
	return num_nodes;
}

/* does this edge already exist*/
bool edge_exists(const node* const n1, const node* const n2)
{
	int i;
	edge* edge;

	for(i = 0; i < n1->num_edges; i++){
		edge = n1->edges[i];
		if(edge->a == n2 || edge->b == n2){
			return TRUE;
		}
	}
	return FALSE;
}

/* create  a new edge in this universe :) */
static edge* new_edge(graph* const graph,
					  node* const n1, 
					  node* const n2, 
					  const int weight)
{
	edge* edge = ckalloc(sizeof(struct edge_st));

	edge->a      = n1;
	edge->b      = n2;
	edge->weight = weight;
	edge->condemned = FALSE;
	sladdhead(&graph->edgelist, edge);
	return edge;
}

static void add_new_edge(node* const node, edge* const edge)
{
	/*do i need more memory ? */
	if(node->num_edges == node->num_alloced){
		int diff = node->num_alloced;
		int old_size = node->num_alloced * sizeof(struct edge_st*);
		node->num_alloced +=  NUM_EDGES;
		int new_size = node->num_alloced * sizeof(struct edge_st*);
		node->edges = ckrealloc(node->edges, new_size);
		memset(node->edges+diff, 0, (new_size-old_size));
	}

    node->edges[node->num_edges++] = edge;
}

/* make the edge between the two nodes for this graph */
edge* make_edges(graph* const graph, 
				 const char* const name1,
				 const char* const  name2,
				 const int weight)
{
	node* n1 = find_value_hashtable(graph->nodes, name1, strlen(name1));
	node* n2 = find_value_hashtable(graph->nodes, name2, strlen(name2));
	edge* edge = NULL;

	if(!edge_exists(n1, n2)){
		/* create a new edge and add it to the egde list of the two nodes */
		edge = new_edge(graph, n1, n2, weight);
		add_new_edge(n1, edge);
		add_new_edge(n2, edge);
	}
	return edge;
}


/* free the resources for this graph */
void free_graph(graph** pgraph)
{
	if(*pgraph){
		graph* graph = *pgraph;
		slfreelist(&graph->edgelist);
		
		node* iter;
		for(iter = graph->nodelist; iter; iter = iter->next){
			ckfree(iter->edges);
			if(iter->sequence != NULL){
				ckfree(iter->sequence);
			}
		}
		slfreelist(&graph->nodelist);
		free_hashtable(&graph->nodes);
		ckfree(graph);
		*pgraph = NULL;
	}
}

/* clear all the visit flags for all the nodes */
static void clear_visit_flags(graph* const graph)
{
	node* iter;
	for(iter = graph->nodelist; iter; iter = iter->next){
		iter->visited = FALSE;
	}
}

static void depthfirst(node* const root, const int identifier)
{
	int i;
	node* next;
	root->visited = TRUE;
	assert(root->component == -1);
	root->component = identifier;

	for(i = 0; i < root->num_edges; i++){
		if(root->edges[i]->condemned == TRUE){
			continue;
		}
		next = root->edges[i]->a == root 
 								  ? root->edges[i]->b : root->edges[i]->a;

		if(next->visited == FALSE){
			depthfirst(next, identifier);	
		}
	}
}

/* find and mark the connected comonents in this graph. Return the number of
 * connected components in the graph */
int find_connected_components(graph* const graph)
{
	clear_visit_flags(graph);	

	int componentId = 0;

	node* iter;
	for(iter = graph->nodelist; iter; iter = iter->next){
		if(iter->visited == FALSE){
			/*mark all the nodes that you can reach from this node*/
			depthfirst(iter, componentId);

			componentId++;
		}
	}

	return componentId;
}

/* remove all the edges connected to this node in this graph */
void remove_all_edges(graph* const graph __attribute__((unused)), 
					  node* const node)
{
	int i;
	for(i = 0; i < node->num_edges; i++){
		node->edges[i]->condemned = TRUE;
	}
}

/* go through each of the nodes. If a node has more than 'maxhits edges', then
 * all those edges should be deleted. The node should be marked as a repeat.
 * Return the number of nodes filtered. Write those nodes in the repeatfile */
int filter_nodes(graph* const graph, FILE* const fp, const int maxhits)
{
	int num = 0;
	
	node* iter;
	for(iter = graph->nodelist; iter; iter = iter->next){
		if(iter->num_edges > maxhits){
			remove_all_edges(graph, iter);
			iter->condemned = TRUE;
			num++;
			if(iter->sequence != NULL){
				fprintf(fp,">%s\n", iter->name);
				fprintf(fp,"%s\n",  iter->sequence);
			}
		}
	}
	return num;
}

/* sort the nodes based on the component they are in */
int component_sort(const void* const a, const void* const b)
{
	node* n1 = *((node**)a);
	node* n2 = *((node**)b);

	return n1->component - n2->component;
}

/* print the clusters from this graph */
void print_graph_clusters(graph* const graph, 
						  const char* const clusterfile,
					      const int index,
					      const char* const individual,
						  const int maxhits)
{
	node* iter;
	node* next;
	int i, id = index + 1;

	FILE* fp = ckopen(clusterfile,"aw");

	for(iter = graph->nodelist; iter; iter = iter->next){
		if(iter->sequence == NULL){
			continue;
		}
		
		/* print just the sequence for the next iteration */
		fprintf(fp,">Cluster_%d\n", id);
		fprintf(fp,"%s\n", iter->sequence);
	
		/* print the cluster gory details */
		printf(">Cluster_%d\n", id++);
		if(iter->num_edges > maxhits){
			printf("??");
		}
		printf("%s\n", iter->sequence);
		printf("%s %s\n", iter->name, individual);

		if(iter->num_edges > maxhits){
			iter->num_edges = 0;
		}
		for(i = 0; i < iter->num_edges; i++){
			if(iter->edges[i]->condemned == FALSE){
				next = iter->edges[i]->a == iter ? 
										iter->edges[i]->b : iter->edges[i]->a;
				printf("%s", next->name);
				if(next->sequence != NULL){
					printf(" %s", individual);
				}
				printf("\n");
			}
		}
		printf("\n");
	}

	fclose(fp);
}

/* for debuggin only: Print the names of the nodes of the graph */
void print_graph(const graph* const graph)
{
	node* iter;
	for(iter = graph->nodelist; iter; iter = iter->next){
		printf("%s\n", iter->name);
	}
}

