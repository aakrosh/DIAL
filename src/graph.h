#ifndef GRAPH_H
#define GRAPH_H

#include "common.h"
#include "utilities.h"
#include "sequences.h"
#include "slinklist.h"
#include "hashtable.h"

#define NUM_EDGES 5

typedef struct graph_st
{
	hashtable*       nodes;
	struct node_st*  nodelist;
	struct edge_st*  edgelist;
}graph;

typedef struct node_st
{
	struct node_st* next;			/* next node in the list */
	char*            name; 			/* name of the node. normally a ptr to 
									 * the name in the hashtable */
	uchar*           sequence;		/* the sequence of the node */
	bool 	         condemned;		/* does this node have more than MAX_HITS*/
	bool             visited;		/* for DFS. has this been visited ? */
	int 		     component;		/* which connected component is this from*/

	int              num_alloced;	/* number of edges allocated for */
	int              num_edges;		/* number of edges */
	struct edge_st** edges;			/* the edges */ 
}node;

typedef struct edge_st
{
	struct edge_st* next;
	node*      a;					/* one of the nodes */
	node*      b;					/* the second node */
	uint16_t   weight;				/* weight of the edge */
	bool       condemned;			/* should traversals consider this edge*/
}edge;

/* return a new graph */
graph* new_graph();

/* add a new node to the graph */
node* add_node(graph* const graph, 
	   	       const char* const name,
			   const uchar* const seq);

/* add the name of the nodes to the graph. Return the number of nodes added */
int add_nodes(graph* const graph, const char* const seqfile);

/* does this edge already exist*/
bool edge_exists(const node* const n1, const node* const n2);

/* make the edge between the two nodes for this graph */
edge* make_edges(graph* const graph, 
				 const char* const name1,
				 const char* const  name2,
				 const int overlap);

/* free the resources of the graph */
void free_graph(graph** pgraph);

/* find and mark the connected comonents in this graph. Return the number of
 * connected components in the graph */
int find_connected_components(graph* const graph);

/* remove all the edges connected to this node in this graph */
void remove_all_edges(graph* const graph, node* const node);

/* go through each of the nodes. If a node has more than maxhits edges, then
 * all those edges should be deleted. The node should be marked as a repeat.
 * Return the number of nodes filtered. Write those nodes in the repeatfile */
int filter_nodes(graph* const graph, FILE* const fp, const int maxhits);

/* sort the nodes based on the component they are in */
int component_sort(const void* const a, const void* const b);

/* sort the nodes based on the connected component they are in. Then print the
 * clusters. Every node in a component will act as a starting point for the
 * cluster. This is so that we dont have to deal with extending the consensus
 * sequence of the clutser. This will also allow us to keep on using the
 * MAX_HITS property for each new read in subsequence. */
void print_graph_clusters(graph* const graph, 
						  const char* const clusterfile,
					      const int index,
					      const char* const individual,
						  const int maxhits);

/* for debuggin only: Print the names of the nodes of the graph */
void print_graph(const graph* const graph);
#endif
