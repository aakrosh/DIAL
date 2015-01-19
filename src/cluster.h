#ifndef CLUSTER_H
#define CLUSTER_H

#include "common.h"
#include "utilities.h"
#include "sequences.h"

#define MAX_DIFFS 15

typedef struct diffs_st
{
	int position;
	char reference;
	char target;
	int reads;
}diffs;

typedef struct cluster_st
{
	FILE* fp;			/* the file from which we are reading */
	char* fptr;         /* the buffer for reading */
	size_t n;           /* the size of the buffer */

	char* name;         /* name of the cluster */
	uchar* sequence;	/* the sequence of the cluster */
	uint slen;			/* length of the sequence */
	
	uint num_alloced;	/* number of members allocated for */
	uint num_members;	/* number of members of the cluster */
	char** members;		/* the members themselves */
	
	uint num_diffs;					/* number of differences in the cluster */
	diffs differences[MAX_DIFFS+1];	/* the differences themselves */
}cluster;

/* open the file corresponding to this name */
cluster* open_cluster_file(const char* const file);

/* open and read the first cluster from this file */
cluster* read_cluster_file(const char* const file);

/* read the next clutser from the file */
cluster* read_next_cluster(cluster* const cp);

/* deep copy this whole cluster */
cluster* copy_cluster(cluster* const cp);

/* free all the resources */
void close_cluster_file(cluster* cp);

/* add this member to the cluster */
void addmember(cluster* const cp, char* const name, const char* const person);

/* add this difference to the cluster */
void adddifference(cluster* const cp, const diffs difference);

/* print this cluster information */
void print_cluster(const cluster* const cp, 
				   const char* const person, 
				   const uint maxhits, 
				   FILE* const fp);
#endif
