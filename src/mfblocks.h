#ifndef MFBLOCKS_H
#define MFBLOCKS_H

#include "utilities.h"
#include "slinklist.h"
#include "hashtable.h"
#include "maf.h"

/* a list of mafblocks associated with a sequence */
typedef struct mfblocks_st
{
	struct mfblocks_st* next;
	mafblock* block;
}mfblocks;

/* create a new mafblocks structure to add to a list */
mfblocks* new_mfnode(mafblock* const mp);

/*sort the current link list of mafblocks using this function*/
int sort_mafblocks(const void* const a, const void* const b);

/* free  all the resources associated with the mfblocks*/
void free_mfblocks(mfblocks** plist);
#endif
