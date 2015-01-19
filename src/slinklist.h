#ifndef SLINKLIST_H
#define SLINKLIST_H

#include "utilities.h"

/*abstraction for singly linked lists*/
typedef struct slinklist_st
{
	struct slinklist_st* next;
}slinklist;

#define sladdhead(plist, node) ((node)->next = *(plist), *(plist) = node)

/*return the number of elements in the list*/
int slcount(const void* const list);

/*reverse order of  a list*/
void slreverse(void* plist);

/*remove the item from the list*/
void* slremove(void* plist, void* const node);

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist);

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, 
		    int(*compare)(const void* const elem1, const void* const elem2));
#endif
