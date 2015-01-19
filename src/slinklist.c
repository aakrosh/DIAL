#include "slinklist.h"

/*return the number of elements in the list*/
int slcount(const void* const list)
{
	int count = 0;
	slinklist* iter = (slinklist*) list;

	for(; iter; ++count, iter = iter->next);
	return count;
}

/*reverse order of  a list*/
void slreverse(void* plist)
{
	slinklist** ppt = (slinklist**)plist;
	slinklist* newlist = NULL;
	slinklist* el = NULL;
	slinklist* next = NULL;

	next = *ppt;
	while(next != NULL){
		el = next;
		next = el->next;
		el->next = newlist;
		newlist = el;
	}
	*ppt = newlist;
}

/*remove the item from the list*/
void* slremove(void* plist, void* const node)
{
	slinklist* iter = *((slinklist**) plist);
	slinklist* t = (slinklist*) node;
	slinklist* pt = (slinklist*) node;

	if(iter == t){
		*(slinklist**)plist = t->next;
		return t;
	}
	for(; iter && (iter != t); pt = iter, iter = iter->next);	
	pt->next = t->next;
	return t;
}

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist)
{
	slinklist** ppt = (slinklist**)plist;
	slinklist* next = *ppt;
	slinklist* el = NULL;

	while(next != NULL){
		el = next;
		next = el->next;
		ckfree((char*)el);
	}
	*ppt = NULL;
}

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, 
		    int(*compare)(const void* const elem1, const void* const elem2))
{
	slinklist** pl = (slinklist**)plist;
	slinklist* list = *pl;

	int count = slcount(list);
	if(count > 1){
		slinklist* el = NULL;
		slinklist** array = NULL;
		int i = 0;

		array = ckalloc(count * sizeof(*array));
		for(el = list, i=0; el != NULL; el = el->next, i++){
			array[i] = el;
		}
		qsort(array, count, sizeof(array[0]), compare);
		list = NULL;
		for( i = 0; i < count; i++){
			array[i]->next = list;
			list = array[i];
		}
		ckfree(array);
		slreverse(&list);
		*pl = list;
	}
}
