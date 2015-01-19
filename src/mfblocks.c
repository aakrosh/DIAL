#include "mfblocks.h"

/* create a new mafblocks structure to add to a list */
mfblocks* new_mfnode(mafblock* const mp)
{
	mfblocks* node = ckallocz(sizeof(struct mfblocks_st));
	node->block = copy_mafblock(mp);
	return node;
}

/*sort the current link list of mafblocks using this function*/
int sort_mafblocks(const void* const a, const void* const b)
{
	mfblocks* m1 = *((mfblocks**)a);
	mfblocks* m2 = *((mfblocks**)b);

	int com1 = strcmp((char*)m1->block->name1, (char*)m2->block->name1);

	if(0 == com1){
		int com2 =  m2->block->o1 - m1->block->o2;
		if(0 == com2){
			return m1->block->s1 - m2->block->s1;
		}
		return com2;
	}

	return com1;
}

/* free  all the resources associated with the mfblocks*/
void free_mfblocks(mfblocks** plist)
{
	if(*plist != NULL){
		mfblocks* blocks = *plist;
		mfblocks* iter;
		for(iter = blocks; iter; iter = iter->next){
			ckfree(iter->block->name1);
			ckfree(iter->block->name2);
			ckfree(iter->block->aln1);
			ckfree(iter->block->aln2);
			ckfree(iter->block);
		}
		slfreelist(plist);
		*plist = NULL;
	}
}
