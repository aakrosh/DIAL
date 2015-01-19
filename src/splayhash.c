#include "splayhash.h"

/*allocate a new splayhash*/
splayhash* new_splayhash(const int po2size)
{
	int max_size = po2size;
	if(po2size > 24){
		fprintf(stderr,"hash po2size should not exceed 24, using 24\n");
		max_size = 24;
	}

	splayhash* splayhash = ckallocz(sizeof(struct splayhash_st));
	splayhash->po2size = max_size;
	splayhash->size = (1 << splayhash->po2size);
	splayhash->mask = splayhash->size - 1;

	splayhash->bins  = ckallocz(splayhash->size*sizeof(struct bin_st*));
	return splayhash;
}

/*The splay function which is used from
 * http://www.link.cs.cmu.edu/link/ftp-site/splaying/top-down-splay.c*/
static splaybin* splay(const kmer name, splaybin* tree)
{
	assert(tree != NULL);
    splaybin N, *l, *r, *y;
    N.left = N.right = NULL;
    l = r = &N;

    for (;;) {
	if (name < tree->name) {
	    if (tree->left != NULL && name < tree->left->name) {
		y = tree->left; tree->left = y->right; y->right = tree; tree = y;
	    }
	    if (tree->left == NULL) break;
	    r->left = tree; r = tree; tree = tree->left;
	} else if (name > tree->name) {
	    if (tree->right != NULL && name > tree->right->name) {
		y = tree->right; tree->right = y->left; y->left = tree; tree = y;
	    }
	    if (tree->right == NULL) break;
	    l->right = tree; l = tree; tree = tree->right;
	} else break;
    }
    l->right=tree->left; r->left=tree->right; 
	tree->left=N.right; tree->right=N.left;
    return tree;
}

/*look up the hash table and return the bin if it does exist*/
splaybin* lookup_splayhash(splayhash* const hash, 
				           const kmer name)
{
	uint32_t index =  superfasthash((char*)&name, sizeof(name));
	index = index & hash->mask;

	if(hash->bins[index] == NULL){
		return NULL;
	}

	hash->bins[index] = splay(name, hash->bins[index]);
	if(hash->bins[index]->name == name){
		return hash->bins[index];
	}

	return NULL;
}

/*add the following and return the bin corresponding to it*/
splaybin* add_splayhash(splayhash* const hash,  	/*the hashtable*/
					    const kmer name,			/*the string*/
				   		const int counts,	    	/*the index (0-based)*/
						const int score)			/*the score(0-based)*/
{
	splaybin* bin = ckalloc(sizeof(struct splaybin_st));
	bin->name = name;
	bin->counts = counts;
	bin->score = score;
	bin->left = NULL;
	bin->right = NULL;

	uint32_t index =  superfasthash((char*)&name, sizeof(name));
	index = index & hash->mask;

	splaybin* tree = hash->bins[index];
	if(tree != NULL){	
		tree = splay(name, tree);
		if(name < tree->name){
			bin->left = tree->left;
			bin->right = tree;
			tree->left = NULL;
		}else if(name > tree->name){
			bin->right = tree->right;
			bin->left = tree;
			tree->right = NULL;
		}else{
			fatalf("should not add this to hash again");
		}
	}
	hash->bins[index] = bin;	
	hash->elcount++;

	return bin;
}

static void free_bins(splaybin* bin)
{
	if(bin->left != NULL){
		free_bins(bin->left);
	}
	if(bin->right != NULL){
		free_bins(bin->right);
	}
	bin->left = NULL;
	bin->right = NULL;
	ckfree(bin);
}

/* free all the resources held by the splayhash */
void free_splayhash(splayhash** phash)
{
	splayhash* hash = *phash;

	int i;
	for(i = 0; i < hash->size; i++){
		if(hash->bins[i] != NULL){
			free_bins(hash->bins[i]);
		}
	}
	ckfree(hash->bins);
	ckfree(hash);
	*phash = NULL;
}

