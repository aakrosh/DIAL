/* Routines for reading and writing maf files */

#include "maf.h"

/* open the maf file corresponding to this name */
mafblock* open_maf_file(const char* const file)
{
	mafblock* mp = ckallocz(sizeof(mafblock));
	FILE* fp = strcmp(file,"/dev/stdin") == 0 ? stdin : ckopen(file,"r"); 

	mp->fp = fp;
	mp->n = 1;
	mp->fptr = ckalloc(mp->n+1);

	return mp;
}

/* read the block */
static rt_status read_block(mafblock* const mp)
{
	/* is it the end of te file */
	if(feof(mp->fp) || ferror(mp->fp)){
		return ENDFILE;
	}

	/* some temporary holders */
	char score[128];
	char name[128];

	/* this should point to the next block */
	assert(mp->fptr[0] == 'a');
	if(sscanf(mp->fptr,"a %s\n", score) != 1){
		fatalf("error in reading the score for this block:%s", mp->fptr);
	}
	mp->score = atoi(score+6);
	
	/* read the stuff for the first sequence */
	if(getline(&mp->fptr, &mp->n, mp->fp) == -1){
		fatalf("error in reading the first sequence details:%s", mp->fptr);
	}
	assert(mp->fptr[0] == 's');

	if(sscanf(mp->fptr,"s %s %d %d %c %d ", 
					    name, &mp->s1, &mp->o1, &mp->strand1, &mp->len1) != 5){
		fatalf("error in parsing details of the first sequence:%s", mp->fptr);
	}
		
	mp->name1 = ckrealloc(mp->name1, strlen(name)+1);
	strcpy((char*)mp->name1, name);

	char* ptr = strrchr(mp->fptr, ' ');
	ptr++;
	do{
		append(&mp->aln1, &mp->amax, &mp->alen, *ptr);
		ptr++;
	}while(*ptr != '\n');
	append(&mp->aln1, &mp->amax, &mp->alen, 0);
	mp->alen--;
	
	/*read the stuff for the second sequence*/
	if(getline(&mp->fptr, &mp->n, mp->fp) == -1){
		fatalf("error in reading the second sequence details:%s", mp->fptr);
	}
	assert(mp->fptr[0] == 's');

	if(sscanf(mp->fptr,"s %s %d %d %c %d ", 
					   name, &mp->s2, &mp->o2, &mp->strand2, &mp->len2) != 5){
		fatalf("error in parsing details of the second sequence:%s", mp->fptr);
	}

	mp->name2 = ckrealloc(mp->name2, strlen(name) + 1);
	strcpy((char*)mp->name2, name);

	ptr = strrchr(mp->fptr, ' ');
	ptr++;
	do{
		append(&mp->aln2, &mp->bmax, &mp->blen, *ptr);
		ptr++;
	}while(*ptr != '\n');
	append(&mp->aln2, &mp->bmax, &mp->blen, 0);
	mp->blen--;

	/*move to the next block*/
	while(getline(&mp->fptr, &mp->n, mp->fp) != -1 && mp->fptr[0] != 'a'){
	}

	return SUCCESS;
}

/* open the maf file and return the first block */
mafblock* read_maf_file(const char* const file)
{
	/* open the maf file */
	mafblock* mp = open_maf_file(file);

	/*skip the header for the maf file */
	while(getline(&mp->fptr, &mp->n, mp->fp) != -1 &&
		  mp->fptr[0] == '#'){
	}

	/* read the block */
	if(read_block(mp) == ENDFILE){
		return NULL;
	}

	return mp;
}

/* get the next maf block from the file */
mafblock* get_next_block(mafblock* const mp)
{
	mp->alen = 0;
	mp->blen = 0;

	if(read_block(mp) == ENDFILE){
		return NULL;
	}

	return mp;
}

/*clone a maf block*/
mafblock* copy_mafblock(const mafblock* const mp)
{
	mafblock* mafblock = ckallocz(sizeof(struct mafblock_st));

	mafblock->score = mp->score;

	mafblock->name1   = (uchar*)copy_string((char*)mp->name1);
	mafblock->s1      = mp->s1;
	mafblock->o1      = mp->o1;
	mafblock->strand1 = mp->strand1;
	mafblock->len1    = mp->len1;
	mafblock->aln1    = (uchar*)copy_string((char*)mp->aln1);
	mafblock->amax    = mp->alen;
	mafblock->alen    = mp->alen;

	mafblock->name2   = (uchar*)copy_string((char*)mp->name2);
	mafblock->s2      = mp->s2;
	mafblock->o2      = mp->o2;
	mafblock->strand2 = mp->strand2;
	mafblock->len2    = mp->len2;
	mafblock->aln2    = (uchar*)copy_string((char*)mp->aln2);
	mafblock->bmax    = mp->blen;
	mafblock->blen    = mp->blen;

	return mafblock;
}

/* just print this mafblock */
void print_mafblock(const mafblock* const mp)
{
	printf("a score=%d\n", mp->score);
	printf("s %s %d %d %c %d %s\n", mp->name1, 
									mp->s1, 
									mp->o1, 
									mp->strand1, 
									mp->len1, 
									mp->aln1);
	printf("s %s %d %d %c %d %s\n", mp->name2, 
									mp->s2,
									mp->o2, 
									mp->strand2, 
									mp->len2, 
									mp->aln2);
	printf("\n");
}

/* free all the resources used by this maf block */
void close_maf_file(mafblock* mp)
{
	ckfree(mp->fptr);
	ckfree(mp->name1);
	ckfree(mp->aln1);
	ckfree(mp->name2);
	ckfree(mp->aln2);
	if(mp->fp != stdin){
		fclose(mp->fp);
	}
	ckfree(mp);
	*(&mp) = NULL;
}
