/* Routines for reading and writing clusters */

#include "cluster.h"

/* open the file corresponding to this name */
cluster* open_cluster_file(const char* const file)
{
	cluster* cp = ckallocz(sizeof(cluster));
	FILE* fp = strcmp(file,"/dev/stdin") == 0 ? stdin : ckopen(file,"r"); 
	
	cp->fp = fp;
	cp->n = 1;
	cp->fptr = ckalloc(cp->n+1);

	cp->name = NULL;
	cp->sequence = NULL;
	cp->slen = 0;

	cp->num_alloced = 0;
	cp->num_members = 0;
	cp->num_diffs = 0;
	return cp;
}

static rt_status read_cluster(cluster* const cp)
{
	/* lets find the next header */
	while(getline(&cp->fptr, &cp->n, cp->fp) != -1 && cp->fptr[0] != '>');

	/* is this the end of the file */
	if(feof(cp->fp) || ferror(cp->fp)){
		return ENDFILE;
	}

	/* remove the trailing \n */
	int len = strlen(cp->fptr);
	cp->fptr[len-1] = '\0';
	len--;

	/* read the name of the cluster */
	cp->name = ckrealloc(cp->name, len);
	strcpy(cp->name, cp->fptr + 1);
	
	/* read the sequence of the cluster */
	if(getline(&cp->fptr, &cp->n, cp->fp) == -1){
		fatalf("error in reading the sequence of the cluster : %s %s", 
														cp->name, cp->fptr);
	}

	/* remove the trailing \n */
	len = strlen(cp->fptr);
	cp->fptr[len-1] = '\0';
	len--;
	
	cp->sequence = ckrealloc(cp->sequence, len + 1);
	strcpy((char*)cp->sequence, cp->fptr);
	cp->slen = len;

	/* read the members of this cluster */
	uint membercount = 0;
	char* ptr = NULL;

	while(getline(&cp->fptr, &cp->n, cp->fp) != -1 &&
	      cp->fptr[0] != ' ' && 
		  cp->fptr[0] != '\n'){
		if(membercount == cp->num_alloced){
			cp->num_alloced++;
			cp->members = ckrealloc(cp->members, cp->num_alloced*sizeof(char*));			cp->members[cp->num_alloced - 1] = NULL;
		}

		ptr = strchr(cp->fptr, '\n');
		len = ptr - cp->fptr;
		ptr = ckrealloc(cp->members[membercount], len + 1);
		strncpy(ptr, cp->fptr, len);
		ptr[len] = '\0';

		cp->members[membercount] = ptr;
		membercount++;
	}
	cp->num_members = membercount;
	assert(strncmp(cp->members[0], cp->name, strlen(cp->name)) == 0);

	/* do we need to read the differences as well ? */
	int numdiffs = 0, position, reads;
	char reference, target;

	if(cp->fptr[0] == ' '){	
		if(sscanf(cp->fptr, " %d %c %c %d\n", 
							&position, &reference, &target, &reads) == 4){
			cp->differences[numdiffs].position = position;
			cp->differences[numdiffs].reference = reference;
			cp->differences[numdiffs].target = target;
			cp->differences[numdiffs].reads = reads;
			numdiffs++;
		}
		while(getline(&cp->fptr, &cp->n, cp->fp) != -1 && 
			  cp->fptr[0] != '\n' && 
			  sscanf(cp->fptr, " %d %c %c %d\n", 
			  				   &position, &reference, &target, &reads) == 4){
			cp->differences[numdiffs].position = position;
			cp->differences[numdiffs].reference = reference;
			cp->differences[numdiffs].target = target;
			cp->differences[numdiffs].reads = reads;
			numdiffs++;
			if(numdiffs > MAX_DIFFS){
                fprintf(stderr,"The number of diffs cannot exceed MAX_DIFFS\n");
				return FAILURE;
			}
		}
	}
	cp->num_diffs = numdiffs;
	return SUCCESS;
}

/* open and read the first cluster from this file */
cluster* read_cluster_file(const char* const file)
{
	/* open the cluster file */
	cluster* cp = open_cluster_file(file);

	/* read the first cluster */
	rt_status status;
	if((status = read_cluster(cp)) == ENDFILE){
		if(status == FAILURE){
			fatalf("A cluster cannot have more than MAX_DIFFS(%d) members", 
				   MAX_DIFFS);
		}
        fprintf(stderr, "%s has no clusters in it\n", file);
		return NULL;
	}

	return cp;
}

/* read the next clutser from the file */
cluster* read_next_cluster(cluster* const cp)
{
	rt_status status;
	if((status = read_cluster(cp)) == ENDFILE){
		if(status == FAILURE){
			fatalf("A cluster cannot have more than MAX_DIFFS(%d) members", 
				   MAX_DIFFS);
		}

		return NULL;
	}

	return cp;
}

/* deep copy this whole cluster */
cluster* copy_cluster(cluster* const cp)
{
	cluster* copy = ckalloc(sizeof(cluster));

	copy->fp = NULL;
	copy->fptr = NULL;
	copy->n = 0;

	copy->name = copy_string(cp->name);
	copy->sequence = (uchar*)copy_string((char*)cp->sequence);
	copy->slen = cp->slen;

	copy->num_alloced = cp->num_members;
	copy->num_members = cp->num_members;
	copy->members = ckalloc(copy->num_members*sizeof(char*));

	uint i;
	for(i = 0; i < copy->num_members; i++){
		copy->members[i] = copy_string(cp->members[i]);
	}

	copy->num_diffs = cp->num_diffs;
	for(i = 0; i < copy->num_diffs; i++){
		copy->differences[i].position = cp->differences[i].position;
		copy->differences[i].reference = cp->differences[i].reference;
		copy->differences[i].target = cp->differences[i].target;
		copy->differences[i].reads = cp->differences[i].reads;
	}
	return copy;
}

/* free all the resources */
void close_cluster_file(cluster* cp)
{
	ckfree(cp->name);
	ckfree(cp->sequence);
	ckfree(cp->fptr);
	if(cp->fp != stdin){
		fclose(cp->fp);
	}
	ckfree(cp->members);
	ckfree(cp);
	*(&cp) = NULL;
}

/* add this member to the cluster */
void addmember(cluster* const cp, char* const name, const char* const person)
{
	if(cp->sequence[0] == '?'){
		return;
	}

	/* does this member exist already ? */
	uint i;
	for(i = 0; i < cp->num_members; i++){
		if(strncmp(cp->members[i], name, strlen(name)) == 0){
			return;
		}
	}

	if(cp->num_members == cp->num_alloced){
		cp->num_alloced++;
		cp->members = ckrealloc(cp->members, cp->num_alloced*sizeof(char*));
		cp->members[cp->num_alloced - 1] = NULL;
	}
	char* ptr = ckrealloc(cp->members[cp->num_members], strlen(name) + strlen(person) +  8);
	sprintf(ptr, "%s\t%s", name, person);

	cp->members[cp->num_members++] = ptr;
}

/* add this difference to the cluster */
void adddifference(cluster* const cp, const diffs difference)
{
	if(cp->sequence[0] == '?'){
		return;
	}

	/* does this difference exist already ? */
	uint i = 0;
	for(i = 0; i < cp->num_diffs; i++){	
		if(cp->differences[i].position == difference.position &&
		   cp->differences[i].reference == difference.reference &&
		   cp->differences[i].target == difference.target){
			cp->differences[i].reads++;
			return;
		}
	}

	/* add this new difference */
	if(i > MAX_DIFFS){
		return;
	}

	cp->differences[i].position = difference.position;
	cp->differences[i].reference = difference.reference;
	cp->differences[i].target = difference.target;
	cp->differences[i].reads = 1;
	cp->num_diffs++;
}

/* print this cluster information */
void print_cluster(const cluster* const cp, 
				   const char* const person, 
				   const uint maxhits, 
				   FILE* const fp)
{
	/* name */
	fprintf(fp,">%s\n", cp->name);

	/* sequence */
	if(cp->num_members > maxhits || cp->num_diffs > MAX_DIFFS){
		fprintf(fp,"??");
	}
	fprintf(fp,"%s\n", cp->sequence);

	/* members of the cluster */
	if(cp->num_members <= maxhits && cp->num_diffs <= MAX_DIFFS){
		uint i;
		for(i = 0; i < cp->num_members; i++){
			fprintf(fp,"%s\n", cp->members[i]);
		}
	}else{
		fprintf(fp,"%s\t%s\n", cp->name, person);
	}
	
	/* differences */
	if(cp->num_members <= maxhits && cp->num_diffs <= MAX_DIFFS){
		uint i;
		for(i = 0; i < cp->num_diffs; i++){
			fprintf(fp," %d %c %c %d\n", 
				         cp->differences[i].position, 
					 	 cp->differences[i].reference, 
					 	 cp->differences[i].target, 
					 	 cp->differences[i].reads);
		}
	}
	fprintf(fp,"\n");
}
