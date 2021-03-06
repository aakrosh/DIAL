/*! Routines for 
 * 		Error reporting 
 * 		Memory allocation/deallocation
 * 		File handling
 * 		String operations
 */

#include "utilities.h"

char* argv0;

/**---------------------Error reporting routines --------------------------*/

/*print the name of the program*/
void print_argv0()
{	
	if(argv0){
		char* p = strrchr(argv0,'/');
		fprintf(stderr,"%s: ", p ? p+1 : argv0);
	}
}

/*format the message, print it on stderr and die*/
void fatalf(const char* const fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	fflush(stdout);
	print_argv0();
	vfprintf(stderr, fmt, ap);
	fputc('\n', stderr);
	va_end(ap);
	exit(EXIT_FAILURE);
}

/*print the message on stderr and die*/
void fatal(const char* const msg)
{
	fatalf("%s", msg);
}

/**-------------------Memory handling routines -----------------------------*/

/*allocate the desired memory.*/
void* ckalloc(const size_t size)
{
	void* ptr;
	
	if ((long)size < 0){                               /* was "<= 0" -CR */
		fatal("ckalloc: request for negative space."); 
	}

	if((ptr = malloc(size)) == NULL){
		fatalf("failed to allocate %zu bytes", size);
	}

	return ptr;
}

/*allocate the desired memory and fill it with 0's*/
void* ckallocz(const size_t size)
{
	void* ptr = ckalloc(size);
	memset(ptr,0,size);
	return ptr;
}

/*reallocate the memory to the given size*/
void *ckrealloc(void * p, const size_t size)
{
	p = p ? realloc(p, size) : malloc(size);
    if (!p){
		fatal("ckrealloc failed");
	}
	return p;
}

/*free the resources. A wrapper around free*/
void ckfree(void* const p)
{
	if(p){
		free(p);
	}
}

/*allocate resources for heavy recursion */
void allocate_resources()
{
	struct rlimit rl;
	if(getrlimit(RLIMIT_STACK, &rl) != 0){
		fatal("Getting values for the system failed");
	}
	/*set the limit on the STACK*/
	rl.rlim_cur = rl.rlim_max;
	if(setrlimit(RLIMIT_STACK, &rl) != 0){
			fatal("Setting of stack failed");
	}
}

/**------------------File handling routines---------------------------------*/

/*open this file and return a file pointer*/
FILE* ckopen(const char* const name, const char* const mode)
{
	FILE* fp;
	
	if((fp = fopen(name, mode)) == NULL){
		fatalf("error in opening the file %s: %s", name, strerror(errno));
	}
	return fp;
}

/**-------------------String handling routines------------------------------*/
char* copy_string(const char* const str)
{
	char* copy = ckalloc(strlen(str)+1);
	strcpy(copy,str);
	return copy;
}

/*NOTE: this is compatible with strncmp */
int compare_names(const char* const s1, const char* const s2, const int len)
{
	int i;
	
	for(i = 0; i < len; i++){
		if(*(s1+i) != *(s2+i)){
			return *(s1+i) - *(s2+i);
		}
	}
	return 0;
}

void append(uchar** parray, 	/*the pointer to the array to be used*/
            uint* const pmax, 	/*how many bytes have been allocated*/
            uint* const plen, 	/*how many bytes have been used*/
			const int ch)		/*the byte to be added*/
{
	uint max = *pmax;
	uint len = *plen;

	/*do we need more memory*/
	if(len >= max){
		*parray = ckrealloc(*parray, max + CHUNK);
		*pmax = *pmax + CHUNK;
	}
	(*parray)[len] = ch;
	*plen = *plen + 1;
}

/*drop in replacement for 'getline' function*/
signed long getline(char** lineptr, size_t* max, FILE* stream)
{
	int ch;
	unsigned long size = 0;
	char* ptr;
	size_t sz;

	while((ch = fgetc(stream)) != EOF){
		if(size >= *max){
			sz = size + (size >> 5) + 16;
			*max = sz;
			if((ptr = ckrealloc(*lineptr, *max)) == NULL){
				return -1;
			}
			*lineptr = ptr;
		}
		(*lineptr)[size++] = ch;
		if(ch == '\n'){
			break;
		}
	}

	if(size != 0){
		if(size >= *max){
			sz = size + (size >> 5) + 16;
			*max = sz;
			if((ptr = ckrealloc(*lineptr, *max)) == NULL){
				return -1;
			}
			*lineptr = ptr;
		}
		(*lineptr)[size] = '\0';
	}
	if(0 == size || ch == EOF){
		return -1;
	}

	return size;
}

/*run the passed command in the shell*/
void run_command(const char* const command)
{
	int status;
	status = system(command);

	if(-1 == status){
		fatalf("error invoking fork during %s with %s", 
											command, strerror(errno));
	}else if(127 == status){
		fatalf("execution of the shell failed during %s", command);
	}else if(status != EXIT_SUCCESS){
		fatalf("error in execution of command:%s", command);
	}
	fprintf(stderr,"%s\n", command);
}

