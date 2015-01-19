/*! Routines for 
 *  	Error reporting
 *  	Memory allocation
 *  	File Handling
 *  	String operations
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#include <sys/resource.h>
#ifndef __USE_BSD
#define __USE_BSD
#endif
#include <stdlib.h>		
#include <stdio.h>		/*io functions*/
#include <string.h>		/*string functions*/
#include <stdarg.h>		/*va_list and other structures*/
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#undef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))

/*the memory is incremented in chunks of bytes specified by this macro*/
#define CHUNK 60

/*the return value from a routine*/
typedef enum rt_status_st {ENDFILE,SUCCESS,FAILURE}rt_status;

typedef enum boolean{FALSE = 0, TRUE} bool;

typedef unsigned int uint;

typedef unsigned char uchar;

/*the program using this header*/
extern char* argv0;

/*print the name of the program*/
void print_argv0();

/*error reporting routines*/
void fatal(const char* const msg);
void fatalf(const char* const fmt, ...);

/*memory allocation routines*/
void* ckalloc(const size_t size);
void* ckallocz(const size_t size);
void *ckrealloc(void * p, const size_t size);
void ckfree(void* const p);

/*resource allocation routine*/
void allocate_resources();

/*file handling routines*/
FILE* ckopen(const char* const name, const char* const mode);

/*string handling routines*/
char* copy_string(const char* const str);
int compare_names(const char* const s1,const char* const s2, const int len);
void append(uchar** parray, 	/*the pointer to the array to be used*/
            uint* const pmax, 	/*how many bytes have been allocated*/
		    uint* const plen, 	/*how many bytes have been used*/
			const int ch);		/*the byte to be added*/

/*drop in replacement for the getline function*/
signed long getline(char** lineptr, size_t* max, FILE* stream);

/*run the passed command in the shell*/
void run_command(const char* const command);
#endif
