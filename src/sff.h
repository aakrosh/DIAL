#ifndef SFF_H
#define SFF_H

#include <wchar.h>

#include "utilities.h"
#include "slinklist.h"
#include "hashtable.h"

typedef struct common_header_st
{
	uint64_t 	index_offset;
	uint32_t	index_length;
	uint32_t	num_reads;
	uint16_t	header_length;
	uint16_t	key_length;
	uint16_t	num_flows_per_read;
	char* 		flow_chars;	 /*This will be read as is and not written as is*/
	char*		key_sequence;/*This will be read as is and not written as is*/
}common_header;

typedef struct read_header_st
{
	uint16_t 	header_length;
	uint16_t 	name_length;
	uint32_t 	num_bases;
	uint16_t 	clip_qual_left;
	uint16_t	clip_qual_right;
	uint16_t	clip_adapter_left;
	uint16_t	clip_adapter_right;
	uint16_t    num_flows_per_read;
	char*		name;
}read_header;

typedef struct read_data_st
{
	uint16_t* 	flowgram_values;
	uint8_t*	flow_index_per_base;
	char* 		bases;
	uint8_t*	quality_scores;
}read_data;

typedef struct position_st
{
	char seq;
	char qual;
}position;

typedef struct ancilliary_st
{
	char* header;		/*the accession value of the read*/
	bool is_used;		/*has this read been used in this contig*/
	bool complement;	/*is the read complemented*/
	int offset;			/*0-based offset of the read*/
	unsigned int begin;	/*0-based start position on read used in alignment*/
	unsigned int end;	/*0-based end position on read used in alignment*/
	unsigned int num_bases;		/*number of bases in the read*/
	position* sequence;			/*the sequence and the quality values*/
}ancilliary;

typedef struct sff_read_st
{
	struct sff_read_st* next;
	read_header*        header;
	read_data*          data;
	ancilliary*         information;
	char*               individual;
}sff_read;

typedef struct sff_st
{
	common_header common_header;
	sff_read*     reads;
}sff;

#define SFF_MAGIC   0x2E736666
#define SFF_VERSION 0x00000001

/*read the sff file. If the hash is NULL, then just return*/
common_header* sff_read_file(const char* const name, hashtable* const hash);

/* fill in the info about the next read */
void read_sff_rinfo(sff_read** const pread, 
					FILE* const fp,     
					const int flows_per_read);

/*read the common header section of the file specified by fp*/
void read_common_header(FILE* const fp, common_header* const header);

/*clear the is_used flag in the sff_read structure*/
void clear_used_flags(sff_read* const reads);

/*write a sff file with this structure*/
void sff_write_file(const char* const file_name, 
					const common_header common_header, 
					const sff_read* const reads);

/*check the integrity of the sff file*/
void check_integrity(const char* const file_name, const sff_read* const reads);

/*free the resources allocated for this read structure*/
void sff_read_free(sff_read** pread);

/*free resources allocated for this SFF data structure*/
void sff_free(sff** psff);

#endif
