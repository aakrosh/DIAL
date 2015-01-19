#include "sff.h"

bool need_swapped = FALSE;

inline static uint16_t BE2(const uint16_t number)
{
	uint16_t value = number;
	if(need_swapped){
		value = ((number & 0x00FF) << 8) + ((number & 0xFF00) >> 8);
	}
	return value;
}

inline static uint32_t BE4(const uint32_t number)
{
	uint32_t value = number;
	if(need_swapped){
		value = ((number & 0x000000FF) << 24)	+ 	
				((number & 0x0000FF00) << 8)	+ 	
				((number & 0x00FF0000) >> 8)	+ 	
				((number & 0xFF000000) >> 24);
	}
	return value;
}

inline static uint64_t BE8(const uint64_t number)
{
	uint64_t value = number;
	if(need_swapped){
		value = ((number & 0x00000000000000FF) << 56) + 
			    ((number & 0x000000000000FF00) << 40) + 
			    ((number & 0x0000000000FF0000) << 24) + 
			    ((number & 0x00000000FF000000) << 8) + 
			    ((number & 0x000000FF00000000) >> 8) + 
			    ((number & 0x0000FF0000000000) >> 24) + 
			    ((number & 0x00FF000000000000) >> 40) + 
			    ((number & 0xFF00000000000000) >> 56);
	}
	return value;
}

static uint16_t read_two_bytes(FILE* const fp)
{
	uint16_t holder;
	if(fread(&holder, 2, 1, fp) != 1){
		if(feof(fp) != 0){
			fatal("error in reading 2 bytes from the stream, eof reached");
		}else if(ferror(fp) != 0){
			fatal("error in reading 2 bytes from the stream, error");
		}
	}
	holder = BE2(holder);
	return holder;
}

static uint32_t read_four_bytes(FILE* const fp)
{
	uint32_t holder;
	if(fread(&holder, 4, 1, fp) != 1){
		if(feof(fp) != 0){
			fatal("error in reading 4 bytes from the stream, eof reached");
		}else if(ferror(fp) != 0){
			fatal("error in reading 4 bytes from the stream, error");
		}
	}
	holder = BE4(holder);
	return holder;
}

static uint64_t read_eight_bytes(FILE* const fp)
{
	uint64_t holder;
	if(fread(&holder, 8, 1, fp) != 1){
		if(feof(fp) != 0){
			fatal("error in reading 8 bytes from the stream, eof reached");
		}else if(ferror(fp) != 0){
			fatal("error in reading 8 bytes from the stream, error");
		}
	}
	holder = BE8(holder);
	return holder;
}

/*given the char array, its size and a file pointer, this function reads size 
 * number of bytes from the file and puts it into the array. The reading is 
 * done byte by byte and since the sff files are big endian, we can directly 
 * put the stuff in the array*/
static void read_array(char* const array, const uint16_t size, FILE* const fp)
{	
	uint8_t one;
	uint16_t i;
	for(i = 0; i < size; i++){
		if(fread(&one, 1, 1, fp) != 1){
    		if(feof(fp) != 0){
    			fatal("error in reading 1 byte from the stream, eof reached");
    		}else if(ferror(fp) != 0){
    			fatal("error in reading 1 byte from the stream, error");
    		}
        }
		array[i] = one;
	}
}

static void read_large_array(char* const array, const uint32_t size, FILE* const fp)
{	
	uint8_t one;
	uint32_t i;
	for(i = 0; i < size; i++){
		if(fread(&one, 1, 1, fp) != 1){
    		if(feof(fp) != 0){
    			fatal("error in reading 1 byte from the stream, eof reached");
    		}else if(ferror(fp) != 0){
    			fatal("error in reading 1 byte from the stream, error");
    		}
        }
		array[i] = one;
	}
}

void sff_read_free(sff_read** pread)
{
	ckfree((*pread)->header->name);
	ckfree((*pread)->data->flowgram_values);
	ckfree((*pread)->data->flow_index_per_base);
	ckfree((*pread)->data->bases);
	ckfree((*pread)->data->quality_scores);
	ckfree((*pread)->header);
	ckfree((*pread)->data);
	ckfree((*pread)->information->sequence);
	ckfree((*pread)->information);
}

/*free resources allocated for this SFF data structure*/
void sff_free(sff** psff)
{
	ckfree((*psff)->common_header.flow_chars);
	ckfree((*psff)->common_header.key_sequence);
	
	/*free all the resources used by the reads*/
	sff_read* read = NULL;
	for(read = (*psff)->reads; read; read = read->next){
		sff_read_free(&read);
	}

	slfreelist(&(*psff)->reads);
	ckfree(*psff);
}

void read_common_header(FILE* const fp, common_header* const header)
{
	uint32_t magic, version;

	/*generic holders for the various sizes*/
	uint8_t one;

	/*check the magic number for the header*/
	if(fread(&magic,4,1,fp) != 1){
		if(feof(fp) != 0){
			fatal("error in reading 4 bytes from the stream, eof reached");
		}else if(ferror(fp) != 0){
			fatal("error in reading 4 bytes from the stream, error");
		}
    }
	need_swapped = (magic == SFF_MAGIC) ?  FALSE : TRUE;
	magic = BE4(magic);
	if(magic != SFF_MAGIC){
		wprintf(L"%" PRIX32 "\n",magic);
		fatalf("wrong magic number in the file");
	}

	/*check the version number*/
	if((version = read_four_bytes(fp)) != SFF_VERSION){
		wprintf(L"%" PRIX32 "\n",version);
		fatalf("wrong version number in the file");
	}

	/*put the remainder of the values in the header*/
	header->index_offset = read_eight_bytes(fp);
	header->index_length = read_four_bytes(fp);
	header->num_reads = read_four_bytes(fp);
	header->header_length = read_two_bytes(fp);
	header->key_length = read_two_bytes(fp);
	header->num_flows_per_read = read_two_bytes(fp);
    if(fread(&one, 1, 1, fp) != 1){
        if(feof(fp) != 0){
    	    fatal("error in reading 1 byte from the stream, eof reached");
    	}else if(ferror(fp) != 0){
    		fatal("error in reading 1 byte from the stream, error");
    	}
    }
	assert(one == 0x01);
	
	header->flow_chars = ckallocz(header->num_flows_per_read + 1);
	read_array(header->flow_chars, header->num_flows_per_read, fp);
	header->flow_chars[header->num_flows_per_read] = '\0';
	assert(header->num_flows_per_read % 8 == 0);

	header->key_sequence = ckallocz(header->key_length + 1);
	read_array(header->key_sequence, header->key_length, fp);
	header->key_sequence[header->key_length] = '\0';

	assert((31 + header->num_flows_per_read + header->key_length) < INT_MAX);
	int length =  31 + header->num_flows_per_read + header->key_length;
	fseek(fp, length % 8 == 0 ? 0 : 8 - (length % 8),SEEK_CUR);
}

void read_sff_rinfo(sff_read** const pread, 
	 			    FILE* const fp, 
					const int flows_per_read)
{
	read_header* header = (*pread)->header;
	header->header_length = read_two_bytes(fp);
	header->name_length = read_two_bytes(fp);
	header->num_bases = read_four_bytes(fp);
	header->clip_qual_left = read_two_bytes(fp);
	header->clip_qual_right = read_two_bytes(fp);
	header->clip_adapter_left = read_two_bytes(fp);
	header->clip_adapter_right = read_two_bytes(fp);

	header->name = ckalloc(header->name_length + 1);
	read_array(header->name, header->name_length, fp);
	header->name[header->name_length] = '\0';

	assert((16 + header->name_length) < INT_MAX);
	int length = 16 + header->name_length;
	fseek(fp, 8 - (length % 8), SEEK_CUR);

	/*read the data section for this read*/
	read_data* data = (*pread)->data;
	int num_bases = header->num_bases;

	data->flowgram_values = ckallocz(flows_per_read*2);
	read_array((char*)data->flowgram_values, flows_per_read*2, fp);

	data->flow_index_per_base = ckallocz(num_bases);
	read_large_array((char*)data->flow_index_per_base, num_bases, fp);

	data->bases = ckallocz(num_bases);
	read_large_array(data->bases, num_bases, fp);

	data->quality_scores = ckallocz(num_bases);
	read_large_array((char*)data->quality_scores, num_bases, fp);

	assert((2*flows_per_read + 3*num_bases) < INT_MAX);
	length = 2*flows_per_read + 3*num_bases;
	fseek(fp, length%8 == 0 ? 0 : 8 - (length % 8), SEEK_CUR);
}

/*read the sff file. If the hash is NULL, then just return */
common_header* sff_read_file(const char* const name, hashtable* const hash)
{
	sff_read* read  = NULL;
	common_header* common_header = ckalloc(sizeof(struct common_header_st));

	FILE* fp = ckopen(name,"r");
	bin* hel = NULL;

	/*read the common header for the sff file*/
	read_common_header(fp, common_header);
	int flows_per_read = common_header->num_flows_per_read;

	/*read the remainder of the data*/
	uint32_t i = 0;
	for(i = 0; i < common_header->num_reads; i++){
		read = ckallocz(sizeof(sff_read));
		read->header = ckallocz(sizeof(read_header));
		read->data = ckallocz(sizeof(read_data));	
		read->information = ckallocz(sizeof(ancilliary));
		read->header->num_flows_per_read = flows_per_read;

		read_sff_rinfo(&read, fp, flows_per_read);
		if((hel = lookup_hashtable(hash,read->header->name,
								  strlen(read->header->name))) != NULL){
			hel->val = read;
		}else{
			sff_read_free(&read);
			ckfree(read);	
		}
	}
	fclose(fp);
	return common_header;
}

static void write_two_bytes(FILE* const fp, const uint16_t number)
{
	uint16_t holder = number;
	if(need_swapped){
		holder = BE2(holder);
	}
	fwrite(&holder, 2, 1, fp);
}

static void write_four_bytes(FILE* const fp, const uint32_t number)
{
	uint32_t holder = number;
	if(need_swapped){
		holder = BE4(holder);
	}
	fwrite(&holder, 4, 1, fp);
}

static void write_eight_bytes(FILE* const fp, const uint64_t number)
{
	uint64_t holder = number;
	if(need_swapped){
		holder = BE8(holder);
	}
	fwrite(&holder, 8, 1, fp);
}

/*write the common header for the sff file into this file*/
static void write_common_header(common_header header, 
                                FILE* const fp, 
								const uint32_t num_reads, 
                                const uint64_t indexoffset)
{
	uint8_t code = 0x01;
	uint8_t zero = 0x00;

	write_four_bytes(fp,SFF_MAGIC);
	write_four_bytes(fp,SFF_VERSION);
	write_eight_bytes(fp, indexoffset);
	write_four_bytes(fp, 0x00000000);
	write_four_bytes(fp, num_reads);
	write_two_bytes(fp, header.header_length);
	write_two_bytes(fp, header.key_length);
	write_two_bytes(fp, header.num_flows_per_read);
	fwrite(&code, 1, 1, fp);
    fprintf(fp, "%s", header.flow_chars);
    fprintf(fp, "%s", header.key_sequence);

	/*pad the zeroes*/
	int length = 31 + header.num_flows_per_read + header.key_length;
	int zeroes = header.header_length - length;
	assert(zeroes == ((length % 8 == 0) ? 0 : 8 - (length % 8)));

	int i;
	for(i = 0; i < zeroes; i++){
		fwrite(&zero, 1, 1, fp);
	}
}

/*write the read header followed by the read data into the file*/
static void write_sff_rinfo(const sff_read* const read, FILE* const fp,
							const common_header cheader)
{
	uint8_t zero = 0x00;
	read_header* header = read->header;
	read_data*   data   = read->data;

	/*write the header section*/
	write_two_bytes(fp, header->header_length);
	write_two_bytes(fp, header->name_length);
	write_four_bytes(fp, header->num_bases);
	write_two_bytes(fp, header->clip_qual_left);
	write_two_bytes(fp, header->clip_qual_right);
	write_two_bytes(fp, header->clip_adapter_left);
	write_two_bytes(fp, header->clip_adapter_right);
	fwrite(header->name, 1, header->name_length, fp);

	/*pad the zeroes*/
	int length = 16 + header->name_length;
	int zeroes = header->header_length - length;
	assert(zeroes == ((length % 8 == 0) ? 0 : 8 - (length % 8)));
	int i;
	for(i = 0; i < zeroes; i++){
		fwrite(&zero, 1, 1, fp);
	}

	/*write the data section*/
	fwrite(data->flowgram_values, 1, 2*cheader.num_flows_per_read, fp);
	fwrite(data->flow_index_per_base, 1, header->num_bases, fp);
	fwrite(data->bases, 1, header->num_bases, fp);
	fwrite(data->quality_scores,1,header->num_bases, fp);

	length = 2*cheader.num_flows_per_read + 3*header->num_bases;
	zeroes = (length % 8 == 0) ? 0 : 8 - (length % 8);
	for(i = 0; i < zeroes; i++){
		fwrite(&zero, 1, 1, fp);
	}
}

/*write a SFF file with this structure*/
void sff_write_file(const char* const file_name,
					const common_header common_header,
					const sff_read* reads)
{
	if(reads == NULL){
		return;
	}
	FILE* fp = ckopen(file_name, "w");

	/*common-header*/
    uint64_t offset = common_header.header_length;  /* offset of the index */
	uint32_t num_reads = 0;
    uint64_t datalength = 0;
	const sff_read* read;

	for(read = reads; read; read = read->next){
		num_reads++;

        offset += read->header->header_length;

        datalength  = 2*common_header.num_flows_per_read;
        datalength += 3*read->header->num_bases;

        offset += datalength;
        offset += (datalength % 8 == 0) ? 0 : 8 - (datalength % 8);
	}

	write_common_header(common_header, fp, num_reads, offset);

	/*write the remainder of the data*/
	for(read = reads; read; read = read->next){
		write_sff_rinfo(read, fp, common_header);
	}

	fclose(fp);
}

/*clear the is_used flags for the SFF_READ in this list*/
void clear_used_flags(sff_read* const reads)
{
	if(reads == NULL){
		return;
	}
	sff_read* read = NULL;
	for(read = reads; read ; read = read->next){
		read->information->is_used = FALSE;
	}
}

/*check the integrity of the sff file*/
void check_integrity(const char* const file_name, const sff_read* const reads)
{
	if(reads == NULL){
		return;
	}

	sff_read* read  = NULL;
	FILE* fp = ckopen(file_name,"r");
	common_header* common_header = ckallocz(sizeof(struct common_header_st));

	uint32_t num_reads = 0;
	const sff_read* iter;
	for(iter = reads; iter; iter = iter->next){
		num_reads++;
	}
	
	/*read the common_header*/
	read_common_header(fp, common_header);
	assert(common_header->num_reads == (uint32_t)num_reads);
	int flows_per_read = common_header->num_flows_per_read;

	uint32_t i = 0;
	for(i = 0; i < common_header->num_reads; i++){
		read = ckallocz(sizeof(sff_read));
		read->header = ckallocz(sizeof(read_header));
		read->data = ckallocz(sizeof(read_data));
		read->information = ckallocz(sizeof(ancilliary));

		read_sff_rinfo(&read, fp, flows_per_read);
		sff_read_free(&read);
		ckfree(read);	
	}

	fclose(fp);
	ckfree(common_header->flow_chars);
	ckfree(common_header->key_sequence);
	ckfree(common_header);
}


