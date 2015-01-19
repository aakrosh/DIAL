/* Simple module that makes the first cluster file from a file of reads */

#include "utilities.h"
#include "sequences.h"

#define USE "make_template reads.fa individual"

static void make_template(const char* const seqfile, const char* const name)
{
	sequence* sp;
	if((sp = read_fasta_sequence(seqfile)) == NULL){
		fatalf("error in reading the sequences from %s", seqfile);
	}

	while(sp){
		printf(">%s\n", sp->header);

		printf("%s\n", sp->sequence);
		printf("%s\t%s\n\n", sp->header, name);

		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);
}

int main(int argc, char** argv)
{
	if(argc != 3){
		fatal(USE);
	}

	make_template(argv[1], argv[2]);

	return EXIT_SUCCESS;
}
