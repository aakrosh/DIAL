//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: build_sff_index.c
//
//----------
//
// build_sff_index--
//	Support for indexing genomic sequences across several sff files, using a
//	hash-based techinique.
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable arg list stuff
#include "sff_index.h"			// sff index file stuff
#include "sff.h"				// aakrosh's sff file stuff
#include "utilities.h"

//----------
//
// global variables (private to this module)
//
//----------

#define maxSffFiles 255

int   numFilenames = 0;
char* inputFilenames[maxSffFiles];

//----------
//
// prototypes for private functions
//
//----------

int main (int argc, char** argv);

static int  string_to_int (const char* s);
static void suicide       (const char* message);
static void suicidef      (const char* format, ...);

//----------
//
// build_sff_index--
//  Main program
//
//----------

int main
   (int		argc,
	char**	argv)
	{
	int		argIx;
	char*	arg;
	uint		numBuckets, numSeqs;
	sffft*	ft;
	sffbt*	bt;
	sffst*	st;
	FILE*	fp;
	uint		fIx, seqIx;
	char*	seqName;
	u32		h, bucket;
	f31		filePos;

	// parse command line

	numBuckets = 0;

	for (argIx=1 ; argIx<argc ; argIx++)
		{
		arg = argv[argIx];
		if (strstr(arg,"--buckets=") == arg)
			{
			numBuckets = string_to_int (1+strchr(arg,'='));
			if (numBuckets < 1)
				suicidef ("gotta have at least one bucket (not %d)", numBuckets);
			}
		else if (strstr(arg,"--") == arg)
			suicidef ("unknown argument: \"%s\"", arg);
		else
			{
			if (numFilenames == maxSffFiles)
				suicidef ("too many input sff files (max is %d)", maxSffFiles);
			inputFilenames[numFilenames++] = arg;
			}
		}

	if (numBuckets < 1)
		suicide ("gotta tell me how many hash buckets to use (--buckets=<num>)");
	if (numFilenames < 1)
		suicide ("gotta give me at least one input sff file");

	// read the info for the index

	ft = sff_index_filename_table (numFilenames, inputFilenames);
	if (ft == NULL)
		suicidef ("can't create table for %d filenames", numFilenames);

	st = sff_index_sequence_table (0, 24);	// $$$ de-hardwire this
	if (st == NULL)
		suicide ("can't create empty sequence index table");

	numSeqs = 0;
	for (fIx=0 ; fIx<ft->num ; fIx++)
		{
		sff_read* read = NULL;
		common_header* common_header;
		uint32_t i;
		int flows_per_read;

		fp = fopen (ft->name[fIx], "rb");
		if (fp == NULL)
			suicidef ("failed to open %s", ft->name[fIx]);

		// read the common_header

		common_header = ckallocz(sizeof(struct common_header_st));
		read_common_header (fp, common_header);
		flows_per_read = common_header->num_flows_per_read;

		i = 0;
		for (i=0 ; i<common_header->num_reads ; i++)
			{
			// $$$ need to watch filePos for overflow
			filePos = ftell (fp);

			read              = ckallocz (sizeof(sff_read));
			read->header      = ckallocz (sizeof(read_header));
			read->data        = ckallocz (sizeof(read_data));
			read->information = ckallocz (sizeof(ancilliary));
			read_sff_rinfo (&read, fp, flows_per_read);

			seqName = read->header->name;
			if (strlen(seqName) != 14)
				suicidef ("sequence name is not to our liking: %s", seqName);

			numSeqs++;
			h = sff_index_hash (seqName) % numBuckets;
			st = sff_index_add_sequence (st, h, seqName, fIx, filePos);
			if (st == NULL)
				suicidef ("unable to add sequence #%d to sequence index table", numSeqs);

			//fprintf (stderr, "%08X %s %03d %08X\n", h, seqName, fIx, filePos);

			sff_read_free (&read);
			ckfree (read);
			}

		ckfree (common_header->flow_chars);
		ckfree (common_header->key_sequence);
		ckfree (common_header);

		fclose (fp);
		}

	sff_index_sort_sequences (st);

	// write the index
	// $$$ need to watch filePos for overflow
	// $$$ need check for errors

	bt = sff_index_bucket_table (numBuckets);
	if (bt == NULL)
		suicidef ("can't create table for %d buckets", numBuckets);

	filePos =  sff_index_write_header (stdout, ft, bt, st);
	filePos += sff_index_write_filenames (stdout, ft);

	bucket = 0;
	for (seqIx=0 ; seqIx<st->num ; seqIx++)
		{
		h = sff_index_sequence_record_h (st, seqIx);
		while (bucket < h) // fill entries for intervening (empty) buckets
			bt->bucket[bucket++] = f31_empty + filePos;
		if (h == bucket)   // fill entry for first sequence in a bucket
			bt->bucket[bucket++] = filePos;
		// $$$ need to watch filePos for overflow
		filePos += sff_index_write_sequence (stdout, st, seqIx);
		}

	while (bucket <= numBuckets) // fill entries for final (empty) buckets
		bt->bucket[bucket++] = f31_empty + filePos;

	filePos += sff_index_write_sequences_pad (stdout, st);
	filePos += sff_index_write_buckets (stdout, bt);

	free (ft);
	free (bt);
	free (st);

	return EXIT_SUCCESS;	// $$$ return proper error code in error cases
	}

//----------
//
// string_to_int--
//	Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in fatality.
//
//----------

static int string_to_int
   (const char*	s)
	{
	int			v;
	char		extra;

	if (sscanf (s, "%d%c", &v, &extra) != 1)
		suicidef ("\"%s\" is not an integer", s);

	return v;
	}

//----------
//
// suicide, suicidef--
//	Cause program fatality, after pushing a message out to the user.
//
//----------
//
// Arguments for suicide():
//	const char*	message:	The message to write to stderr before death.  This
//							.. may be NULL.
//
// Arguments for suicidef():
//	const char*	format:		A format string, as per printf.  This may be NULL.
//	...:					(same as for printf)
//
// Returns:
//	(nothing;  it does not return).
//
//----------

static void suicide
   (const char*	message)
	{
	if (message == NULL) suicidef (NULL, NULL);
	                else suicidef ("%s", message);
	}

static void suicidef
   (const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);

	fflush  (stdout);
	fprintf (stderr, "FAILURE: ");
	if (format != NULL)
		{
		vfprintf (stderr, format, args);
		fprintf  (stderr, "\n");
		}

	va_end (args);

	exit (EXIT_FAILURE);
	}

