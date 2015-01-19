/* filter:  impose constraints to Aakrosh's Newbler output
*
* Command line:

filter aakrosh.out [options]

* where aakrosh.out is the output of the assemble_exons program, and the
* options can have the following forms (integer constants can vary).
*	most-diffs=2	most differences found in the contig
*	most-reads=4	most reads in the contig
*	min-side=50	min distance of SNP from either end of the contig
*	min-score=20	lowest score for an allele to be counted
*	most-units=2	most difference "units", e.g. "cedric:C=2(24,21)";
*			"most-units=2 two-names" finds homozygous SNPs.
*	most-depth=4	most reads at the SNP position
*	homopolymer=5,13 disallow a homopolymer of lengt >= 5 in a 13-bp window
*			centered on the SNP
*	two-names	require that two different names appear in the SNP
*	contig		print Newbler's "contig" line in the output
*	good		mark good SNPs as "good"
*
*  The input for this program consists of blocks looking like this:
>EITVI5O01EZ9T6
>contig00001  length=259  numreads=6
AGCTCAaCAACATTGATCAGTTAGCAGTCCGGTTCAGTACAACATTGATCAGGTAGCAAG
CCGGTTtGGCAGTGCCAGTCGGAGACACAGATGGGGACAATCgAAGAAGATGCAGGTCAG
GACAAGGATGCAGACATGTACCTGTCACTACCATGGACTTCAGCAGCTtCTTAGCAGGTT
TATGATTAAAGAAGAAAGACATCGTTGATATGTAAGCATGTTACAaTTAATTTGCTGTTT
AGTTATTAGCAAACTACAT
66      67      celeste:C=3(25,28,28) samuel:T=3(16,13,13) 
102     103     celeste:A=3(17,15,26) samuel:G=3(26,28,26)

* (ending with a blank line). The contig is followed by zero or more
* differences. If there are multiple contigs associated with the same read,
* they are listed with no separating blank line.
*
* In the output, the letter at each differing position is placed in square
* brackets. In addition, cases violating any of the following conditions
* are not printed.
* 1. There can only be one contig for the cluster.
* 2. The number of reads in the contigs is at most MOST_READS
* 3. The number of differences is (a) at least 1 and (b) at most MOST_DIFFS.
* 4. One of the differences much be "good" in thes sense of meeting the
*    conditions a-g:
* 	(a) there are at most 2 alleles
* 	(b) each allele is supported by at least 2 reads
* 	(c) there are at least MIN_SIDE positions on each side of the SNP,
* 	(d) alleles with score less than MIN_SCORE are ignored,
* 	(e) the number of "name:nuc=many(score1,score2,..)" units is at most
*	    MOST_UNITS,
*	(f) the number of reads at the SNP is at most MOST_DEPTH
*	(g) there is no homopolymer of length HOMO_LEN in a window of
	    HOMO_WINDOW bp centered at the SNP,
* 	(h) optionally require that at least 2 individuals be represented.
*/

#include <ctype.h>

#include "utilities.h"

// most reads in the contigs
#define MOST_READS INT_MAX

// most SNPs predicted by Newbler in a contig, default
#define MOST_DIFFS 3

// minimum number of bp on each side of the SNP, default
#define MIN_SIDE 40

// minimum total score for each allele at a "good" difference, default
#define MIN_SCORE 20

// maximum number of "name:nuc=many(score)" units, default
#define MOST_UNITS 4

// maximum number of reads at the SNP
#define MOST_DEPTH INT_MAX

// homopolymer length and window diameter for rejecting a SNP, defaults
#define HOMO_LEN 0
#define HOMO_WINDOW 0

int two_names = 0, contig = 0;

char *options = "[most-reads=?] [most-diffs=?] [min-side=?]\n[min-score=?] [most-units=?] [homopolymer=?,?] [two-names] [contigs] [good]";

// longest line of input
#define LONGEST 1000
// most lines in an contig block
#define MOST 1000
// longest contig
#define CONT_MAX 10000

char x[MOST][LONGEST], y[CONT_MAX];
int is_good[MOST];

int main(int argc, char **argv) {
	
	argv0="filter";
	FILE *fp;
	int nline, i, j, k, n, kk, p, pos, blocks = 0, no_diffs = 0, units,
	  no_diffs2 = 0, excess_diffs = 0, good = 0, contig2 = 0, nhomo = 0,
	  many[4], most_reads = MOST_READS, min_score = MIN_SCORE,
	  most_diffs = MOST_DIFFS, most_units = MOST_UNITS, min_side = MIN_SIDE,
	  most_depth=MOST_DEPTH, homo_len = HOMO_LEN, homo_window = HOMO_WINDOW,
	  names, length, radius = -1, excess_reads = 0, do_good = 0;
	char c, *s, *t, *g, *h;

	if (argc < 2)
		fatalf("args : Aakrosh-style-snps %s", options);
	for (i = 2; i < argc; ++i) {
		if (strncmp(s = argv[i], "most-reads=", 11) == 0)
			most_reads = atoi(s+11);
		else if (strncmp(s, "most-diffs=", 11) == 0)
			most_diffs = atoi(s+11);
		else if (strncmp(s, "min-side=", 9) == 0)
			min_side = atoi(s+9);
		else if (strncmp(s, "min-score=", 10) == 0)
			min_score = atoi(s+10);
		else if (strncmp(s, "most-units=", 11) == 0)
			most_units = atoi(s+11);
		else if (strncmp(s, "most-depth=", 11) == 0)
			most_depth = atoi(s+11);
		else if (strncmp(s, "homopolymer=", 12) == 0) {
			homo_len = atoi(s+12);
			if ((t = strchr(s, ',')) == NULL)
				fatalf("no comma: %s", s);
			homo_window = atoi(t+1);
			radius = homo_window/2;
			if (radius + radius == homo_window)
				fatal("homopolymer window size must be odd");
		} else if (strcmp(s, "two-names") == 0)
			two_names = 1;
		else if (strcmp(s, "contig") == 0)
			contig = 1;
		else if (strcmp(s, "good") == 0)
			do_good = 1;
		else
			fatalf("improper option: %s", s);
	}
	fprintf(stderr,
	  "most-reads = %d, most-diffs = %d, min-side = %d, min-score = %d\n",
	  most_reads, most_diffs, min_side, min_score);
	fprintf(stderr,
	  "most_units = %d,  homo-len = %d, most-depth = %d, homo-win = %d, two-names = %d\n",
	  most_units, homo_len, most_depth, homo_window, two_names);
	fp = ckopen(argv[1], "r");
	while (fgets(x[nline = 0], LONGEST, fp)) {
		if (x[0][0] != '>')
			fatalf("block starts with %s", x[0]);
		++blocks;
		// get the next block
		while (nline < MOST && fgets(x[++nline], LONGEST, fp) &&
		   x[nline][0] != '\n')
			is_good[nline] = 0;; 
		if (nline == MOST)
			fatalf("too many lines: %s", x[0]);
		for (i = 2; i < nline && x[i][0] != '>'; ++i)
			;
		// condition 1
		if (i < nline) { // multiple contigs
			++contig2;
			continue;
		}
		if (x[1][0] == '\n')
			continue;	// no contigs

		// condition 2
		// limit the number of reads in the contig
		if ((s = strstr(x[1], "numreads=")) == NULL)
			fatalf("numreads: %s", x[0]);
		if (atoi(s+9) > most_reads) {
			++excess_reads;
			continue;
		}

		// find the first line of the differences
		for (k = nline; isdigit(x[k-1][0]); --k)
			;
		// condition 3a
		if (k == nline) {
			++no_diffs;
			continue;
		}
		// condition 3b
		// limit the number of differences
		if (nline - k > most_diffs) {
			++excess_diffs;
			continue;
		}

		if ((s = strstr(x[1], "length=")) == NULL)
			fatalf("length: %s", x[0]);
		length = atoi(s+7);

		// collect the contig
		y[0] = '\0';
		s = y;
		for (i = 2; i < k; ++i) {
			strcpy(s, x[i]);
			while (*s && *s != '\n')
				++s;
			if (*s != '\n')
				fatalf("line size is inadequate");
			*s = '\0';
		}
		if (length != (s-y))
			fatalf("%s   lengths %d %d\n", x[0], length, s-y);

		// check for a good diff
		// 132 133 cedric:C=2(20,28) spirit:T=2(19,30) spirit:C=1(33)
		for (i = k; i < nline; ++i) {
			// keep looking for a good SNP
			pos = atoi(x[i]);
			// condition 4c
			if (pos < min_side || pos > length-min_side)
				continue;
			for (j = 0; j < 4; ++j)
				many[j] = 0;
			g = NULL;
			// count units, names and number of adequate-scoring
			// reads per allele
			for (names = units = 0, s = x[i];
			     (s = strchr(s,'=')) != NULL; ++units) {
			     	if (g == NULL) {
					g = s-2;
					if (*g != ':')
						fatalf("colon: %s", x[i]);
					while (g > x[i] && isalnum(g[-1]))
						--g;
					names = 1;
				} else {
					h = s-2;
					if (*h != ':')
						fatalf("colon: %s", x[i]);
					while (h > g && isalnum(h[-1]))
						--h;
					while (*g == *h && *g != ':') {
						++g;
						++h;
					}
					if (*g != ':' || *h != ':')
						names = 2;
				}
				if ((t = strchr("ACGT",s[-1])) == NULL)
					fatalf("nucleotide: %c", s[-1]);
				kk = t - "ACGT";
				for ( ; ; ) {
					while (!strchr("(),", *s))
						++s;
					if (*s == ')')
						break;
			// condition 4d
					if (atoi(++s) >= min_score)
						++many[kk];
				}
			}
/*
fprintf(stderr, "%s", x[i]);
fprintf(stderr, " A:%d C:%d G:%d T:%d\n", many[0], many[1], many[2], many[3]);
*/
			// condition 4h
			if (two_names && names < 2)
				continue;
			// condition 4e
			if (units > most_units)
				continue;
			// how many alleles are represented?
			// how many alleles have at least 2 reads
			for (k = j = 0; j < 4; ++j)
				if (many[j] > 0)
					++k;
			// conditing 4f
			if (many[0]+many[1]+many[2]+many[3] > most_depth)
				continue;
			
			// condition 4a
			if (k != 2)
				continue;
			// neither has only one example
			for (j = 0; j < 4; ++j)
				if (many[j] == 1)
					break;
			// condition 4b
			if (j < 4)
				continue;
			// check for nearby homopolymer
			for (n = j = pos-radius; j <= pos+radius; ++j) {
				if (toupper(y[n]) != toupper(y[j]))
					n = j;
				if (j-n == homo_len - 1)
					break;
			}
			// condition 4g
			if (j <= (pos+radius)) {
			//fprintf(stderr, "homopolymer %s%s\n\n", x[i], y);
				++nhomo;
				continue;
			}
			is_good[i] = 1;
		}
		for (i = k; i < nline && is_good[i] == 0; ++i)
			;
		if (i == nline) {
			++no_diffs2;
			continue;
		}

		++good;
		pos = atoi(x[kk = k]);
		printf("%s", x[0]);
		// skip "numreads" line unless the "contig" flag is given
		if (contig)
			printf("%s", x[1]);
		for (p = 0, i = 2; i < kk; ++i) {
			for (j = 0; (c = x[i][j]) != '\n'; ++j, ++p)
				if (p == pos) {
					printf("[%c]", c);
					pos = atoi(x[++k]);
				} else
					putchar(c);
				putchar('\n');
		}
		for (k = kk; k <= nline; ++k) {
			if ((s = strchr(x[k], '\n')) != NULL)
				*s = '\0';
			printf("%s", x[k]);
			if (do_good && is_good[k])
				printf(" good");
			putchar('\n');
		}
	}
	fprintf(stderr,
	  "%d blocks, %d good, %d no diffs, %d no good diffs, %d excess diffs,\n%d excess reads, %d multiple contigs, %d homopolymer\n",
	  blocks, good, no_diffs, no_diffs2, excess_diffs, excess_reads, contig2, nhomo);
	return 0;
}
