DIAL
====
De novo Identification of Alleles

## REQUIREMENTS
DIAL should work on any standard 64 bit Linux environment with gcc and python.
However we have tested the compilation on the following platforms (output of
gcc --version and uname -a):

- gcc (GCC) 4.1.2 20080704 (Red Hat 4.1.2-46)
  Linux 2.6.18-164.el5 ##1 SMP x86\_64 GNU/Linux

- gcc (Debian 4.3.2-1.1) 4.3.2
  Linux 2.6.18-4-amd64 ##1 SMP x86\_64 GNU/Linux

- powerpc-apple-darwin8-gcc-4.0.1 (GCC) 4.0.1 (Apple Computer, Inc. build 5367)
  Darwin 8.11.0 Darwin Kernel Version 8.11.0: root:xnu-792.24.17~1/RELEASE\_PPC Power Macintosh powerpc


DIAL uses LASTZ (http://bx.psu.edu/miller_lab) for aligning the sequences to
each other. Please install LASTZ and add the LASTZ binary in your
executable/binary search path before installing DIAL. 

## SUMMARY
DIAL (De novo Identification of Alleles) is a collection of programs to automate
the discovery of alleles for a species where we lack a reference sequence. The
SNPs/alleles are specifically selected for a low error rate in genotyping 
assays. See the details of the methods here: http://www.biomedcentral.com/1471-2105/11/130/ 

## INSTALLATION
The programs can be compiled by following the following recipe (see INSTALL for
greater detail):

'''
  % configure --prefix=/usr/local (or whatever installation path you prefer)
  % make
  % make install
'''

This complies all the components of the pipeline and puts the binaries in the
folder $prefix/bin.  For more in depth instructions, consult the INSTALL file.

Please add the $prefix/bin folder to your executable/binary search path to
complete the installation.

## DESCRIPTION
The pipeline is run using the 'DIAL' binary. The format of commands issued to
DIAL is as follows:

DIAL ${action} ${arguments}

where ${action} is one of the following:
a) create
b) add
c) update

We now describe each of these actions and the related arguments:

a) create- This action is used to create and establish a folder structure for a
new project. The format of the command is as follows: 
		
```
DIAL create project_name 454/illumina 
```
where project\_name is the unique identifier and the name of the parent folder
for the project. The folder is created in the directory where the command is
issued. Alternatively, one can specify a full path to an alternative
location. The last argument is '454' or 'illumina' depending on the sequencing
technology used to generate the reads used in this project. It also creates a
'status.txt' file in the project folder, which should be edited, to make it
suitable for the specific project. The default 'status.txt' file looks as
follows for a 454 dataset:
 
```
Dataset: 454
454 Binaries: /usr/local/rig/bin
Expected size of the genome: 3000000000
Runs added: 0
Bases added: 0
Bases used: 0
```

The expected size of the genome must be changed to reflect the expected genome
size in the current project.  The location corresponding to "454 Binaries"
should point to directory with 'sffinfo', 'Newbler', and other binaries that are
shipped with the 454 sequencer. For a illumina dataset, the location should
point to the directory with the velvetg, velveth
(http://www.ebi.ac.uk/~zerbino/velvet/) binaries.

b) add- This action is used to add a single run/lane of sequences to the
		project. The format of the command is as follows:

```		
DIAL add project_name sff_file/lane_sequence name_individual [-transcript]
```

where project\_name is the unique identifier and the name of the parent folder
for the project. It should exist in the directory that the command is issued in.
Alternatively, one can list the full path for the directory. We use SFF files
for a 454 dataset and the complete path of the SFF file must be supplied as as
input. "fastq" files for a lane must be used if the dataset under consideration
has been generated using Illumina read technology. 'name\_individual' refers to a
unique identifier/name of the sample this particular lane/run belongs to. The
"-transcript" switch should be added if the input is from a transcriptome
dataset.

c) update- This action is used to call SNPs from the runs in the project. The format of the command is as follows:

```		
DIAL update project_name
```
where project\_name is the same identifier/folder name used in the earlier
commands. The filtered subset of high-confidence calls can be seen in the file
project\_names/alleles/snps.txt. The file shows the assembled contigs, with the
position of the variant alleles. It also shows the number of reads supporting
the allele and the quality value of the reads at that position.

## TEST-DATASET (An example for 454 reads and another one for Illumina reads)
A sample test-dataset is included with this archive to demonstrate the use of
DIAL. Please read and run the run454\_test.bash and/or runillumina\_test.bash in 
the folder 'test\_data' to see a sample run of this pipeline.
 
The 454 test should result in 8 SNPs (when Newbler 2.0 is used), whereas the 
Illumina test case should result in 18 SNPs

## Application Note for Illumina Sequences
The number of base-pairs/lane has increased significantly in the past year. DIAL
was designed to work with significantly fewer sequences/lane as compared to that
can be generated today. For better results, users should split their Illumina
datasets into multiple files and then "add" them to the project one after the 
other. For example if they have a fastq sequence file s\_1\_1\_sequence.txt, then
they should do the following:

```
# split the sequences into multiple files with names xaaa,xaab,xaac and so
# on. Each of those files will have 500,000 fastq sequences.
split -a 3 -l 2000000 s_1_1_sequence.txt x

# now create the project "foo", and then modify the foo/status.txt file
DIAL create foo illumina

# add the files (xaaa,xaab....) one after the other
DIAL add foo /home/bar/xaaa individual1
DIAL add foo /home/bar/xaab individual1
DIAL add foo /home/bar/xaac individual1
.
.
.

# finally call the required SNPs
DIAL update foo
```

## CHANGE-HISTORY
Jun 04, 2011:
* Fixed a bug in calculate\_maxhits which led to an underestimation of maximum
  allowed hits for high coverages.
* Fixed a bug in assemble\_illumina, which led to incorrect format in
  alleles/report.txt.
* Added patches to support compilation on Solaris. Thanks to Nathan Weeks (Iowa
  State University) for the patch.

Jul 31, 2013:
* Made changes to convertFastqFasta to handle new Illumina reads which have more
  than one space delimited token in their name. If that is the case,
  convertFastqFasta now converts the name to the old format.
* Made changes to "DIAL" to handle the new read names
* Made changes to update\_clusters so that a "N" in one of the reads is no longer
  treated as a difference.
* Made changes to remove\_clones to remove the deprecation message for md5
* Added a warning message to assemble.c to warn the user if the "flows per
  cycle" in the input is not yet supported.
