#!/usr/bin/env bash

set -e 

# This script walks you through the steps to run DIAL on a illumina dataset. The
# script fetches 4 simulated fastq files, and then runs the whole pipeline on 
# them.
#
# This script assumes that 'lastz' (http://bx.psu.edu/miller_lab/) is installed
# and the binary is in your $PATH. It also assumes that DIAL is complied and the
# binaries are in the './../bin' directory (relative to it). 

# Please type
#	./runillumina_test.bash step1
# This will create the project and folder heirarchy. Then please go and update
# the watterson/status.txt file to point it to the folder which has velveth and
# velvetg; both are part of the Velvet assembler 
# (http://www.ebi.ac.uk/~zerbino/velvet/). Also please update the expected size
# of the genome to be 8000, as the reads in this experiment were generated from
# a synthetic reference of that size.

# Then please type
#	./runillumina_test.bash step2
# This will add the first fastq file to the project. The file 
# 'watterson/status.txt' should reflect the addition of one run to the project. 
# The fasta file and quality file for the run are in watterson/reads. The 
# folder watterson/clusters holds the clusters that were formed when this run 
# was added.

# Then proceed to type:
#	./runillumina_test.bash step3
# Again please notice the changes to the various files in the project.
# Then add the last two runs to the project by typing:
#	./runillumina_test.bash step4
#   ./runillumina_test.bash step5
# 
# At this point, we have included all the runs in the project. Now to call 
# SNPs and select a high-confidence subset from them, please type the following:
#	./runillumina_test.bash step6
# The resultant high confidence SNPs are now listed in watterson/alleles/snps.txt
# in the following format:
#>Cluster_15
#TtCTTTTtGCGCTTCTTCAGAaGTTGACTCTTTTGGCCCTTTGGTCTTCTATACaCATTT
#TAGAAATGCTTTGTtGAGGACTAAGAGGAaTGCTAAGaTTTTGATAGGAATTTCATTGAA
#TTTTGAGTATATTcgCATGCTACAATGGTTAGTGCTTTATACATGAAAATAaTATATCCC
#TTCCTCTTTTCCTAGTATCATGAGATGTTTGTAGGCAGACATGAATaTTGAGTTGTATCA
#AATGTGGTTTt
#133     134     calvin:C=2(28,24) hobbes:G=2(27,27)
#
#listing the assembled contig, the variant position, the alleles, number of
#reads supporting that allele, and the quality value of the reads at the variant
#position

C=./../bin

#lets get the full path for this directory
current=$(pwd)

if [ $# -ne 1 ]; then 
    echo "runillumina_test.bash step1/step2/step3/step4/step5/step6"
    exit 1
fi

if [ ! -f reads11.fq ]; then
	wget http://www.bx.psu.edu/~ratan/data/reads11.fq
fi
if [ ! -f reads12.fq ]; then
	wget http://www.bx.psu.edu/~ratan/data/reads12.fq
fi
if [ ! -f reads21.fq ]; then
	wget http://www.bx.psu.edu/~ratan/data/reads21.fq
fi
if [ ! -f reads22.fq ]; then 
	wget http://www.bx.psu.edu/~ratan/data/reads22.fq
fi

if [ "$1" == "step1" ]; then
	#create the project 
	$C/DIAL create watterson illumina
fi

if [ "$1" == "step2" ]; then
	#add the first fastq file to the project
	$C/DIAL add watterson ${current}/reads11.fq calvin
fi

if [ "$1" == "step3" ]; then
	#add the second fastq file to the project
	$C/DIAL add watterson ${current}/reads12.fq calvin
fi

if [ "$1" == "step4" ]; then
	#add the third fastq file to the project
	$C/DIAL add watterson ${current}/reads21.fq hobbes
fi

if [ "$1" == "step5" ]; then
    #add the last fastq file to the project
    $C/DIAL add watterson ${current}/reads22.fq hobbes
fi

if [ "$1" == "step6" ]; then
	#call SNPs
	$C/DIAL update watterson
fi
