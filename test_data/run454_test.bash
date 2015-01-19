#!/usr/bin/env bash

set -e 

# This script walks you through the steps to run DIAL on a 454 dataset. The
# script fetches three SFF files from the Watson dataset using wget, and then
# runs the whole pipeline on them. This uses wget to fetch the SFF files; if
# that fails the files should be downloaded from
# http://www.bx.psu.edu/~ratan/data before this script is run.
#
# This script assumes that 'lastz' (http://bx.psu.edu/miller_lab/) is installed
# and the binary is in your $PATH. It also assumes that DIAL is complied and the
# binaries are in the './../bin' directory (relative to it). 
# Please type
#	./run454_test.bash step1
# This will create the project and folder heirarchy. Then please go and update
# the watson/status.txt file to point it to the folder which has Newbler,
# sffinfo and other binaries that are shipped with the 454 sequencing system.
# Then please type
#	./run454_test.bash step2
# This will add the first SFF file to the project. The file 'watson/status.txt'
# should reflect the addition of one run to the project. The fasta file and 
# quality file for the run are in watson/reads and soft links to the 
# SFF file are found in watson/sff. The folder watson/clusters holds the
# clusters that were formed when this run was added.
# Then proceed to type:
#	./run454_test.bash step3
# Again please notice the changes to the various files in the project.
# Then add the last run to the project by typing:
#	./run454_test.bash step4
# 
# At this point, we have included all the runs in the project. Now to call 
# SNPs and select a high-confidence subset from them, please type the following:
#	./run454_test.bash step5
# The resultant high confidence SNPs are now listed in watson/alleles/snps.txt
# in the following format:
#>Cluster_15
#TtCTTTTtGCGCTTCTTCAGAaGTTGACTCTTTTGGCCCTTTGGTCTTCTATACaCATTT
#TAGAAATGCTTTGTtGAGGACTAAGAGGAaTGCTAAGaTTTTGATAGGAATTTCATTGAA
#TTTTGAGTATATTcgCATGCTACAATGGTTAGTGCTTTATACATGAAAATAaTATATCCC
#TTCCTCTTTTCCTAGTATCATGAGATGTTTGTAGGCAGACATGAATaTTGAGTTGTATCA
#AATGTGGTTTt
#133     134     watson:C=2(28,24) watson:G=2(27,27)
#
#listing the assembled contig, the variant position, the alleles, number of
#reads supporting that allele, and the quality value of the reads at the variant
#position

C=./../bin

#lets get the full path for this directory
current=$(pwd)

if [ $# -ne 1 ]; then
    echo "./run454_test.bash step1/step2/step3/step4/step5"
    exit 1
fi

if [ ! -f SRR000375.sff ]; then
	wget http://www.bx.psu.edu/~ratan/data/SRR000375.sff
fi
if [ ! -f SRR000376.sff ]; then
	wget http://www.bx.psu.edu/~ratan/data/SRR000376.sff
fi
if [ ! -f SRR000377.sff ]; then
	wget http://www.bx.psu.edu/~ratan/data/SRR000377.sff
fi

if [ "$1" == "step1" ]; then
	#create the project 
	$C/DIAL create watson 454
fi

if [ "$1" == "step2" ]; then
	#add the first SFF file to the project
	$C/DIAL add watson ${current}/SRR000375.sff jwatson
fi

if [ "$1" == "step3" ]; then
	#add the second SFF file to the project
	$C/DIAL add watson ${current}/SRR000376.sff jwatson
fi

if [ "$1" == "step4" ]; then
	#add the last SFF file to the project
	$C/DIAL add watson ${current}/SRR000377.sff jwatson
fi

if [ "$1" == "step5" ]; then
	#call SNPs
	$C/DIAL update watson
fi
