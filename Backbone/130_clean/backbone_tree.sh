#!/bin/bash

#BSUB -W 4:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=1000]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 12                	# Where X is in the set {1..X}
#BSUB -J search     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}


module load MAFFT/7.313
module load noisy/1.5.12
module load iqtree/1.6.3


#construct database
cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
database=combined_onelines.fa

echo "collecting sequences..."
input_file=cannonical_IDs.txt
	grep -w -A 1 -f $input_file --no-group-separator $database | awk '{print $1}' > $input_file.seqs.fa

	echo "Aligning..."
	mafft $input_file.seqs.fa >  $input_file.align.fasta
echo "Filtering with noisy"
noisy -s $input_file.align.fasta #clean with noisy

echo "launching IQtree"
	iqtree -s *out.fas -bb 1000 -nt 12 -m LG+F+R7


rm $database

#remove iqtree log files
rm *splits.nex
rm *ckp.gz
rm *bionj
rm *treefile
rm *mldist
rm *.log
rm *.iqtree
rm *model.gz
rm *uniqueseq.phy
rm *idx.txt
rm *sta.gr
rm *typ.eps

