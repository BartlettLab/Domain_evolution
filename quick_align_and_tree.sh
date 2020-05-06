#!/bin/bash

#BSUB -W 4:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=500]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 17                	# Where X is in the set {1..X}
#BSUB -J quick     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22
module load noisy/1.5.12
module load iqtree/1.6.3
module load R/3.6.1_packages/tidyverse/1.3.0 gcc/8.1.0

run=FLS2

mafft --thread 16 $run.hits_and_scaffold_seqs.fa > $run.hits_and_scaffold_align.fasta #quick version

echo "running noisy to clean up alignment..."
	noisy -s --noconstant $run.hits_and_scaffold_align.fasta
	#remove noisy log files
	rm *typ.eps
	rm *sta.gr
	rm *idx.txt



echo "Running IQtree in fast mode with auto find model of evo and inference of site rate evolution.."
	iqtree -s *out.fas -fast -nt 16 -m LG+F+R7 #quick version
	#iqtree -s *out.fas -bb 1000 -wsr -nt 16 -m LG+F+R7 #full version, infer site rates


echo "cleaning up.."
		#remove iqtree log files
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
	rm *.log
	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy





exit


#this will produce a .treefile tree, with no support values. Use this to extract the target clade.





