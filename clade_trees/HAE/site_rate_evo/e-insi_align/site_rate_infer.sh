#!/bin/bash

#BSUB -W 48:00             	# How much time does your job need (HH:MM)
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -R rusage[mem=1000]	# How much memory
#BSUB -n 16                  # Where X is in the set {1..X}
#BSUB -J flanks          	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file


module load noisy/1.5.12 
module load iqtree/1.6.3
module load MAFFT/7.313




#linsi --thread -10 *seqs.fa > linsi.align.fasta #linsi was found to be be the best of all mafft algorithms





#perform the site evolution rate test 
#iqtree -s *align.fasta -t *.contree -wsr -nt 14 # if already have a tree can specify with the -t, otherwise it will infer a new one

iqtree -s *align.fasta -wsr -nt 14
