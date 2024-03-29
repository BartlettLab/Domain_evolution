#!/bin/bash

##BSUB -W 96:00             	# How much time does your job need (HH:MM)
##BSUB -R rusage[mem=500]	# How much memory
##BSUB -R span[hosts=1]		# Keep on one CPU cluster
##BSUB -n 17                	# Where X is in the set {1..X}
##BSUB -J refine         	# Job Name
##BSUB -o out.%J           	# Append to output log file
##BSUB -e err.%J           	# Append to error log file
##BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22
module load noisy/1.5.12
module load iqtree/1.6.3
module load R/3.6.1_packages/tidyverse/1.3.0 gcc/8.1.0

################################################################################
#       1. This script will collect the clade hits from initial search, remove
#			genes less than 85% of the target length, then build a tree and 
#			infer the evolution rate at each site using IQtree
#       2. It will then add the other clade XI scaffold genes to the alignment,
#           and remake the tree.
#       3.It will then collect the nucleotide sequences for the same genes
################################################################################
#  to run directly: bsub < /home/jm33a/domain_evolution/scripts_and_resources/1KP_followup_analyses.sh

export run=${1} #this is the argument passed in from directing script
#run=IKU2 #use this if need to just run one analysis, but make sure to un-comment the BSUBs above
echo "This will run follow-up analyses for $run clade"

scripts_dir=/home/jm33a/domain_evolution/scripts_and_resources
cd /home/jm33a/domain_evolution/clade_trees/$run/tree_from_hits

mkdir follow_up_analyses


######### 1. SITE RATES #########
echo "***************************************************************"
echo "Step 1. Infer site-by-site rate of evolution for $run clade"
echo "***************************************************************"

mkdir follow_up_analyses/site_rates/
	target_length=$(grep -A1 AT.*_$run $run.hits_and_scaffold_seqs.fa | tail -1 | wc -m) #get target seq length
	echo "$run sequence length is $target_length"
	min_length=$(awk -v target_length="${target_length}" -v percent=".85" 'BEGIN{print (target_length*percent)}') #get 85% of that length
echo "only genes longer than $min_length will be collected and used for analysis"
	fgrep -w --no-group-separator -A 1 -f *geneIDs.txt $run.hits_and_scaffold_seqs.fa > follow_up_analyses/site_rates/$run.clade_seqs.fa #search in results for clade genes, extract those
	cd follow_up_analyses/site_rates/
echo "removing genes less than $min_length residues (85% of full length $run gene which is $target_length residues)"	
	awk -v min_length="${min_length}" 'BEGIN {RS = ">" ; ORS = ""} length($2) >= min_length {print ">"$0}' $run.clade_seqs.fa > $run.full_length_clade_seqs.fa #collect only sequences that length or longer
echo "$(grep ">" $run.full_length_clade_seqs.fa | wc -l) full length $run sequences were retained (out of $(grep ">" $run.clade_seqs.fa | wc -l) total $run 1KP sequences detected from initial search)" #report on how many sequences were too short

echo "aligning full-length $run clade sequences"
	linsi --thread 16 $run.full_length_clade_seqs.fa > $run.full_length_clade_align.fasta 

echo "Running IQtree with auto find model of evo and inference of site rate evolution.."
iqtree -s $run.full_length_clade_align.fasta -bb 1000 -wsr -nt 16 -m JTT+R9 #full version, infer evo model and site rates

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

echo "Finished with $run site rate tree and rate analysis."



######### 2. FULL TREE WITH OUTGROUPS #########
echo "***************************************************************"
echo "Step 2. Infer full $run tree with outgroups from peptide sequences"
echo "***************************************************************"
    cd /home/jm33a/domain_evolution/clade_trees/$run/tree_from_hits
    mkdir follow_up_analyses/full_tree_with_outgroups/

echo"collecting sequences..."
    #remove scaffold genes from inputs list and store new input IDs file as the refined input gene IDs
    grep -v "\." *geneIDs.txt > follow_up_analyses/full_tree_with_outgroups/$run.1kp_clade_geneIDs.txt #scaffold genes all have a period (.), this line removes these genes
    #collect sequences of the 1KP genes in clade
    fgrep -w --no-group-separator -A 1 -f follow_up_analyses/full_tree_with_outgroups/$run.1kp_clade_geneIDs.txt $run.hits_and_scaffold_seqs.fa > follow_up_analyses/full_tree_with_outgroups/$run.clade_all_seqs.fa #search in results for clade genes, extract those
    #move to folder and add scaffold sequences
    cd follow_up_analyses/full_tree_with_outgroups/
    cat $scripts_dir/scaffold_seqs.fa >> $run.clade_all_seqs.fa #add the backbone scaffold sequences to tree

echo "aligning $run clade sequences"
	linsi --thread 16 $run.clade_all_seqs.fa > $run.clade_all_align.fasta 

echo "running noisy to clean up alignment..."
	noisy -s *align.fasta
	#remove noisy log files
	rm *typ.eps
	rm *sta.gr
	rm *idx.txt

echo "Running IQtree with auto find model of evo and inference of site rate evolution.."
    iqtree -s *out.fas -bb 1000 -nt 16 #full version, infer evo model and site rates

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

echo "Finished with $run full tree with outgroups"



######### 3. COLLECT NUCLEOTIDE SEQS #########
echo "***************************************************************"
echo "Step 3. Collect full-length nucleotide sequences of $run clade genes"
echo "***************************************************************"
cd /home/jm33a/domain_evolution/clade_trees/$run/tree_from_hits/
mkdir follow_up_analyses/nucleotide_seqs


grep ">" follow_up_analyses/site_rates/$run.full_length_clade_seqs.fa | grep -v "\." | cut -c 2- > follow_up_analyses/nucleotide_seqs/$run.1kp_clade_geneIDs.txt #grab full length 1KP gene IDs from site rates analysis


cd follow_up_analyses/nucleotide_seqs


species_to_scan=$(wc -l <  $scripts_dir/flowering_plants_species_IDs) #the number of flowering plant species handles from 1KT
count=1#counter for progress
while read species_ID
do
    echo "Processing species $species_ID ($count out of $species_to_scan)"
    directory=$(find /project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs -name "*$species_ID*" -type d) #find directory of the species
    nucleotide_database_path=$directory/*.dnas.out #this is now the path to species' nucleotide database
   
    test -f $nucleotide_database_path #test if file exists. Write 0 to $? if yes, 1 if no
    if [ $? -eq 0 ] #only perform the search when a database exists, otherwise report the error
    then
        active_genes=$(grep $species_ID $run.1kp_clade_geneIDs.txt) #return only the genes from clade that are in that species 
        #next step is to collect the sequnce as a fasta, but format it in the same style ([speciesID]_[geneNumber])
        for each_gene in $active_genes
        do
            gene_number=$(echo $each_gene | cut -d "_" -f2) #strip the species handle and underscore
            echo ">$each_gene" >> $run.clade_1KP_nucleotides_seqs.fa #write the full gene name to the seqs file in fasta format
            fgrep -w --no-group-separator -A 1 $gene_number $nucleotide_database_path | tail -n 1 >> $run.clade_1KP_nucleotides_seqs.fa #pull the sequence of that gene from the database and write to seqs file
        done
    else
        echo "$species_ID does not have a nucleotide database"
    fi
    count=`expr $count + 1` #add to count for progress report
done < $scripts_dir/flowering_plants_species_IDs #the input for flowering plant species handles from 1KT

echo "adding $run scaffold sequences"
    fgrep -A 1 -f /home/jm33a/domain_evolution/clade_trees/$run/*geneIDs.txt --no-group-separator $scripts_dir/scaffold_CDS_seqs.fa | awk '{print $1}' >> $run.clade_1KP_nucleotides_seqs.fa # add scafold seqs from clade
    fgrep -A 1 outgroup $scripts_dir/scaffold_CDS_seqs.fa  >> $run.clade_1KP_nucleotides_seqs.fa # add outgroup to polarize
    

echo "aligning $run clade sequences"
#    linsi --thread 16 $run.clade_1KP_nucleotides_seqs.fa > $run.clade_1KP_nucleotides_align.fasta #note this is not by codon, still don't know how to do this from command line

