#!/bin/bash

#BSUB -W 96:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=500]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 17                	# Where X is in the set {1..X}
#BSUB -J search     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22
module load noisy/1.5.12
module load iqtree/1.6.3
module load R/3.6.1_packages/tidyverse/1.3.0 gcc/8.1.0

################################################################################
#       This script will take a file of input search genes [*geneIDs.txt]
#       and perform an HMM search into all angiosperm 1KT databases.
#       It will collect the top 2 hits from each species, then remove hits that
#       do not contain both an LRR and an RLK domain.
#       A custom R script sort_full_LRR_RLKs_olnly.R is used for filtering,
#       then the top hits are used to infer a gene tree. 
#       The gene tree can then be used to find the clade and extract only genes
#       that fall within the ratiation after divergence with the Amborella copy.
################################################################################
#run as bsub < /home/jm33a/domain_evolution/scripts_and_resources/1KP_initial_search_and_tree.sh

scripts_dir=/home/jm33a/domain_evolution/scripts_and_resources

run=${PWD##*/} #name of current directory without path, should be exactly the same as the gene name in the gene ID
    #set up target gene length and cutoffs
	target_length=$(grep -A1 AT.*_$run $scripts_dir/scaffold_seqs.fa | tail -1 | wc -m) #get target seq length
	echo "$run sequence length is $target_length"
	min_length=$(awk -v target_length="${target_length}" -v percent=".85" 'BEGIN{print (target_length*percent)}' | awk '{print int($1+0.5)}') #get 85% of that length	
	max_length=$(awk -v target_length="${target_length}" -v percent="1.15" 'BEGIN{print (target_length*percent)}'| awk '{print int($1+0.5)}') #get 115% of that length
echo "only genes of length between $min_length and $max_length will be collected and used for analysis"



##this section will take a list of genes, make an HMM, search all 1kt and n best hit collect sequences, then scan them for domain content
cd /home/jm33a/domain_evolution/clade_trees/$run/

<<MULTILINE_COMMENT
echo "constructing database for input sequences, collecting and aligning inputs..."
    cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa #source databases for input genes
    grep -w -A 1 -f *geneIDs.txt --no-group-separator combined_onelines.fa | awk '{print $1}' > $run.input_seqs.fa #collect input clade seqs
    rm combined_onelines.fa
    linsi $run.input_seqs.fa > $run.input_align.fasta
MULTILINE_COMMENT

echo "make $run hmm from input alignment..."
	hmmbuild -o hmmout.txt $run.model.hmm $run.input_align.fasta
	rm hmmout.txt


species_to_scan=$(wc -l <  $scripts_dir/flowering_plants_species_IDs) #the number of flowering plant species handles from 1KT
#species_ID=LSKK
echo "collect the top 20 hits for each of the $species_to_scan species listed in flowering_plants_species_IDs"
count=1 #counter for progress
while read species_ID
do
    echo "Processing species $species_ID ($count out of $species_to_scan)"
    directory=$(find /project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs -name "*$species_ID*" -type d) #find directory of the species
    protein_database_path=$directory/*.prots.out #this is now the path to species' protein database
   
    test -f $protein_database_path #test if file exists. Write 0 to $? if yes, 1 if no
    if [ $? -eq 0 ] #only perform the search when a database exists, otherwise report the error
    then #this bock will collect only sequences that meet length criteria and add the top 2 hits to the top_hits list
        hmmsearch -o hmmout.txt --noali --tblout table.output.txt $run.model.hmm $protein_database_path #scan all sequences from a species with the HMM model
        for hit_number in {4..23} #top 20 hits pre-length thresholding
        do  #check sequence length
            hit=$(sed -n "${hit_number}p" < table.output.txt | awk '{print $1}') #gene ID from hit number line
            hit_length=$(grep -w -A 1 $hit $protein_database_path | tail -n 1 | wc -m)
            if [ $hit_length -gt $min_length ]
            then
                if [ $hit_length -lt $max_length ]
                then #add to the temp list if both these criteria are met
                    hit_full_ID=">${species_ID}_${hit}" #concatenate the species ID and the geneID with an underscore between, using fasta format >
                    echo $hit_full_ID >> temp_hits.fa #add gene full gene ID to temp doc
                    fgrep -w -A 1 ">${hit}" $protein_database_path | tail -n 1  >> temp_hits.fa #add gene's sequence to running list
                fi
            fi
        done
 
        test -f temp_hits.fa #if any hits were found, add top two hits to the top hits sequences list and reset the temp list
        if [ $? -eq 0 ]
        then
            cat temp_hits.fa | head -4 >> $run.top_hits.seqs.fa
            rm temp_hits.fa #reset the temp list
        fi
    fi
    count=`expr $count + 1` #add to count for progress report
done < $scripts_dir/flowering_plants_species_IDs #the input for flowering plant species handles from 1KT

#clean up a bit before proceeding
rm hmmout.txt
rm table.output.txt

echo "removing hits with sequences that were not found"
    grep -v "#" $run.top_hits.seqs.fa > top_hits_clean.fa
    rm $run.top_hits.seqs.fa
    cat top_hits_clean.fa > $run.top_hits.seqs.fa
    rm top_hits_clean.fa


#this section to remove any genes for which both an LRR and RLK domain are not found
echo "Building Pfam domain table from genes" 	##in case need to regenerate pfam searchable database: $ hmmpress Pfam-A.hmm top_hits.seqs.fa
	hmmscan --noali -o pfamout.temp --cut_tc --tblout $run.pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm $run.top_hits.seqs.fa 
	rm pfamout.temp
echo "use R script to process this list to only genes with both domains, and only the gene IDs. output is $run.both_domains_IDs.txt"
    Rscript $scripts_dir/sort_full_LRR_RLKs_olnly.R $run



echo "use the filtered list of gene IDs to collect sequence, make an alignment and tree"
fgrep -w --no-group-separator -A 1 -f $run.hits_with_both_domains_IDs.txt $run.top_hits.seqs.fa > $run.hits_with_both_domains_seqs.fa #pull only seqs from R screening
mkdir tree_from_hits
cp $run.hits_with_both_domains_seqs.fa tree_from_hits
cd tree_from_hits
cat $scripts_dir/scaffold_seqs.fa $run.hits_with_both_domains_seqs.fa >> $run.hits_and_scaffold_seqs.fa #add the backbone scaffold sequences to tree

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





