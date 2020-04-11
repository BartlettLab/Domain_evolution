#!/bin/bash

#BSUB -W 72:00             	# How much time does your job need (HH:MM)
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



#grep -A1 -f extract_node_1567_geneIDs.txt ../tree_seqs.fa --no-group-separator > 1567_tree_seqs.fa
#echo ">AT2G20850.1_SRF1_outgroup" >> 1567_tree_seqs.fa
#echo "MRSMRSGRDNNICFLGFLSFALISLPSLSLALTNPDDVAAINSLFLALESPLLPGWVASGGDPCGESWQGVLCNASQVETIILISANLGGELGVGLNMFTSLKAMDFSNNHIGGSIPSTLPVSLQNLFLSGNNFTGTIPESLSSLKSLSVMSLNNNLLSGKIPDVFQDLGLMINIDLSSNNLSGPLPPSMQNLSTLTSLLLQNNHLSGELDVLQDLPLKDLNVENNLFNGPIPEKLLSIPNFIKGGNLFNVTIAPSPSPETPPSPTSPKRPFFGPPSPNASAGHGQAHVRSPPSDHHPSRPTPQGKEDSFTSKRIIWISILGAFSFVVLALVCLLCGRKCLRKREDSEQLSKPHLTSEYGRAREGSRSNASMLPPSNTFNKDKEARPKERVGGASKLHGGAERSVGSESKQESHEIDMNGNAMDLMHPSSIPPIKRVIAKATEPAEASLKRTTSKSHGPLTAVKHFTVASLQQHTNSFSHENLIGTGMLGSVYRAELPGGKLFAVRKLDKKSPNHEEEGKFLELVNNIDRIRHANIVQLVGFCSEHSQRLLIHEYCRNGTLHDLLHIDDRLKIELSWNVRVRIALEAAKALEYLHEICDPPSIHRNFKSANILLDDDIRVHVSDCGLAPLISSGAVSQLSGQLLAAYGYGAPEFEYGIYTMKCDVYSFGVVMLELLTGRKSYDKKRDRGEQFLVRWAIPQLHDIDALAKMVDPSLKGDYPAKSLSHFADVISRCVQSEPEYRPLMSEVVQDLSDMIQREHRRNDSNGDNQYTGRR*" >> 1567_tree_seqs.fa




echo "Aligning..."
    mafft 1567_tree_seqs.fa > 1567_tree_alignment.fasta
    
echo "Filtering with noisy..."
    noisy -s 1567_tree_alignment.fasta


COMMENT

echo "inferring tree..."
    iqtree -s *out.fas -bb 1000 -nt 16 -m LG+F+R7
#    iqtree -s *out.fas -fast -nt 16 -m LG+F+R7
    


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



