#!/bin/bash

#BSUB -W 0:01             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=100]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 1                	# Where X is in the set {1..X}
#BSUB -J refine         	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}

runs="CLV1 GSO HAE HBD PEPR"
echo "Director script, this will execute follow-up analyses for the following clades: $runs"
for run in $runs
do
    echo "Now submitting job for $run clade"
    cd /home/jm33a/domain_evolution/clade_trees/$run/tree_from_hits/
    bsub -q long -J "$run_followup" -W 72:00 -R rusage[mem=500] -n 17 -R span[hosts=1] -o out.%J -e err.%J "sh /home/jm33a/domain_evolution/scripts_and_resources/1KP_followup_analyses.sh $run"
    echo "Finished submitting job for $run clade"
done
echo "Finished submitting jobs for all clades"
exit
