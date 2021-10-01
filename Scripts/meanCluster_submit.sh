#!/bin/bash -l 
 
# Set SCC Project 
#$ -P bunuc 
# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time. 
#$ -l h_rt=12:00:00 
 
# Reuest a node with at least 4 GB of memory per core 
#$ -l mem_per_core=8G 
 
# Request a paralell environemtn with _ cores 
#$ -pe omp 2
 
# Give job a name 
#$ -N mm_meanCluster
 
# Combine output and error files into a single file 
#$ -j y 
 
# Specify the output file name 
#$ -o results/meanCluster.txt
 
# Keep track of information related to the current job 
 
module load python3
 
python3 meanCluster.py
 
