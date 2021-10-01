#!/bin/bash -l 
 
# Set SCC Project 
#$ -P bunuc 
# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time. 
#$ -l h_rt=120:00:00 
 
# Send an email when the job finishes or if it is aborted (by default no email is sent). 
#$ -m ea 
 
# Reuest a node with at least 4 GB of memory per core 
#$ -l mem_per_core=8G 
 
# Request a paralell environemtn with _ cores 
#$ -pe omp 1 
 
# Give job a name 
#$ -N sand_5
 
# Combine output and error files into a single file 
#$ -j y 
 
# Specify the output file name 
#$ -o sand_output.txt 
 
# Keep track of information related to the current job 
 
declare -i timer=5000000
 
declare -i L=32
 
disp=0
 
drive=1
 
declare -i rand=5
 
echo -e "$timer\n$L\n$disp\n$drive\n$rand" |  ./sandpile_manna 
