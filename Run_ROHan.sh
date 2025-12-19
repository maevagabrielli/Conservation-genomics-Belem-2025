#!/bin/bash

### Here we work with one chromosome, and one individual for each group
# chromosome SUPER_17, 40M
# waglerianus: DB23461 LC: LC03 ST: ST02

# define parameters of the analysis
size=100000 # ROH length of 1Mb
cores=4 # 4 threads for each run
mkdir -p window$size #create the output directory if not already present
for ind in DB23461 LC03 ST02 # loop over the 3 individuals of interest
do # open the loop
bam="$ind"_SUPER_17.bam # individual bam file
sbatch --time=2:00:00 --nodes=1 --ntasks=1 --mem=4G --cpus-per-task="$cores" -J $ind -e ROHan_"$ind".e -o ROHan_"$ind".o 2a_ROHan.sh $ind $bam $size $cores #launch ROHan with the parameters defined in the loop. As there are 3 individuals, this will launch 3 jobs
done # close the loop
