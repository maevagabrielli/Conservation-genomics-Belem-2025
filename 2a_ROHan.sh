#!/bin/bash

# recover parameters defined in Run_ROHan.sh
ind=$1 # first parameter
bam=$2 # second parameter
size=$3 # third parameter
cores=$4 # fourth parameter

# load ROHan
module load bioinfo/ROHan/1.0.1
# Run ROHan
rohan -t $cores -o window$size/$ind --size $size --rohmu 1e-05 --map ../data/rPodRaf1.pri.cur.mt.softmasked.norepeats.SUPER_17.bed --auto ../data/chr.list ../data/SUPER_17.fasta ../data/$bam
