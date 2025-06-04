#!/bin/bash


for infile in /home/abedoya/river_phylogeo/*-pe.sorted_TE.bam
        do
        base=$(basename ${infile} -pe.sorted_TE.bam)
        samtools idxstats ${infile} | awk '$3 > 0 { print $1 }' > ${base}_mapped_contigs.txt
        done