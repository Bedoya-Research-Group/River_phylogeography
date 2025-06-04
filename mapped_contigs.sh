#!/bin/bash

#for infile in /home/abedoya/river_phylogeo/*_R1_trimmed.fq.gz
#        do
#        base=$(basename ${infile} _R1_trimmed.fq.gz)
#        bwa mem CRR103268_03_mutated.fasta ${infile} ${base}_R2_trimmed.fq.gz -t 20 > ${base}-pe_mapped.sam
#        done


for infile in /home/abedoya/river_phylogeo/*-pe.sorted_TE.bam
        do
        base=$(basename ${infile} -pe.sorted_TE.bam)
        samtools idxstats ${infile} | awk '$3 > 0 { print $1 }' > ${base}_mapped_contigs.txt
        done