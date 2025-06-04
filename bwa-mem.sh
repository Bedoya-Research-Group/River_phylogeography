#!/bin/bash

#for infile in /home/abedoya/river_phylogeo/*_R1_trimmed.fq.gz
#        do
#        base=$(basename ${infile} _R1_trimmed.fq.gz)
#        bwa mem CRR103268_03_mutated.fasta ${infile} ${base}_R2_trimmed.fq.gz -t 20 > ${base}-pe_mapped.sam
#        done


for infile in /home/abedoya/river_phylogeo/*_R1_trimmed.fq.gz
        do
        base=$(basename ${infile} _R1_trimmed.fq.gz)
        bwa mem Target_Sequences.fasta ${infile} ${base}_R2_trimmed.fq.gz -t 20 > ${base}-pe_mapped_TE.sam
        done