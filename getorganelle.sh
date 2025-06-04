#!/bin/bash

#for infile in /home/abedoya/river_phylogeo/*_R1_trimmed.fq.gz
#        do
#        base=$(basename ${infile} _R1_trimmed.fq.gz)
#        bwa mem CRR103268_03_mutated.fasta ${infile} ${base}_R2_trimmed.fq.gz -t 20 > ${base}-pe_mapped.sam
#        done


for infile in /home/abedoya/river_phylogeo/*_R1_trimmed.fq.gz
        do
        base=$(basename ${infile} _R1_trimmed.fq.gz)
        get_organelle_from_reads.py -1 ${infile} -2 ${base}_R2_trimmed.fq.gz -o ${base}_plastome_output -R 15 -k 21,45,65,85,105 -F embplant_pt
        done