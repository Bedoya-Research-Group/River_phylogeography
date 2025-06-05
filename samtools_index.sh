#!/bin/bash

for infile in /home/abedoya/river_phylogeo/*-pe.sorted_TE.bam
        do
        samtools index ${infile}
        done
