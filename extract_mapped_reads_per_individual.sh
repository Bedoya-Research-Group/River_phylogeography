#!/bin/bash

for bam in *mapped.bam; do
    id=$(basename "$bam" .bam)
    echo "Extracting mapped reads for $id"
    samtools view -b -F 4 "$bam" > "${id}.mapped_reads.bam"
done
