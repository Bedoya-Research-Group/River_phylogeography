#!/bin/bash

for infile in /home/abedoya/river_phylogeo/*-pe.sorted_TE.bam
        do
        sample=$(basename ${infile} -pe.sorted_TE.bam)

        samtools idxstats "$infile" | awk '$3 > 0 {print $1}' > tmp_contigs.txt

        while read contig
                do
                echo "Extracting $contig from $sample"

                mkdir -p extracted_TE_contigs/${contig}

                samtools view -h "${infile}" "$contig" | samtools fastq -1 extracted_TE_contigs/${contig}/${sample}_R1.fq -2 extracted_TE_contigs/${contig}/${sample}_R2.fq -0 /dev/null -s /dev/null -n -
                done < tmp_contigs.txt
        done
