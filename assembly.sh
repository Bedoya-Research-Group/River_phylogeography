#!/bin/bash


for target_dir in extracted_TE_contigs/*
	do
	target=$(basename "$target_dir")

	for r1 in "$target_dir"/*_R1.fq
	do
	r2="${r1/_R1.fq/_R2.fq}"
    	sample=$(basename "$r1" _R1.fq)

	outdir="${target_dir}/${sample}_assembly"

    	echo "Running SPAdes for $sample in $target_dir"

    	#Run SPAdes
    	spades.py -1 "$r1" -2 "$r2" -o "$outdir" --careful --sc  -t 8
  	done
done
