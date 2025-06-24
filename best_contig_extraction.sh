#!/bin/bash

for dir in */*_assembly/; do
  # Remove trailing slash and get the folder name
  base=$(basename "$dir")
  
  # Remove the _assembly suffix to get the target
  target="${base%_assembly}"
  
  # Copy the file
  cp "${dir}contigs.fasta" "${target}_contig.fasta"
done
