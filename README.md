Riverweed phylogeography

## 1. PREPARING REFERENCE GENOME


1. Convert BAM to fasta with samtools
`samtools fasta CRR103268.bam > CRR103268.fasta`

2. Index reference file
`bwa index CRR103268.fasta`

3. Map M. utile reads to reference
`bwa mem CRR103268.fasta SRR12956182_Mutile_R1.trim.fq.gz SRR12956182_Mutile_R2.trim.fq.gz -t 22 > mutile_aligned.sam`

4. Convert SAM to BAM
`samtools view -@ 20 -b SRR12956182_Mutile-pe.sam > SRR12956182_Mutile-pe.bam`

5. Sort BAM
`samtools sort -@ 20 -o SRR12956182_Mutile-pe.sorted01.bam SRR12956182_Mutile-pe.bam`

6. bcftools pileup
`bcftools view -Oz -i 'QUAL>20 || DP>5' SRR12956182_Mutile-pe.sorted01.bcf > SRR12956182_Mutile-pe.filtered01.vcf.gz`

7. index
`tabix SRR12956182_Mutile-pe.filtered01.vcf.gz`

8. Normalize
`bcftools norm -f CRR103268.fasta -m +any SRR12956182_Mutile-pe.filtered01.vcf.gz -Oz -o SRR12956182_Mutile.norm.filtered01.vcf.gz`

9. Index
`tabix SRR12956182_Mutile.norm.filtered01.vcf.gz`

10. Mutate consensus
`cat CRR103268.fasta | bcftools consensus -s SRR12956182_Mutile-pe.sorted01.bam -H R SRR12956182_Mutile.norm.filtered01.vcf.gz > CRR103268_01.fasta`

11. Repeat 2 times


## 2. Genome-skimming data analysis


1.Check Novogene sequence download
`for dir in */; do   if [[ -f "$dir/MD5.txt" ]];` then     `echo "Checking MD5 in $dir";     (cd "$dir" && md5sum -c MD5.txt);   else     echo "No MD5.txt found in $dir";   fi; done`

2. Check data quality `bash fastqc.sh`

3. Trim reads `bash fastp.sh`

4. Change names to correspond to voucher numbers and taxa in the study `bash name_change.sh` You will need the name_map.tsv file



## 3. ASSEMBLY OF TE DATA FROM GENOME SKIMMING


1. Map reads to Bedoya et al., 2021 target file ('Target_sequences.fasta') `bash bwa-mem.sh`

2. Convert sam to bam `bash samtobam.sh`

3. Sort bam `bash bam_sort.sh`

4. Extract contigs names (target sequence loci names) `bash mapped_contigs_names.sh`

5. Index bam `bash samtools_index.sh`



## 4. PLASTOME ASSEMBLY

1. `bash getorganelle.sh`



## 5. WHOLE GENOME DATA ANALYSIS

1. Map reads to mutated reference genome `bash bwa-mem.sh`

2. Convert sam to bam `bash samtobam.sh`

3. Sort bam `bash bam_sort.sh`