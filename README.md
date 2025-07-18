Riverweed phylogeography


## 1. PREPARING REFERENCE GENOME

1.1. Convert BAM to fasta with samtools
`samtools fasta CRR103268.bam > CRR103268.fasta`

1.2. Index reference file
`bwa index CRR103268.fasta`

1.3. Map M. utile reads to reference
`bwa mem CRR103268.fasta SRR12956182_Mutile_R1.trim.fq.gz SRR12956182_Mutile_R2.trim.fq.gz -t 22 > mutile_aligned.sam`

1.4. Convert SAM to BAM
`samtools view -@ 20 -b SRR12956182_Mutile-pe.sam > SRR12956182_Mutile-pe.bam`

1.5. Sort BAM
`samtools sort -@ 20 -o SRR12956182_Mutile-pe.sorted01.bam SRR12956182_Mutile-pe.bam`

1.6. bcftools pileup
` bcftools mpileup -Ou -f CRR103268.fasta SRR12956182_Mutile-pe.sorted01.bam | bcftools call -mv -Oz -o SRR12956182_Mutile-pe.sorted01.bcf`

1.7. Filter variants
`bcftools view -Oz -i 'QUAL>20 || DP>5' SRR12956182_Mutile-pe.sorted01.bcf > SRR12956182_Mutile-pe.filtered01.vcf.gz`

1.8. index
`tabix SRR12956182_Mutile-pe.filtered01.vcf.gz`

1.9. Normalize
`bcftools norm -f CRR103268.fasta -m +any SRR12956182_Mutile-pe.filtered01.vcf.gz -Oz -o SRR12956182_Mutile.norm.filtered01.vcf.gz`

1.10. Index
`tabix SRR12956182_Mutile.norm.filtered01.vcf.gz`

1.11. Mutate consensus
`cat CRR103268.fasta | bcftools consensus -s SRR12956182_Mutile-pe.sorted01.bam -H R SRR12956182_Mutile.norm.filtered01.vcf.gz > CRR103268_01.fasta`

1.12. Repeat 2 times



## 2. Genome-skimming data prep

2.1. Check Novogene sequence download
`for dir in */; do   if [[ -f "$dir/MD5.txt" ]];` then     `echo "Checking MD5 in $dir";     (cd "$dir" && md5sum -c MD5.txt);   else     echo "No MD5.txt found in $dir";   fi; done`

2.2. Check data quality `bash fastqc.sh`

2.3. Trim reads `bash fastp.sh`

2.4. Change names to correspond to voucher numbers and taxa in the study `bash name_change.sh` You will need the name_map.tsv file



## 3. WHOLE GENOME DATA ANALYSIS

3.1. Map reads to mutated reference genome `bash bwa-mem.sh`

3.2. Convert sam to bam `bash samtobam.sh` and sort bam `bam_sort.sh`

3.3. Call SNPs `bcftools mpileup -Ou -f CRR103268_03_mutated.fasta *-pe.sorted.bam | bcftools call -mv -Oz -o all_samples_raw.vcf.gz`

3.4 Filter SNPs `vcftools --gzvcf all_samples_raw.vcf.gz --remove-indels --minQ 30 --minDP 20 --maf 0.05 --max-missing 0.80 --thin 10000 --recode --out all_samples_filtered.vcf.gz`

3.5 Remove individuals with >50% missing data (88 indiv left) `vcftools --vcf all_samples_filtered.vcf.gz.recode.vcf --missing-indv --out missing_data_individuals`

`vcftools --vcf all_samples_filtered.vcf.gz.recode.vcf --remove individuals_to_remove_from_all.txt --recode --recode-INFO-all --out all_samples_filtered_final`


Data subsets:

For Marathrum all (87 indiv left):
`vcftools --vcf all_samples_filtered.vcf.gz.recode.vcf --remove individuals_to_remove_from_marathrum.txt --max-missing 0.80 --recode --recode-INFO-all --out marathrum_samples_filtered_final`

Filtering missing data for admixture analysis:
`vcftools --vcf all_samples_filtered.vcf.gz.recode.vcf --remove individuals_to_remove_from_marathrum.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_samples_filtered_final_095`

For Marathrum Panama (69 indiv left):
`vcftools --vcf marathrum_samples_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_panama.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_panama_filtered_final`

For Colon Province (57 indiv left):
`vcftools --vcf marathrum_panama_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_colon.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_colon_filtered_final`

For Rio Diego (26 indiv left):
`vcftools --vcf marathrum_colon_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_diego.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_diego_filtered_final`

For Rio Aguacate (20 indiv left):
`vcftools --vcf marathrum_colon_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_aguacate.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_aguacate_filtered_final`

For Rio Cocle del Norte (11 indiv left):
`vcftools --vcf marathrum_colon_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_cocle_norte.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_cocle_norte_filtered_final`

For Chiriqui Province (12 indiv left):
`vcftools --vcf marathrum_panama_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_chiriqui.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_chiriqui_filtered_final`

For Paraiso Province (10 indiv left):
`vcftools --vcf marathrum_chiriqui_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_paraiso.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_paraiso_filtered_final`

For Aguacate+Diego (46 individuals):
`vcftools --vcf marathrum_colon_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_aguacate_diego.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_aguacate_diego_filtered_final`

For Aguacate+Cocle (31 individuals):
`vcftools --vcf marathrum_colon_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_aguacate_cocle.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum_aguacate_cocle_filtered_final`

For Diversity stats (gene flow within rivers) (74 indiv left):
`vcftools --vcf marathrum_panama_filtered_final.recode.vcf --remove individuals_to_remove_from_marathrum_diversity.txt --max-missing 0.95 --recode --recode-INFO-all --out marathrum4diversity_stats`

3.6 Convert vcf to plink:
`vcftools --vcf each_vcf_above.vcf --plink --out each_vcf_above_plink`

3.7 Convert to binary plink:
`plink --file each_resulting_file_above_plink --out each_resulting_file_above_bin_plink --make-bed --noweb`

3.8 Run admixture:

`bash admixture.sh`
*****

## 4. ASSEMBLY OF TE DATA FROM GENOME SKIMMING

4.1. Map reads to Bedoya et al., 2021 target file ('Target_Sequences_single_exon.fasta'; previously indexed) `bash bwa-mem.sh`

4.2. Convert sam to bam `bash samtobam.sh`
****************
I am trying this alternative
4.3 Extract mapped reads for each individual `bash extract_mapped_reads_per_individual.sh`

4.4 Convert mapped BAMs to Fastq `bash mapped_bams2fastq.sh`

4.5 De novo assemble `bash assembly.sh`

4.6 For each individual, map reads to the resulting scaffolds (previously indexed) `bash mapping_round2.sh`

4.7. Repeat step 4.2 to 4.6 X2 times
***************

4.3. Sort bam `bash bam_sort.sh`

4.4. Extract contigs names (target sequence loci names) `bash mapped_contigs_names.sh`

4.5. Index bam `bash samtools_index.sh`

4.6. Extract contigs from Target reference `bash extract_contigs.sh`

4.7. Assemble scaffolds with SPAdes `bash assembly.sh`

4.8. Extract the best contig `bash best_contig_extraction.sh`

4.9. Replace reference with assembled contigs


Alternatively, run Hybpiper and the PPD. PPD:
1. `bash putative_paralog/Step1.sh . namelist.txt`



## 5. PLASTOME ASSEMBLY

5.1. `bash getorganelle.sh`

5.2. Check assemblies in Bandage

5.3. Rename *.fq files (`bash files_folders_name_changes.sh`) so they can be move to a separate folder for Plastome SNP calls

5.4. Annotate plastomes `perl PGA/PGA.pl -r /home/abedoya/river_phylogeo/plastome/ -t assembled_plastomes/`



## 6. Plastid genomes SNP analysis

6.1. Map reads and call SNPs like done with the mutate reference genome (above).

6.2. Use bedtools to mask IRb in reference as defined by Bedoya et al., 2019
`bedtools maskfasta -fi Mutile_MN165814_plastome.fasta -bed IR_mask.bed -fo Mutile_MN165814_plastome_masked.fasta`

6.3. Index Masked fasta `samtools faidx Mutile_MN165814_plastome_masked.fasta`

6.4. Call SNPs `bcftools mpileup -Ou -f ../Mutile_MN165814_plastome_masked.fasta *sorted_plastome.bam | bcftools call -mv -Oz -o all_samples_raw.vcf.gz`

6.5. Filter SNPs `vcftools --gzvcf all_samples_raw.vcf.gz --remove-indels --minQ 30 --minDP 20 --maf 0.05 --max-missing 0.95 --recode --out all_samples_filtered.vcf`

6.6. Plot results

3. Sort bam `bash bam_sort.sh`


## 7. Gene flow and genetic diversity across space

See scripts in spatial_gene_flow_gendiv.R


## 8. Demographic modeling+phylogenetic inference

8.1. Phylogenetic inference
`vcftools --vcf ../all_samples_filtered_final.recode.vcf --remove individuals_to_remove_for_tree_inference.txt --max-missing 0.80 --recode --recode-INFO-all --out marathrum_for_tree_inference`
`python vcf2phylip.py -i marathrum_for_tree_inference.recode.vcf`
`raxml-ng --all --msa marathrum_for_tree_inference.recode.min4.phy --model GTR+G --bs-metric fbp --bs-trees autoMRE --threads 10`

8.2. Down projection preview `easySFS.py -i marathrum_colon_filtered_final.recode.vcf -p pops_file.txt  --preview`

8.3. Down projection `easySFS.py -i marathrum_colon_filtered_final.recode.vcf -p pops_file.txt  --proj 48,36,20 -a -v`
