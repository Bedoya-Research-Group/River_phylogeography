River-weed phylogeography

1. Converted bam to fasta with samtools
`samtools fasta CRR103268.bam > CRR103268.fasta`


2. Index reference file
`bwa index CRR103268.fasta -b -t 4`


3. Map M. utile reads to reference
`bwa mem -t 8 CRR103268.fasta SRR12956182_Mutile_R1.trim.fq.gz SRR12956182_Mutile_R2.trim.fq.gz > mutile_aligned.sam`
