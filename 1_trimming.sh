for infile in /Volumes/Bedoya/Genome/*R1.fastq.gz
	do
	base=$(basename ${infile} R1.fastq.gz)
	trimmomatic PE -threads 15 ${infile} ${base}R2.fastq.gz ${base}R1.trim.fq.gz ${base}R1_un.trim.fq.gz ${base}R2.trim.fq.gz ${base}R2_un.trim.fastq.gz ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done

