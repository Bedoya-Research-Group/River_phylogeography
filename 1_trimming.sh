for infile in /media/ambedoya/extradrive1/Lagomarsino/Palicourea/RAW/*R1_001.fastq.gz
	do
	base=$(basename ${infile} R1_001.fastq.gz)
	trimmomatic PE -threads 15 ${infile} ${base}R2_001.fastq.gz ${base}R1_001.trim.fq.gz ${base}R1_001un.trim.fq.gz ${base}R2_001.trim.fq.gz ${base}R2_001un.trim.fastq.gz ILLUMINACLIP:all_adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done

