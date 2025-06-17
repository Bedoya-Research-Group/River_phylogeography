for infile in /home/abedoya/river_phylogeo/*_1.fq.gz
	do
	base=$(basename ${infile} _1.fq.gz)
	fastp -i ${infile} -I ${base}_2.fq.gz -o ${base}_R1_trimmed.fq.gz -O ${base}_R2_trimmed.fq.gz --detect_adapter_for_pe -q 20 -l 50 -h ${base}_fastp.html -j ${base}_fastp.json	
	done
