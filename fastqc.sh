for infile in /home/abedoya/river_phylogeo/*_1.fq.gz
        do
        base=$(basename ${infile} _1.fq.gz)
        fastqc -t 10 -o fastqc_reports ${infile} ${base}_2.fq.gz
done
~    
