#!/bash/bin/


for mapped_bam in *.mapped_reads.bam; do
    id=$(basename "$mapped_bam" .mapped_reads.bam)
    echo "Converting $mapped_bam to FASTQ"
    samtools fastq "$mapped_bam" \
        -1 "${id}_R1.fq" \
        -2 "${id}_R2.fq" \
        -0 /dev/null -s /dev/null -n
done
