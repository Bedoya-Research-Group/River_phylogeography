#!/usr/bin/env bash
set -euo pipefail

for d in *-pe_mapped3_TE_assembly; do
    # skip if not a directory
    [[ -d "$d" ]] || continue

    name=$(basename "$d" -pe_mapped3_TE_assembly)   # strip suffix
    infile="$d/contigs.fasta"
    outfile="${name}_contigs_filtered.fasta"

    if [[ -f "$infile" ]]; then
        echo "Filtering $infile → $outfile"
        seqkit seq -m 200 -M 10000 "$infile" > "$outfile"
    else
        echo "Warning: $infile not found" >&2
    fi
done
