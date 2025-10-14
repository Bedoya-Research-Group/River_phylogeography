#!/usr/bin/env bash
# align_all_loci_with_mafft.sh
# Align per-locus multi-FASTAs (from 01_data) using MAFFT.

set -euo pipefail

BASE="${1:-/path/to/HybPhaser/exons/01_data}"   # where <locus>.fasta live
OUTDIR="${2:-$(dirname "$BASE")/02_alignments}"
JOBS="${JOBS:-8}"      # number of loci to align in parallel (processes)
THREADS_PER_JOB="${THREADS_PER_JOB:-1}"   # MAFFT threads per locus (keep modest if using parallel jobs)

mkdir -p "$OUTDIR/logs"

# Find all locus files (top-level in 01_data)
mapfile -t LOCUS_FASTAS < <(find "$BASE" -maxdepth 1 -type f -name "*.fasta" | sort)

if (( ${#LOCUS_FASTAS[@]} == 0 )); then
  echo "No locus FASTAs found in: $BASE"
  exit 1
fi

# Function to align one file
align_one() {
  in="$1"
  out="$2"
  log="$3"

  # skip files with <2 sequences
  nseq=$(grep -c '^>' "$in" || true)
  if (( nseq < 2 )); then
    echo "[SKIP] $(basename "$in"): only $nseq sequence(s)" | tee -a "$log"
    return 0
  fi

  # MAFFT: --auto chooses a good strategy; --adjustdirectionaccurately safeguards strand
  # If you KNOW all sequences are same strand already, you can drop --adjustdirectionaccurately to go faster.
  mafft --auto --thread "$THREADS_PER_JOB" --adjustdirectionaccurately \
    "$in" > "$out" 2> "$log"
}

export -f align_one
export OUTDIR THREADS_PER_JOB

# Run in parallel (POSIX xargs). Adjust -P to taste or set JOBS env var.
printf "%s\0" "${LOCUS_FASTAS[@]}" \
| xargs -0 -I{} -P "$JOBS" bash -c '
  in="$1"
  locus=$(basename "${in%.fasta}")
  out="'"$OUTDIR"'/${locus}.aln.fasta"
  log="'"$OUTDIR"'/logs/${locus}.log"
  align_one "$in" "$out" "$log"
' _ {}

echo "Done. Alignments are in: $OUTDIR"
echo "Logs are in:             $OUTDIR/logs"
echo "Tip: set JOBS=16 THREADS_PER_JOB=1 (or JOBS=8 THREADS_PER_JOB=2) for your CPUs."
