#!/usr/bin/env bash
# Create per-locus multi-FASTA files from HybPhaser 01_data/*/consensus/*.fasta

set -euo pipefail

BASE="${1:-.}"
OUTDIR="${2:-$BASE}"

shopt -s nullglob

# 1) Gather all loci present across samples
mapfile -t LOCI < <(find "$BASE" -mindepth 2 -maxdepth 2 -type d -name consensus \
  -exec bash -c 'for f in "$1"/*.fasta; do [[ -e "$f" ]] && basename "${f%.fasta}"; done' _ {} \; \
  | sort -u)

if (( ${#LOCI[@]} == 0 )); then
  echo "No consensus FASTAs found under: $BASE" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

# 2) Build one multi-FASTA per locus
for locus in "${LOCI[@]}"; do
  out="$OUTDIR/${locus}.fasta"
  : > "$out"   # truncate/create

  # loop over samples
  for sdir in "$BASE"/*/ ; do
    sample="$(basename "$sdir")"
    f="$sdir/consensus/${locus}.fasta"
    [[ -s "$f" ]] || continue

    # Reheader first line to >SAMPLE-LOCUS, keep the rest as-is
    # (HybPhaser writes one sequence per file; if multiple appear, this keeps them all under the same header)
    awk -v hdr=">${sample}-${locus}" '
      NR==1 && $0 ~ /^>/ {print hdr; next}
      {print}
    ' "$f" >> "$out"
  done

  # Remove empty files (no samples had this locus)
  if ! grep -q "^>" "$out"; then rm -f "$out"; fi
done

echo "Done. Wrote per-locus multi-FASTAs to: $OUTDIR"
