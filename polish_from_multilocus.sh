#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ===================== USER SETTINGS =====================
MULTI_DIR=${MULTI_DIR:-/home/abedoya/river_phylogeo/assemblies03/per_locus_multifasta/*.multi.fasta}   # input: *.multi.fasta (one per locus)
READS_DIR=${READS_DIR:-/home/abedoya/river_phylogeo/assemblies03}                  # where FASTQs live
OUT_DIR=${OUT_DIR:-/home/abedoya/river_phylogeo/assemblies03/per_locus_polished}         # output multi-FASTAs (polished)
TMP_DIR=${TMP_DIR:-/home/abedoya/river_phylogeo/assemblies03/tmp/polish_work}                    # temp work

# Read file naming: SAMPLE must match the FASTA header prefix before the first '|'
# e.g., header ">AMB1311_Aguacate|NODE_..." -> sample="AMB1311_Aguacate"
# Set your R1/R2 patterns:
R1_SUFFIX=${R1_SUFFIX:--pe_mapped3_TE_R1.fq}
R2_SUFFIX=${R2_SUFFIX:--pe_mapped3_TE_R2.fq}

THREADS=${THREADS:-8}
# bcftools / samtools / bwa binaries (auto-find if in PATH)
BWA=${BWA:-bwa}
SAMTOOLS=${SAMTOOLS:-samtools}
BCFTOOLS=${BCFTOOLS:-bcftools}

# Variant calling thresholds (bcftools)
MIN_BASEQ=${MIN_BASEQ:-20}      # -Q (mpileup)
MIN_MAPQ=${MIN_MAPQ:-20}        # -q (mpileup)
PLOIDY=${PLOIDY:-2}
MIN_DP=${MIN_DP:-10}            # filter low depth
# =========================================================

mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "== Settings =="
echo "MULTI_DIR:  $MULTI_DIR"
echo "READS_DIR:  $READS_DIR"
echo "OUT_DIR:    $OUT_DIR"
echo "TMP_DIR:    $TMP_DIR"
echo "THREADS:    $THREADS"
echo "R1_SUFFIX:  $R1_SUFFIX"
echo "R2_SUFFIX:  $R2_SUFFIX"
echo

command -v "$BWA" >/dev/null || { echo "ERROR: bwa not found"; exit 1; }
command -v "$SAMTOOLS" >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
command -v "$BCFTOOLS" >/dev/null || { echo "ERROR: bcftools not found"; exit 1; }

# Helper: pull one sequence (header -> fasta) from a multi-fasta
extract_one_seq() {
  local multifasta="$1" hdr="$2" outfa="$3"
  awk -v h=">${hdr}" '
    $0==h {print; f=1; next}
    /^>/ && f {exit}
    f {print}
  ' "$multifasta" > "$outfa"
}

# Loop over each locus multi-fasta
for locus_fa in "$MULTI_DIR"/*.multi.fasta; do
  locus_base=$(basename "$locus_fa")
  locus_id="${locus_base%.multi.fasta}"
  echo ">> Locus: $locus_id"

  out_fa="$OUT_DIR/${locus_id}.polished.multi.fasta"
  : > "$out_fa"

  # List headers (one per sample) without the leading '>'
  mapfile -t headers < <(grep '^>' "$locus_fa" | sed 's/^>//')

  for hdr in "${headers[@]}"; do
    # hdr format: SAMPLE|CONTIG_ID ...
    sample="${hdr%%|*}"

    R1="$READS_DIR/${sample}${R1_SUFFIX}"
    R2="$READS_DIR/${sample}${R2_SUFFIX}"
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "  !! Missing reads for $sample (expected $R1 and $R2) — skipping"
      continue
    fi

    # Make tiny per-sample reference
    work="$TMP_DIR/${locus_id}__${sample}"
    mkdir -p "$work"
    ref="$work/ref.fa"
    extract_one_seq "$locus_fa" "$hdr" "$ref"

    # Index ref, map, sort
    "$BWA" index "$ref" >/dev/null 2>&1
    "$BWA" mem -t "$THREADS" "$ref" "$R1" "$R2" \
      | "$SAMTOOLS" view -b - \
      | "$SAMTOOLS" sort -@ "$THREADS" -o "$work/align.sorted.bam"
    "$SAMTOOLS" index "$work/align.sorted.bam"

    # Call variants with bcftools (nucleotide; diploid)
    # mpileup -> call (multiallelic), bgzip+index
    "$BCFTOOLS" mpileup -f "$ref" -q "$MIN_MAPQ" -Q "$MIN_BASEQ" -Ou "$work/align.sorted.bam" \
      | "$BCFTOOLS" call -mv -Ou --ploidy "$PLOIDY" \
      | "$BCFTOOLS" filter -e "FMT/DP<${MIN_DP}" -Ou \
      | "$BCFTOOLS" norm -f "$ref" -Ob -o "$work/var.bcf"
    "$BCFTOOLS" index "$work/var.bcf"

    # Build consensus with IUPAC codes at hets
    cons="$work/consensus.fa"
    "$BCFTOOLS" consensus -f "$ref" -I "$work/var.bcf" > "$cons"

    # Write to per-locus polished multiFASTA with header >SAMPLE|CONTIG_ID
    # Replace original header line with ">SAMPLE|CONTIG_ID" exactly (no description tail)
    seq=$(awk 'NR>1{gsub(/\r/,""); printf "%s",$0} END{printf "\n"}' "$cons")
    printf ">%s\n%s" "$hdr" "$seq" >> "$out_fa"

    # (optional) clean tiny workspace
    rm -rf "$work"
  done

  echo "  -> Wrote $out_fa"
done

echo "All done."

