#!/usr/bin/env bash
set -euo pipefail
ulimit -n 8192 || true
trap 'echo "ERROR on line $LINENO"; exit 1' ERR
shopt -s nullglob

# ================= Tunables =================
TARGET=${TARGET:-Target_Sequences_single_exon.fasta}      # multi-FASTA of loci
IN_GLOB=${IN_GLOB:-*filtered.fasta}           # per-sample contig FASTAs
WORKDIR=${WORKDIR:-$(pwd)/per_locus_work}
OUTDIR=${OUTDIR:-$(pwd)/per_locus_multifasta}
THREADS=${THREADS:-32}

# MMseqs nucleotide search knobs
MIN_ID=${MIN_ID:-0.50}
QCOV=${QCOV:-0.20}
COV_MODE=${COV_MODE:-0}                       # 0=target, 1=query, 2=both
SENS=${SENS:-7.5}                             # sensitivity (higher=more sensitive/slower)
SEARCH_TYPE=${SEARCH_TYPE:-3}                 # 3 = nucleotide vs nucleotide

# ===============================================================
echo "== Settings =="
echo "TARGET:       $TARGET"
echo "IN_GLOB:      $IN_GLOB"
echo "WORKDIR:      $WORKDIR"
echo "OUTDIR:       $OUTDIR"
echo "THREADS:      $THREADS"
echo "MIN_ID:       $MIN_ID"
echo "QCOV:         $QCOV (cov-mode=$COV_MODE)"
echo "SENS:         $SENS"
echo "SEARCH_TYPE:  $SEARCH_TYPE"
echo

command -v mmseqs >/dev/null 2>&1 || { echo "ERROR: mmseqs not in PATH"; exit 127; }
[[ -f "$TARGET" ]] || { echo "ERROR: Target FASTA not found: $TARGET"; exit 2; }

INPUTS=( $IN_GLOB )
(( ${#INPUTS[@]} > 0 )) || { echo "ERROR: No inputs matched: $IN_GLOB"; exit 3; }

mkdir -p "$WORKDIR"/{mmseqs_db,tmp,hits} "$OUTDIR"

# -------- helper: infer sample name from filename --------
infer_sample() {
  local base="${1%.*}"
  [[ "$base" == cluster_* ]] && base="${base#cluster_}"
  for sent in "_contigs" "_scaffolds" "_NODE" "_filtered" "_assembly" "_reads" "_contig"; do
    [[ "$base" == *"$sent"* ]] && { base="${base%%${sent}*}"; break; }
  done
  IFS='_' read -r t1 t2 _ <<< "$base"
  if [[ -n "$t1" && -n "$t2" ]]; then
    local n; n=$(awk -F'_' '{print NF}' <<< "$base")
    (( n >= 3 )) && base="${t1}_${t2}"
  fi
  echo "$base"
}

# -------- 1) Build combined SAMPLE-tagged FASTA --------
COMBINED="$WORKDIR/mmseqs_db/all_samples_tagged.fa"
: > "$COMBINED"

echo ">> Building combined SAMPLE-tagged FASTA from ${#INPUTS[@]} files:"
SAMPLE_LIST_FILE="$WORKDIR/hits/all_samples.txt"
: > "$SAMPLE_LIST_FILE"

for f in "${INPUTS[@]}"; do
  s=$(infer_sample "$(basename "$f")")
  echo "   + $(basename "$f")  -> sample=$s"
  echo "$s" >> "$SAMPLE_LIST_FILE"
  sed -e "s/^>/>${s}|/" "$f" >> "$COMBINED"
done
sort -u "$SAMPLE_LIST_FILE" -o "$SAMPLE_LIST_FILE"

HDRS=$(grep -c '^>' "$COMBINED" || echo 0)
[[ "$HDRS" -gt 0 ]] || { echo "ERROR: Combined FASTA has 0 headers"; exit 4; }

# -------- 2) MMseqs search (nucleotide) --------
QDB="$WORKDIR/mmseqs_db/qDB"
TDB="$WORKDIR/mmseqs_db/tDB"
ALN="$WORKDIR/mmseqs_db/aln"
RAW_HITS="$WORKDIR/hits/raw_hits.tsv"

echo
echo ">> Creating nucleotide DBs"
mmseqs createdb "$COMBINED" "$QDB" --dbtype 2
mmseqs createdb "$TARGET"   "$TDB" --dbtype 2

echo ">> Running MMseqs search"
mmseqs search "$QDB" "$TDB" "$ALN" "$WORKDIR/tmp" \
  --search-type 3 \
  --min-seq-id "$MIN_ID" -c "$QCOV" --cov-mode "$COV_MODE" \
  --threads "$THREADS" -s "$SENS"

echo ">> Converting alignments to TSV"
mmseqs convertalis "$QDB" "$TDB" "$ALN" "$RAW_HITS" \
  --format-output "query,target,pident,alnlen,qlen,tlen,qcov,tcov,evalue,bits"

echo ">> Example hits:"
head -n 5 "$RAW_HITS" || true

# -------- 3) Choose ONE best contig per (target, sample); exclude multi-target contigs --------
MAP_TSV="$WORKDIR/hits/target_sample_to_query.tsv"     # included: target \t sample \t query
WARNINGS="$WORKDIR/hits/multi_target_warnings.log"     # contigs that hit >1 target (excluded)
SUMMARY_TSV="$WORKDIR/hits/summary_per_target.tsv"     # per-target counts + distinct samples + missing

export RAW_HITS MAP_TSV WARNINGS SUMMARY_TSV COMBINED SAMPLE_LIST_FILE
python3 - <<'PY'
import os, csv, sys
from collections import defaultdict

raw = os.environ["RAW_HITS"]
map_tsv = os.environ["MAP_TSV"]
warnfile = os.environ["WARNINGS"]
summary = os.environ["SUMMARY_TSV"]
sample_list_file = os.environ["SAMPLE_LIST_FILE"]

# Load all samples (expected)
with open(sample_list_file) as f:
    all_samples = sorted({line.strip() for line in f if line.strip()})

# Read raw hits and aggregate
# Keep fields we need for ranking
Hit = tuple  # (q,t,pident,qcov,alnlen,bits)
q2t = defaultdict(set)
hits_by_t_s = defaultdict(list)  # (t,sample) -> [Hit...]

with open(raw) as f:
    r = csv.reader(f, delimiter='\t')
    for row in r:
        if not row or len(row) < 10: 
            continue
        q = row[0].split()[0]
        t = row[1].split()[0]
        try:
            pident = float(row[2])
            alnlen = int(row[3])
            qlen   = int(row[4])
            tlen   = int(row[5])
            qcov   = float(row[6])
            tcov   = float(row[7])
            evalue = row[8]  # string is fine
            bits   = float(row[9])
        except Exception:
            continue
        q2t[q].add(t)
        sample = q.split("|",1)[0]
        hits_by_t_s[(t, sample)].append((q,t,pident,qcov,alnlen,bits))

# Identify multi-target contigs (to exclude)
multi_target = {q for q,ts in q2t.items() if len(ts) > 1}

# For each (target, sample), pick the best contig based on:
# 1) highest bitscore, then 2) highest pident, then 3) highest qcov, then 4) longest alnlen
best_choice = {}
for (t,s), hs in hits_by_t_s.items():
    # filter out multi-target queries first
    hs = [h for h in hs if h[0] not in multi_target]
    if not hs:
        continue
    hs.sort(key=lambda x: (x[5], x[2], x[3], x[4]), reverse=True)  # bits, pident, qcov, alnlen
    best_choice[(t,s)] = hs[0]  # keep best

# Write warnings for excluded contigs
excluded = sorted(multi_target)
if excluded:
    with open(warnfile, 'w') as w:
        for q in excluded:
            ts = sorted(q2t[q])
            w.write(f"WARNING: query {q} matched multiple targets {','.join(ts)}; excluded\n")

# Write mapping: target \t sample \t query (one per sample)
with open(map_tsv, 'w', newline='') as g:
    w = csv.writer(g, delimiter='\t')
    for (t,s), (q, *_rest) in sorted(best_choice.items()):
        w.writerow([t, s, q])

# Build per-target summary
t_included = defaultdict(int)      # number of included samples (i.e., sequences) for target
t_missing  = defaultdict(list)     # which samples are missing per target
t_distinct_samples = defaultdict(set)

# All targets observed in hits
targets = sorted({t for (t,_) in best_choice.keys()})
targets_set = set(targets)

# For each target, measure included & missing samples
t_to_samples_present = defaultdict(set)
for (t,s) in best_choice.keys():
    t_to_samples_present[t].add(s)

for t in targets:
    present = t_to_samples_present.get(t, set())
    t_included[t] = len(present)
    missing = [s for s in all_samples if s not in present]
    if missing:
        t_missing[t] = missing
    t_distinct_samples[t] = present

# Write summary with distinct sample count and a comma-separated list of missing samples
with open(summary, 'w', newline='') as s:
    w = csv.writer(s, delimiter='\t')
    w.writerow(["target_id",
                "included_samples",               # number of samples with a chosen sequence
                "distinct_samples",               # same as included_samples (explicit)
                "total_samples_expected",         # |all_samples|
                "missing_samples_list",           # comma-separated list
                "excluded_multitarget_queries"])  # total # contigs excluded globally
    excluded_count = len(excluded)
    total_expected = len(all_samples)
    for t in targets:
        inc = t_included.get(t, 0)
        missing_list = ",".join(t_missing.get(t, []))
        w.writerow([t, inc, inc, total_expected, missing_list, excluded_count])
PY

[[ -s "$WARNINGS" ]] && echo ">> Exclusion log:   $WARNINGS" || echo ">> Exclusion log:   (none)"
echo ">> Summary written to: $SUMMARY_TSV"
echo ">> Map (target,sample -> query): $MAP_TSV"

# -------- 4) Write ONE sequence per (target, sample) into each per-target multi-FASTA --------
export COMBINED OUTDIR MAP_TSV
python3 - <<'PY'
import os, csv
from collections import defaultdict

combined = os.environ["COMBINED"]
map_tsv  = os.environ["MAP_TSV"]
outdir   = os.environ["OUTDIR"]
os.makedirs(outdir, exist_ok=True)

# Read mapping: target \t sample \t query
t_s_to_q = {}
t_to_samples = defaultdict(list)
with open(map_tsv) as f:
    r = csv.reader(f, delimiter='\t')
    for row in r:
        if not row or len(row) != 3:
            continue
        t, s, q = row
        t_s_to_q[(t, s)] = q
        t_to_samples[t].append(s)

# Build target -> set(queries) and the set of all needed query IDs
t_to_queries = defaultdict(set)
for (t, s), q in t_s_to_q.items():
    t_to_queries[t].add(q)
needed = set(t_s_to_q.values())

# Load sequences we need from the combined FASTA
seqs = {}
with open(combined) as f:
    hdr = None
    buf = []
    def flush_current():
        if hdr is None:
            return
        tok = hdr.split()[0].lstrip(">")
        if tok in needed:
            seqs[tok] = "".join(buf)
    for line in f:
        if line.startswith(">"):
            flush_current()
            hdr = line.strip()
            buf = []
        else:
            buf.append(line.strip())
    flush_current()

# Open handles per target and write 1 seq per sample (best already chosen earlier)
handles = {}
def geth(t):
    h = handles.get(t)
    if h is None:
        safe = t.replace("/", "_")
        h = open(os.path.join(outdir, f"{safe}.multi.fasta"), "w")
        handles[t] = h
    return h

for t in sorted(t_to_samples):
    h = geth(t)
    for s in sorted(set(t_to_samples[t])):
        q = t_s_to_q[(t, s)]
        seq = seqs.get(q)
        if seq is None:
            # Should not happen; skip just in case
            continue
        h.write(f">{s}|{q}\n{seq}\n")

for h in handles.values():
    h.close()

print(f"Wrote {len(handles)} per-target multi-FASTAs to {outdir}")
PY

echo
echo ">> Done."
echo "Per-target multi-FASTAs: $OUTDIR/<target_id>.multi.fasta"
echo "Summary table:          $SUMMARY_TSV"
echo "Exclusion warnings:     $WARNINGS (if any)"
