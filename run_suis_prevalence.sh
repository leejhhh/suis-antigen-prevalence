#!/bin/bash
set -e

# --- Configuration ---
QUERY_PROTEIN_FASTA="query_antigens.fasta"
OUTPUT_DIR="suis_prevalence_analysis"
FASTA_DIR="suis_selected"       # 이미 준비된 FASTA 모음
PYTHON_SCRIPT="parse_prevalence.py"
EVALUE="1e-5"
THREADS=$(nproc)
# Allow runtime override of identity/coverage thresholds
: ${MIN_IDENTITY:=70.0}   # default 70%
: ${MIN_COVERAGE:=0.8}    # default 80% (fraction)

# --- Input Validation ---
if [ ! -f "$QUERY_PROTEIN_FASTA" ]; then
  echo "Error: Query FASTA ($QUERY_PROTEIN_FASTA) not found."
  exit 1
fi
if ! command -v makeblastdb &> /dev/null || ! command -v tblastn &> /dev/null; then
  echo "Error: BLAST+ tools (makeblastdb, tblastn) not found in PATH."
  exit 1
fi
if [ ! -d "$FASTA_DIR" ] || [ -z "$(ls -A "$FASTA_DIR"/*.fna 2>/dev/null)" ]; then
  echo "Error: FASTA directory '$FASTA_DIR' not found or empty."
  exit 1
fi

echo "--- Starting S. suis Prevalence Analysis ---"
echo "[1/5] Setting up directories..."
mkdir -p "$OUTPUT_DIR"

# --- 2. Merge FASTA files ---
echo "[2/5] Merging FASTA files from '$FASTA_DIR'..."
MERGED_FASTA="${OUTPUT_DIR}/all_suis_genomes.fna"
find "$FASTA_DIR" -name '*.fna' -exec cat {} + > "$MERGED_FASTA"
echo "Merged into: $MERGED_FASTA"

# --- 3. Count genomes ---
TOTAL_GENOMES=$(ls -1 "$FASTA_DIR"/*.fna | wc -l)
echo "[3/5] Total genomes (FASTA files): $TOTAL_GENOMES"

# --- 4. Create BLAST DB ---
echo "[4/5] Creating BLAST database..."
BLAST_DB_NAME="${OUTPUT_DIR}/suis_db"
makeblastdb -in "$MERGED_FASTA" -dbtype nucl -out "$BLAST_DB_NAME" -parse_seqids -title "S.suis_DB"
echo "DB: $BLAST_DB_NAME.*"

# --- 5. Run tblastn & parse ---
echo "[5/5] Running tblastn & parsing results..."
BLAST_OUT="${OUTPUT_DIR}/blast_results.tsv"
tblastn -query "$QUERY_PROTEIN_FASTA" \
        -db "$BLAST_DB_NAME" \
        -evalue "$EVALUE" \
        -outfmt 6 \
        -num_threads "$THREADS" \
        -out "$BLAST_OUT"
echo "tblastn done: $BLAST_OUT"

PARSER_OUT="${OUTPUT_DIR}/genomes_with_hit_stats.tsv"
python3 "$PYTHON_SCRIPT" \
  -i "$BLAST_OUT" \
  -t "$TOTAL_GENOMES" \
  -q "$QUERY_PROTEIN_FASTA" \
  --min_identity "$MIN_IDENTITY" \
  --min_coverage "$MIN_COVERAGE" \
  -o "$PARSER_OUT"
echo "Parsed stats: $PARSER_OUT"

echo "--- Analysis Complete ---"
