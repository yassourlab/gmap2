#!/bin/bash
# =============================================================================
# HUMAnN3 Pipeline — SLURM array job (one task per sample)
#
# Processes paired-end shotgun metagenomics FASTQ files through HUMAnN3:
#   1. Concatenates R1 + R2 reads (standard HUMAnN3 input)
#   2. Runs HUMAnN3 (taxonomic + functional profiling)
#   3. Normalizes outputs to relative abundance
#
# After all array tasks complete, run the companion script:
#   merge_humann_tables.sh  — joins all per-sample tables into cohort tables
#
# Setup:
#   1. Fill in the CONFIG section below.
#   2. Create a sample list file: one sample name per line (no path, no suffix).
#      Example:  echo -e "sample1\nsample2\nsample3" > samples.txt
#   3. Submit:
#      sbatch --array=0-<N-1>%<MAX_PARALLEL> run_humann_pipeline.sh
#      where N = number of lines in SAMPLE_LIST, e.g. --array=0-49%10
# =============================================================================

#SBATCH --job-name=humann
#SBATCH --time=48:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/humann_%A_%a.out
#SBATCH --error=logs/humann_%A_%a.err

# =============================================================================
# CONFIG — set PROJECT_DIR, then adjust the other paths if needed
# =============================================================================

PROJECT_DIR=   # <-- fill it with path to the dir to the samples

CONDA_INIT="/path/to/miniconda3/etc/profile.d/conda.sh"   # path to conda init script
CONDA_ENV="humann"                                          # conda environment with HUMAnN3

FASTQ_DIR="$PROJECT_DIR/raw_fastq"          # directory containing paired FASTQ files
SAMPLE_LIST="$PROJECT_DIR/samples.txt"       # one sample name per line
OUTPUT_DIR="$PROJECT_DIR/humann_output"      # per-sample output directories go here
SCRATCH_DIR="/tmp/humann_scratch"            # temp dir for concatenated reads

# HUMAnN3 reference databases (these are typically installed system-wide)
CHOCOPHLAN_DB="/path/to/chocophlan"      # nucleotide database (full or DEMO)
UNIREF_DB="/path/to/uniref"             # protein database (e.g. uniref90_diamond)
METAPHLAN_DB="/path/to/metaphlan_db"    # MetaPhlAn bowtie2 index directory

# FASTQ filename pattern: <sample><R1_SUFFIX> / <sample><R2_SUFFIX>
R1_SUFFIX="_R1.fastq.gz"
R2_SUFFIX="_R2.fastq.gz"

THREADS=8   # match --cpus-per-task above

# =============================================================================
# SETUP
# =============================================================================

mkdir -p logs

source "$CONDA_INIT"
conda activate "$CONDA_ENV"

echo "HUMAnN version: $(humann --version 2>&1)"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")

if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "============================================"
echo "Sample:     $SAMPLE"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Start:      $(date)"
echo "Host:       $(hostname)"
echo "============================================"

SAMPLE_OUT="$OUTPUT_DIR/$SAMPLE"
SAMPLE_TMP="$SCRATCH_DIR/$SAMPLE"
mkdir -p "$SAMPLE_OUT" "$SAMPLE_TMP"

# =============================================================================
# STEP 1 — Concatenate R1 + R2
# HUMAnN3 takes a single input file; concatenating R1+R2 is standard practice.
# =============================================================================

R1="$FASTQ_DIR/${SAMPLE}${R1_SUFFIX}"
R2="$FASTQ_DIR/${SAMPLE}${R2_SUFFIX}"
CONCAT="$SAMPLE_TMP/${SAMPLE}_concat.fastq.gz"

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: FASTQ files not found:"
    echo "  R1: $R1"
    echo "  R2: $R2"
    exit 1
fi

echo "[$(date '+%H:%M:%S')] Concatenating R1 + R2..."
cat "$R1" "$R2" > "$CONCAT"

# =============================================================================
# STEP 2 — Run HUMAnN3
# =============================================================================

echo "[$(date '+%H:%M:%S')] Running HUMAnN3..."

humann \
    --input "$CONCAT" \
    --output "$SAMPLE_OUT" \
    --threads "$THREADS" \
    --nucleotide-database "$CHOCOPHLAN_DB" \
    --protein-database "$UNIREF_DB" \
    --metaphlan-options "--bowtie2db $METAPHLAN_DB" \
    --output-basename "$SAMPLE" \
    --remove-temp-output \
    --verbose

if [ $? -ne 0 ]; then
    echo "ERROR: HUMAnN3 failed for sample $SAMPLE"
    exit 1
fi

# =============================================================================
# STEP 3 — Normalize to relative abundance
# =============================================================================

echo "[$(date '+%H:%M:%S')] Normalizing outputs..."

humann_renorm_table \
    --input  "$SAMPLE_OUT/${SAMPLE}_pathabundance.tsv" \
    --output "$SAMPLE_OUT/${SAMPLE}_pathabundance_relab.tsv" \
    --units relab --special n

humann_renorm_table \
    --input  "$SAMPLE_OUT/${SAMPLE}_genefamilies.tsv" \
    --output "$SAMPLE_OUT/${SAMPLE}_genefamilies_relab.tsv" \
    --units relab --special n

# =============================================================================
# CLEANUP
# =============================================================================

rm -rf "$SAMPLE_TMP"

echo "============================================"
echo "Done: $SAMPLE  |  $(date)"
echo "Output: $SAMPLE_OUT"
echo "============================================"
