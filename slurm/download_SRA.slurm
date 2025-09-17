#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J download_SRA
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-4



set -euo pipefail



# --- Configuration ---------------------------------------------------------
# Either: 1) edit the SRR_LIST below and put your four SRR accessions there
#         2) or create a plain text file named `srr_list.txt` (one SRR per line)
#            in the same directory as this script; it will be used instead.
SRR_file='adult_aged_kidney_ATAC/SRR_Acc_List_2.txt'
mapfile -t SRR_LIST < "$SRR_file"

# Path to the unpacked SRA Toolkit (adjust if it's located elsewhere).
# This assumes the toolkit directory is in the directory you submit the job from.
SRATOOLKIT_DIR="sratoolkit.3.2.1-alma_linux64"

# Output directory for FASTQ files
OUTDIR="adult_aged_kidney_snRNAseq"

# Number of threads to pass to fasterq-dump; defaults to cpus-per-task if set
THREADS=${SLURM_CPUS_PER_TASK:-4}

# ---------------------------------------------------------------------------

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2

# Select accession for this array task (SLURM_ARRAY_TASK_ID starts at 1)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
SRR=${SRR_LIST[$IDX]:-}
if [ -z "$SRR" ]; then
	echo "ERROR: No SRR set for array index $SLURM_ARRAY_TASK_ID" >&2
	exit 2
fi

mkdir -p "$OUTDIR"

# Ensure toolkit binaries are on PATH
if [ -d "$SRATOOLKIT_DIR/bin" ]; then
	export PATH="$SRATOOLKIT_DIR/bin:$PATH"
else
	echo "WARNING: SRA Toolkit not found at $SRATOOLKIT_DIR - make sure the path is correct" >&2
fi

echo "Starting download for $SRR -> $OUTDIR"

# Use prefetch to fetch the .sra file into OUTDIR, then convert to FASTQ with fasterq-dump
prefetch "$SRR" -O "$OUTDIR" || { echo "prefetch failed for $SRR" >&2; exit 3; }

SRA_DIR="$OUTDIR/$SRR"
SRA_FILE="$SRA_DIR/${SRR}.sralite"

# Convert to FASTQ (split paired files if present) and write to OUTDIR
fasterq-dump "$SRA_FILE" -O "$SRA_DIR" -e "$THREADS" --split-files || { echo "fasterq-dump failed for $SRA_FILE" >&2; exit 5; }

# Gzip FASTQ outputs to save space
#gzip -f "$SRA_DIR"/*.fastq || echo "Note: gzip failed (maybe no .fastq files present)" >&2

# Optional: remove the .sra to save space
rm -f "$SRA_FILE" || true

echo "Finished $SRR"


