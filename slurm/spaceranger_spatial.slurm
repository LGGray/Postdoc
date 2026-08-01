#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spaceranger_HD
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%A_%a.%N.out
#SBATCH -e out_error/%x.%A_%a.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_tiny
#SBATCH --qos=cm4_tiny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-2

# Visium HD 3' (poly-A capture, no probe set) — Project_1618_Eckstein_Andergassen_Spatial
# Two CytAssist capture areas: adult (9w) and aged (78w) heart.
# Array task 1 = 9w, task 2 = 78w.

module load slurm_setup

export PATH=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/spaceranger-4.1.0:$PATH

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
DATA_DIR=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/andergassen_lab/00_raw_data/Project_1618_Eckstein_Andergassen_Spatial
HISTO_DIR=$DATA_DIR/00_histology
OUT_DIR=$BASE/adult_aged_spatial

# Standard 10x reference for now. For the allelic run, swap to the N-masked
# GRCm39 reference (built with SNPsplit + spaceranger mkref) and re-run.
REF=$BASE/refdata-gex-GRCm39-2024-A
# REF=$BASE/GRCm39/refdata-gex-GRCm39-2024-A-Nmasked-CAST

##### ---------- FILL THESE IN BEFORE SUBMITTING ---------- #####
# Visium slide serial (e.g. H1-XXXXXXX) and capture area (A1/B1/C1/D1).
SLIDE_9W="H1-8H6HPK4"
AREA_9W="A1"
SLIDE_78W="H1-8H6HPK4"
AREA_78W="D1"

# CytAssist 3000x3000 TIFF for each capture area.
CYTA_DIR=$HISTO_DIR/CytAssist_Images
CYTA_9W=$CYTA_DIR/CAVG10444_2025-07-30_12-03-35_2025-07-30_10-46-49_H1-8H6HPK4_A1_Hoelzl_3HD_r1_20250730_19_4106A.tif
CYTA_78W=$CYTA_DIR/CAVG10444_2025-07-30_12-03-35_2025-07-30_10-46-49_H1-8H6HPK4_D1_Hoelzl_3HD_r1_20250730_21_0561A.tif

# High-resolution brightfield H&E taken before CytAssist.
# Confirmed by CytAssist filenames: A1 = 4106A, D1 = 0561A.
IMAGE_9W=$HISTO_DIR/Hoelzl_vorCytassist_HE_RUN1_20250730_DAN4106.tif
IMAGE_78W=$HISTO_DIR/Hoelzl_vorCytassist_HE_RUN1_20250730_DAN0561.tif

# Manual fiducial/image registration exported from Loupe Browser.
LOUPE_9W=$CYTA_DIR/H1-8H6HPK4-A1-DAN-4106-fiducials-image-registration.json
LOUPE_78W=$CYTA_DIR/H1-8H6HPK4-D1-DAN-0561-fiducials-image-registration.json

# Confirm 25L006362 = 9w and 25L006363 = 78w against Project_1618_lims.csv.
SAMPLE_9W=25L006362
SAMPLE_78W=25L006363
##### -------------------------------------------------------- #####

# Extra bin size on top of the default 2/8/16 um outputs. Must be an even
# number between 4 and 100, and not 8 or 16. 32 um ~ one cardiomyocyte width.
CUSTOM_BIN=32

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
  ID=9w
  SAMPLE=$SAMPLE_9W
  FASTQS=$DATA_DIR/Sample_$SAMPLE_9W
  CYTAIMAGE=$CYTA_9W
  IMAGE=$IMAGE_9W
  LOUPE=$LOUPE_9W
  SLIDE=$SLIDE_9W
  AREA=$AREA_9W
elif [ "$SLURM_ARRAY_TASK_ID" -eq 2 ]; then
  ID=78w
  SAMPLE=$SAMPLE_78W
  FASTQS=$DATA_DIR/Sample_$SAMPLE_78W
  CYTAIMAGE=$CYTA_78W
  IMAGE=$IMAGE_78W
  LOUPE=$LOUPE_78W
  SLIDE=$SLIDE_78W
  AREA=$AREA_78W
else
  echo "Unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Fail before burning 24 h on a typo.
for d in "$REF" "$FASTQS"; do
  if [ ! -d "$d" ]; then echo "Missing input directory: $d" >&2; exit 1; fi
done
# These must be files, not directories - spaceranger only errors out after
# the pipeline has already started.
for f in "$CYTAIMAGE" "$IMAGE" "$LOUPE"; do
  if [ ! -f "$f" ]; then echo "Not a file: $f" >&2; exit 1; fi
done
if [ -z "$SLIDE" ] || [ -z "$AREA" ]; then
  echo "SLIDE/AREA not set for $ID - fill them in or switch to --unknown-slide" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR" || exit 1

# spaceranger wants a little headroom below the cgroup limit.
LOCALCORES=$((SLURM_CPUS_PER_TASK - 2))
LOCALMEM=128

spaceranger count \
  --id="$ID" \
  --transcriptome="$REF" \
  --fastqs="$FASTQS" \
  --sample="$SAMPLE" \
  --cytaimage="$CYTAIMAGE" \
  --image="$IMAGE" \
  --slide="$SLIDE" \
  --area="$AREA" \
  --loupe-alignment="$LOUPE" \
  --custom-bin-size="$CUSTOM_BIN" \
  --create-bam=true \
  --localcores="$LOCALCORES" \
  --localmem="$LOCALMEM"

