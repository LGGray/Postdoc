#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J spatial_xist_panels
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00

# Xist expression panels for the annotated Visium HD sections.
#
#   sbatch slurm/spatial_xist_panels.slurm              # 9w and 78w
#   sbatch slurm/spatial_xist_panels.slurm 78w          # one sample
#   sbatch --export=ALL,GENE=Tsix slurm/spatial_xist_panels.slurm
#
# Positional, not environment: --export=NONE above means `GENE=... sbatch` does
# not reach the job, hence the --export=ALL,... form for the overrides.
#
# Memory: one ~400 MB annotated object is read at a time and reduced to seven
# columns before the next, so the peak is a single object plus the plotting
# frame. 64G is generous; 32G would do. Runtime is a few minutes, almost all of
# it rendering the two ~240k-bin maps.

module load slurm_setup

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then SAMPLES=(9w 78w); fi

BASE=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
SPATIAL=$BASE/adult_aged_spatial
BIN=${BIN:-square_008um}
GENE=${GENE:-Xist}
cd "$SPATIAL" || exit 1

# Fail before the queue time is spent rather than after loading an object.
for s in "${SAMPLES[@]}"; do
  f=$SPATIAL/annotation/${s}_${BIN}/${s}_${BIN}_annotated.rds
  if [ ! -s "$f" ]; then
    echo "Missing $f - run spatial_cell_annotation.slurm $s first" >&2
    exit 1
  fi
done

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
# seurat_env, as for spatial_cell_annotation.slurm: the annotated objects are
# Seurat v5 and readRDS needs the same class definitions that wrote them.
conda activate seurat_env

# Bash arrays do not survive into R; join with commas, which the script splits.
SAMPLE_LIST=$(IFS=,; echo "${SAMPLES[*]}")

SPATIAL_DIR="$SPATIAL" BIN="$BIN" GENE="$GENE" SAMPLES="$SAMPLE_LIST" \
  Rscript "$BASE/Postdoc/spatial/spatial_xist_panels.R"
rc=$?

conda deactivate
exit $rc
