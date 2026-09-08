#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J multiome_umap
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=4:00:00

# Clustering, provisional cell typing and UMAPs for the multiome data.
# Depends on multiome_qc_barcodes.slurm having run.
#
# Produces, under adult_aged_multiome/figures/:
#   umap_rna_clusters.pdf   merged UMAP by cluster / cell type / sample / split
#   umap_wnn_9w.pdf, _78w   per-sample RNA | ATAC | WNN joint, side by side
#   marker_dotplot.pdf      the established marker.genes.txt panel
#   marker_featureplots.pdf key markers incl. Xist
#   umap_rna_qc.pdf         QC metrics on the embedding
#   cluster_markers.csv, cluster_sample_composition.csv, nucleus_celltype.csv
#   sinto_*_bycelltype.txt  group maps for the per-celltype allelic route
#
# No Signac: the ATAC side reuses cellranger-arc's own LSI, and WNN comes from
# Seurat's FindMultiModalNeighbors. seurat_env, not RNAseq.

module load slurm_setup

source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate seurat_env

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc || exit 1
Rscript multiome/02_cluster_umap.R
