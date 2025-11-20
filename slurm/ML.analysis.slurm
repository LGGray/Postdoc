#!/bin/bash
#SBATCH --mail-user=lachlan.gray@tum.de
#SBATCH --mail-type=END,FAIL
#SBATCH -J ML.analysis_complete
#SBATCH -D /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2
#SBATCH -o out_error/%x.%j.%N.out
#SBATCH -e out_error/%x.%j.%N.err
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1   
#SBATCH --cpus-per-task=18
#SBATCH --mem=20G
#SBATCH --time=8:00:00
#SBATCH --array=1-48

module load slurm_setup

## Load conda environment ##
source /dss/lrzsys/sys/spack/release/24.4.0/opt/x86_64/miniconda3/24.7.1-gcc-t6x7erm/etc/profile.d/conda.sh
conda activate SLE

cd /dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/SLE_paper/SLE

celltype=$(sed -n "${SLURM_ARRAY_TASK_ID}p" filelist.txt)
# celltype=$(sed -n "${SLURM_ARRAY_TASK_ID}p" predicting.list.txt)

script_dir=/dss/dssfs03/tumdss/pn72lo/pn72lo-dss-0010/go93qiw2/Postdoc/SLE_paper

for i in {1..10}; do
# Feature selection
    python3 $script_dir/boruta_enet.py $celltype $i
done

for i in {1..10}; do
# Logistic regression
  python3 $script_dir/logit.regression.py $celltype $i combined
# Random Forest
  python3 $script_dir/random.forest.py $celltype $i combined
# GBM
  python3 $script_dir/gradient.boosting.machine.py $celltype $i combined
# SVM
  python3 $script_dir/support.vector.machine.py $celltype $i combined
# MLP
  python3 $script_dir/multilayer.perceptron.py $celltype $i combined
# Voting classifier
  python3 $script_dir/voting.classifier.py $celltype $i combined
# Plot all curves
  python3 $script_dir/plot.all.curves.py $celltype $i combined
# print completion message
  echo "Completed run $i for cell type $celltype"
done

for i in {1..10}; do
# Logistic regression
  python3 $script_dir/logit.regression.py $celltype $i intersection
# Random Forest
  python3 $script_dir/random.forest.py $celltype $i intersection
# GBM
  python3 $script_dir/gradient.boosting.machine.py $celltype $i intersection
# SVM
  python3 $script_dir/support.vector.machine.py $celltype $i intersection
# MLP
  python3 $script_dir/multilayer.perceptron.py $celltype $i intersection
# Voting classifier
  python3 $script_dir/voting.classifier.py $celltype $i intersection
# Plot all curves
  python3 $script_dir/plot.all.curves.py $celltype $i intersection
# print completion message
  echo "Completed run $i for cell type $celltype"
done

for i in {1..10}; do
# Logistic regression
  python3 $script_dir/logit.regression.py $celltype $i boruta
# Random Forest
  python3 $script_dir/random.forest.py $celltype $i boruta
# GBM
  python3 $script_dir/gradient.boosting.machine.py $celltype $i boruta
# SVM
  python3 $script_dir/support.vector.machine.py $celltype $i boruta
# MLP
  python3 $script_dir/multilayer.perceptron.py $celltype $i boruta
# Voting classifier
  python3 $script_dir/voting.classifier.py $celltype $i boruta
# Plot all curves
  python3 $script_dir/plot.all.curves.py $celltype $i boruta
# print completion message
  echo "Completed run $i for cell type $celltype"
done

for i in {1..10}; do
# Logistic regression
  python3 $script_dir/logit.regression.py $celltype $i enet
# Random Forest
  python3 $script_dir/random.forest.py $celltype $i enet
# GBM
  python3 $script_dir/gradient.boosting.machine.py $celltype $i enet
# SVM
  python3 $script_dir/support.vector.machine.py $celltype $i enet
# MLP
  python3 $script_dir/multilayer.perceptron.py $celltype $i enet
# Voting classifier
  python3 $script_dir/voting.classifier.py $celltype $i enet
# Plot all curves
  python3 $script_dir/plot.all.curves.py $celltype $i enet
# print completion message
  echo "Completed run $i for cell type $celltype"
done

# for i in {1..10}; do
# python3 $script_dir/predict_independent_SLE.py $celltype $i child
# python3 $script_dir/predict_independent_SLE.py $celltype $i adult
# done