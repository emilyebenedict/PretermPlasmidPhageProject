#!/bin/bash
#===============================================================================
# File Name    : s08_pharokka.sh
# Description  : This script will run prokka in parallel
# Usage        : sbatch s08_pharokka.sh
# Author       : Kailun Zhang. kailun@wustl.edu
# Version      : 2
# Modified     : 2022-07-20 
# Created      : 2024-06-27
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=pharokka
#SBATCH --array=2-2071
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm/pharokka/x_pharokka_%a.out
#SBATCH --error=slurm/pharokka/y_pharokka_%a.err

eval $( spack load --sh /buvzc6u ) # miniconda3@4.10.3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate /home/ai.z/.conda/envs/pharroka_1

basedir="/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3"
indir="/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/combinedviruses/splitviruses"
outdir="/scratch/gdlab/p.sri/250113_P4Pipeline/250204_annotation/pharokka"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250113_P4Pipeline/250204_annotation/pharokkamapping.txt`

set -x

time 
python /home/ai.z/.conda/envs/pharroka_1/bin/pharokka.py -i ${indir}/${sample}.fa -o ${outdir}/${sample} -d /home/ai.z/.conda/envs/pharroka_1/databases -t ${SLURM_CPUS_PER_TASK}
python /home/ai.z/.conda/envs/pharroka_1/bin/pharokka_plotter.py -i ${outdir}/${sample}/phanotate.ffn -n "${sample}" -o ${outdir}/${sample} --interval 8000 --annotations 0.5 --plot_title '${sample}' --truncate 50
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi
