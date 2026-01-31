#!/bin/bash
#===============================================================================
# File Name    : s07_seqtk.sh
# Description  : This script will subsample 2M reads from all samples
#===============================================================================
#SBATCH --job-name=seqtk
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-17
#SBATCH --output=slurm_out/seqtk/x_seqtk_%a.out
#SBATCH --error=slurm_out/seqtk/y_seqtk_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $(spack load --sh seqtk@1.3)

basedir="$PWD"
indir="${basedir}/d01_trimmedreads"
outdir="${indir}/d02_subsampled"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/p4_samples.txt`

# IMPORTANT: random seed (designated by '-s') must be the same for FW/RV reads
set -x

seqtk sample -s7 ${indir}/${sample}_FW_clean.fastq 2000000 > ${outdir}/${sample}_FW_clean.fastq | seqtk sample -s7 ${indir}/${sample}_RV_clean.fastq 2000000 > ${outdir}/${sample}_RV_clean.fastq

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in ${sample}!"
  exit $RC
fi
