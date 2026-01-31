#!/bin/bash
#===============================================================================
#
# File Name    : s07_quast.sh
# Description  : This script will assess the quality of genome assemblies.
# Usage        : sbatch s07_quast.sh
# Author       : Jian Ryou
# Version      : 1.0
# Created On   : Wed Aug 18 12:07:20 CST 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=quast
#SBATCH --array=1-71
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=slurm_out/quast/x_quast_%a.out

eval $( spack load --sh py-quast@5.2.0 )

basedir="$PWD"
indir="${basedir}/d06_QC/tmp_out"
outdir="${basedir}/d06_QC/quast"

#make output directory
mkdir -p ${outdir}

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/barcodes.txt`

set -x

time quast.py ${indir}/${sample}.fasta \
        -o ${outdir}/${sample} \

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
