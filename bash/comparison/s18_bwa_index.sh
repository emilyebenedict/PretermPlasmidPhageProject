#!/bin/bash
#===============================================================================
#
# File Name    : s18_bowtie2_index.sh
# Description  : This script will create an index for a set of fasta files obtained from any source
# Usage        : sbatch s18_bowtie2_index.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Version      : 1.0
# Created On   : Tue Aug  23 CDT 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=bwa
#SBATCH --array=1-71%15
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/bwa/z_bwa_%a.out
#SBATCH --error=slurm_out/bwa/z_bwa_%a.out

eval $( spack load --sh bwa )

basedir="$PWD"
indir="${basedir}/d05_QC/tmp_out"
outdir="${basedir}/d06_binning/bwa_index"

#make output directory
mkdir -p ${outdir}

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/barcodes.txt`

set -x

bwa index -p ${sample} ${indir}/${sample}.fasta ${outdir}

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
