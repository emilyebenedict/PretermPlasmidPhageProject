#!/bin/bash
#==============================================================================
#
# File Name    : s04_filtlong.sh
# Description  : Runs filtlong
# Usage        : sbatch s04_filtlong.sh
# Author       : Mark Gorelik, gmark@wustl.edu
# Version      : 1.0
# Created On   : 2022-10-11
# Last Modified: 2023-11-03 by Emily Benedict, ebenedict@wustl.edu
#==============================================================================
#
#Submission script for HTCF
#SBATCH --time=3-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=flitlong
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=slurm_out/filtlong/s02_filtlong%a.out
#SBATCH --error=slurm_out/filtlong/s02_filtlong%a.err
#SBATCH --array=1-97%12
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --mail-type=END

eval $( spack load --sh filtlong )

basedir="$PWD"
input_dir="${basedir}/d01_dorado/fastqs"
output_dir="${basedir}/d02_filtlong"

mkdir -p ${output_dir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/dorado_demux_bams.txt`

filtlong --min_length 2000 --keep_percent 90 ${input_dir}/${sample}.fastq | gzip > ${output_dir}/${sample}.fastq.gz

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
