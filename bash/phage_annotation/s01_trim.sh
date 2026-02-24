#!/bin/bash
#===================================================================================================
# File Name    : trimM_virome.sh
# Description  : Trim reads for virome analysis using Trimmomatic
# Usage        : sbatch trimM_virome.sh
# Author       : Yunru Peng <p.yunru@wustl.edu>
# Created On   : April 2 2025
# Last Modified: 
#===================================================================================================
# Submission script for HTCF 
#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name=trimM
#SBATCH --output=slurm/01_triminduction/x_trimVirome_%A.out
#SBATCH --error=slurm/01_triminduction/y_trimVirome_%A.err

eval $(spack load --sh trimmomatic@0.39)
# Set up trimmomatic parameters
TRIMMOMATIC_HOME="/ref/gdlab/software/spack-0.18.0/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/trimmomatic-0.39-33pr7ksvkjm2ddshc3po5lmrlwxjdq5n/bin"
adapt="/ref/gdlab/software/spack-0.18.0/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/trimmomatic-0.39-33pr7ksvkjm2ddshc3po5lmrlwxjdq5n/share/adapters/NexteraPE-PE.fa"
export JAVA_ARGS="-Xmx8000M"

basedir="$PWD"
# Read in a csv file of prefixes of paths to the raw FASTQs and sample identifiers
row=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "/scratch/gdlab/p.sri/250414_HVP_induction/samplenames.csv")
IFS=',' read -r prefix SID <<< "$row"
# Make the output directory
outdir="${basedir}/d01a_trimmomatic"
mkdir -p "${outdir}"
QCdir="${basedir}/dxx_QC/${SID}/trimmomatic"
mkdir -p "${QCdir}"

set -x
# Parameters are based on PMID {33612831}
# [ILLUMINACLIP:2:30:10:2:True] allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases). The minimum ad$
# [SLIDINGWINDOW:4:30] scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 30
java -jar $TRIMMOMATIC_HOME/trimmomatic-0.39.jar PE -threads ${SLURM_CPUS_PER_TASK} -phred33 \
  ${prefix}*R1*.gz ${prefix}*R2*.gz \
  ${outdir}/${SID}_FW_trim.fastq.gz ${outdir}/${SID}_FW_trim_UP.fastq.gz \
  ${outdir}/${SID}_RV_trim.fastq.gz ${outdir}/${SID}_RV_trim_UP.fastq.gz \
  ILLUMINACLIP:${adapt}:2:30:10:2:True \
  SLIDINGWINDOW:4:30 \
MINLEN:50 2> ${QCdir}/${SID}_trim.log

set +x

if [ -e ${outdir}/${SID}_RV_trim_UP.fastq.gz ]
then
  echo "[Success] trimmomatic FASTQ array# ${SLURM_ARRAY_TASK_ID}"
else
  echo "[Error] trimmomatic FASTQ array# ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

cat ${outdir}/${SID}_FW_trim_UP.fastq.gz ${outdir}/${SID}_RV_trim_UP.fastq.gz > ${outdir}/${SID}_UP_trim.fastq.gz

if [ -e ${outdir}/${SID}_UP_trim.fastq.gz ]
then
  echo "[Success] zcat UP FASTQs array# ${SLURM_ARRAY_TASK_ID}"
else
  echo "[Error] zcat UP FASTQs array# ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

set -x
rm ${outdir}/${SID}_FW_trim_UP.fastq.gz ${outdir}/${SID}_RV_trim_UP.fastq.gz
set +x
