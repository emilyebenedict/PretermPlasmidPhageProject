#!/bin/bash
#===============================================================================
#
# File Name    : s19_bowtie_map.sh
# Description  : This script will align reads to the indices from s18
# Usage        : sbatch s_bowtie2_human.sh
# Author       : Jian Ryou, j.ryou@wustl.edu
# Version      : 1.0
# Created On   : Thu Apr  7 10:44:04 CDT 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=6-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=bwa
#SBATCH --array=1-16
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/bwa/z_bwa_align-2_%a.out
#SBATCH --error=slurm_out/bwa/z_bwa_align-2_%a.out

# load bwa and samtools
eval $( spack load --sh bwa ) 
eval $( spack load --sh /6p5wlkk )

basedir="$PWD"
reads_in="${basedir}/d03_filt-out"
assembly_in="${basedir}/d04_flye"
outdir="${basedir}/d06_binning/bwa_align"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/bwa_left.txt`
#index is where the output from bowtie2-build is stored
index="${basedir}/d06_binning/bwa_index"
prefix="${sample}"

#make output directory
mkdir -p ${outdir}

#unzipping files
#gunzip ${reads_in}/${sample}.fastq.gz

set -x
bwa mem -t $SLURM_CPUS_PER_TASK ${index}/${sample} ${reads_in}/${sample}.fastq | samtools sort -o ${outdir}/${sample}_alignment.bam

# index BAM file (needed by Metabat2)
samtools index ${outdir}/${sample}_alignment.bam

gzip ${reads_in}/${sample}.fastq

RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured in ${sample}!"
    exit $RC
fi
