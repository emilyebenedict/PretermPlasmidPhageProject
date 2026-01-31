#!/bin/bash
#===============================================================================
#
# File Name    : s03_dorado_fastq.sh
# Description  : This script will convert demultiplexed dorado .bams to .fastqs
# Usage        : sbatch s03_dorado.sh
# Author       : Emily Benedict, ebenedict@wustl.edu
# Version      : 1.0
#===============================================================================
#
#SBATCH --time=72:00:00
#SBATCH --job-name=dorado_fastq
#SBATCH --array=1-97%12
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --output=slurm_out/dorado_fastq/dorado_fastq_%A.out
#SBATCH --error=slurm_out/dorado_fastq/dorado_fastq_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $( spack load --sh samtools )

basedir="$PWD"
indir="${basedir}/d01_dorado/demux_files"
outdir="${basedir}/d01_dorado/fastqs"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/dorado_demux_bams.txt`

samtools bam2fq ${indir}/${sample}.bam > ${outdir}/${sample}.fastq
