#!/bin/bash

#===============================================================================
# File Name    : s24_deconseq.sh
# Description  : This script will remove human contamination from fastq reads
# Usage        : sbatch s24_deconseq.sh
# Author       : Rhiannon Vargas, modified by Jian Ryou j.ryou@wustl.edu
# Version      : 1.0
# Created On   : 2020-04-14
# Last Modified: 2020-07-27
#===============================================================================
# Submission script for HTCF
#SBATCH --job-name=deconseq
#SBATCH --array=1-267%50
#SBATCH --mem=64G
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --output=slurm_out/deconseq/x_deconseq_%a.out
#SBATCH --error=slurm_out/deconseq/y_deconseq_%a.err

eval $( spack load --sh deconseq-standalone ) 

#store the base directory
basedir="$PWD"
#store the input path   
indir="${basedir}/d01_gasparrini_preterm_trimmed"
#store the output path
outdir="${basedir}/d03_preterm_deconseq"

#make the output directory
mkdir -p ${outdir}

#store sample name for this array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/drewg_preterm.txt`

#deconseq command

gunzip ${indir}/${sample}_FW_clean.fastq.gz
time deconseq.pl \
     -f ${indir}/${sample}_FW_clean.fastq \
     -out_dir ${outdir} \
     -id ${sample}_fwd \
     -dbs hs_ref
gzip ${indir}/${sample}_FW_clean.fastq

gunzip ${indir}/${sample}_RV_clean.fastq.gz
time deconseq.pl \
     -f ${indir}/${sample}_RV_clean.fastq \
     -out_dir ${outdir} \
     -id ${sample}_rev \
     -dbs hs_ref
gzip ${indir}/${sample}_RV_clean.fastq
