#!/bin/bash
#===============================================================================
# File Name    : s24_shortbred_quantify.sh
# Description  : This script will qunatify ARGs using pre-processed shotgun sequencing data. 
# Usage        : sbatch s24_shortbred_quantify.sh
# Author       : Bejan Mahmud, modified by Jian Ryou j.ryou@wustl.edu
# Version      : 1.1
# Created On   : 2022-06-05
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=shortbred
#SBATCH --array=1-139%50
#SBATCH --mem=4G
#SBATCH --output=slurm_out/shortbred_quantify/x_shrtbrd-2_%a.out
#SBATCH --error=slurm_out/shortbred_quantify/y_shrtbrd-2_%a.out

#eval $( spack load --sh shortbred@0.9.4 )
eval $( spack load --sh /buvzc6u ) #load miniconda3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate /ref/gdlab/software/envs/shortbred

# Basedir
basedir="$PWD"
indir="${basedir}/d04_abx-naive_fixed"
outdir="${basedir}/d06_shortbred_quantify"

# Make output directory
mkdir -p ${outdir}

# Read in the sample name from the Mapping file
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/sample_names.txt`

# Start debug mode (will send commands to outfile)
set -x

gunzip ${indir}/${sample}_FW_clean.fastq.gz
gunzip ${indir}/${sample}_RV_clean.fastq.gz

time shortbred_quantify.py --markers /ref/gdlab/data/shortbred_markers/ShortBRED_VF_2017_markers.faa\
 --wgs ${indir}/${sample}_FW_clean.fastq ${indir}/${sample}_RV_clean.fastq \
 --results ${outdir}/${sample}_ShortBRED.txt\
 --tmp tmp_shorbred_reads/${sample}_quantify

gzip ${indir}/${sample}_FW_clean.fastq
gzip ${indir}/${sample}_RV_clean.fastq

# Save error code for command
RC=$?
