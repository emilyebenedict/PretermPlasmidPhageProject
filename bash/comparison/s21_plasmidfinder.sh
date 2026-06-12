#!/bin/bash
#===============================================================================
# File Name    : plasmidfinder.sh
# Description  : This script will run plasmidfinder in parallel
# Usage        : sbatch plasmidfinder.sh
# Author       : Lindsey Hall
# Created      : 23720
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=plasmidfinder
#SBATCH --array=1-102%20
#SBATCH --mem=2G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/plasmidfinder/x_plasmidfinder_%a.out
#SBATCH --error=slurm_out/plasmidfinder/y_plasmidfinder_%a.err


eval $(spack load --sh plasmidfinder)

basedir="$PWD"
indir="${basedir}/finalized_assemblies"
outdir="${basedir}/d13_plasmidfinder"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/all_seq_list.txt`

set -x

mkdir -p ${outdir}/${sample}

plasmidfinder.py -i ${indir}/${sample}/scaffolds.fasta -o ${outdir}/${sample} -p /ref/gdlab/data/plasmidfinder_db/20220330/ -t .95 -x
