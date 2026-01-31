#!/bin/bash
#===============================================================================
#
# File Name    : s00_dorado_calling.sh
# Description  : This script will run the GPU dorado basecaller
# Usage        : sbatch testing_dorado_calling.sh
# Author       : Mark Gorelik, modified by Emily Benedict ebenedict@wustl.edu
# Version      : 1.0
#===============================================================================
#
#SBATCH --time=72:00:00
#SBATCH -p gpu
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:2
#SBATCH --array=1-2731%50
#SBATCH --mem=32G
#SBATCH --output=slurm_out/dorado_call/dorado_%A.out
#SBATCH --error=slurm_out/dorado_call/dorado_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $( spack load --sh dorado@0.4.1 )

basedir="$PWD"
indir="${basedir}/d00_raw_lr"
outdir="${basedir}/d01_dorado/calling_output"
doradodir="/scratch/gdlab/ebenedict/preterm_plasmidome"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/round3_pod5_files.txt`

mkdir -p ${outdir}

dorado basecaller ${doradodir}/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 ${indir}/${sample}.pod5 --kit-name "SQK-NBD114-96" > ${outdir}/${sample}.bam 
