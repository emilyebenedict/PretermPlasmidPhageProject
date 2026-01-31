#!/bin/bash
#===============================================================================
#
# File Name    : s01_dorado_merge.sh
# Description  : This script will merge dorado output, step 2
# Usage        : sbatch s01_dorado_merge.sh
# Author     : Emily Benedict, ebenedict@wustl.edu
# Version      : 1.0
#===============================================================================
#
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --output=slurm_out/dorado_merge/dorado_%A.out
#SBATCH --error=slurm_out/dorado_merge/dorado_%A.err
#SBATCH --time=72:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

basedir="$PWD"
dir="${basedir}/d01_dorado"
name="p4_round3_dorado"
indir="${dir}/calling_output"

eval $( spack load --sh samtools )

samtools merge ${dir}/${name}.bam ${indir}/*.bam
