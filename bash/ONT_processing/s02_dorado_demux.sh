#!/bin/bash
#===============================================================================
#
# File Name    : s02_dorado_demux.sh
# Description  : This script will demultiplex output from dorado merge, step 3
# Usage        : sbatch s02_dorado_demux.sh
# Author       : Emily Benedict, ebenedict@wustl.edu
# Version      : 1.0
#===============================================================================
#
#SBATCH --time=72:00:00
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --output=slurm_out/dorado_demux/dorado_%A.out
#SBATCH --error=slurm_out/dorado_demux/dorado_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $( spack load --sh dorado@0.4.1 )

basedir="$PWD"
indir="${basedir}/d01_dorado"
outdir="${indir}/demux_files"
name="p4_round3_dorado"

mkdir -p ${outdir}


dorado demux --output-dir ${outdir} ${indir}/${name}.bam --kit-name SQK-NBD114-96
