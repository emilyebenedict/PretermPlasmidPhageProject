#!/bin/bash

#===============================================================================
# Name         : s25_bbrepair.sh
# Description  : This script will run bbtools repair to repair deconseq reads
# Usage        : sbatch s25_bbrepair.sh
# Author       : Luke Diorio-Toth, modified by Jian Ryou j.ryou@wustl.edu
# Version      : 1.0
# Created On   : 2019_01_21
#===============================================================================


#Submission script for HTCF
#SBATCH --time=0-6:00:00 # days-hh:mm:ss
#SBATCH --job-name=bbtools
#SBATCH --array=1-267
#SBATCH --mem=64G
#SBATCH --output=slurm_out/bbtools/x_fix_%a.out
#SBATCH --error=slurm_out/bbtools/y_fix_%a.err

eval $(spack load --sh bbmap@38.63)

basedir="$PWD"
indir="${basedir}/d03_preterm_deconseq"
outdir="${basedir}/d03_preterm_deconseq/repaired"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/drewg_preterm.txt`


repair.sh ow=t in=${indir}/${sample}_fwd_clean.fq \
                in2=${indir}/${sample}_rev_clean.fq \
                out=${outdir}/${sample}_FW_clean.fastq \
                out2=${outdir}/${sample}_RV_clean.fastq \
                outs=${outdir}/${sample}_singletons_CLEAN.fastq \
                repair=t tossbrokenreads=t 

gzip -f ${outdir}/${sample}_FW_clean.fastq
gzip -f ${outdir}/${sample}_RV_clean.fastq
gzip -f ${indir}/${sample}_rev_clean.fq
gzip -f ${indir}/${sample}_fwd_clean.fq

