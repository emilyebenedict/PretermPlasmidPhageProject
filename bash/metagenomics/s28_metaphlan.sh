#!/bin/bash

#===============================================================================
# File Name    : metaphlan.sh
# Description  : Profiles taxonomic composition of a read set using MetaPhlAn
# Usage        : sbatch metaphlan.sh
# Author       : Luke Diorio-Toth, modified by Jian Ryou j.ryou@wustl.edu
# Version      : 1.2
# Created On   : Wed Jan  6 11:49:47 CST 2021
#===============================================================================

#SBATCH --job-name=metaphlan
#SBATCH --array=1-267%40
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --mem=32G
#SBATCH --output=slurm_out/metaphlan/z_metaphlan_%a_%A.out
#SBATCH --error=slurm_out/metaphlan/z_metaphlan_%a_%A.out

# There are multiple versions of metaphlan installed, which use different dbs
# run "spack find py-metaphlan" to see all of them
eval $( spack load --sh py-metaphlan )

basedir="$PWD"
indir="${basedir}/d03_preterm_deconseq/repaired"
outdir="${basedir}/d04_preterm_metaphlan"
mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/drewg_preterm.txt`

set -x
time metaphlan \
  ${indir}/${sample}_FW_clean.fastq.gz,${indir}/${sample}_RV_clean.fastq.gz \
  --input_type fastq \
  --bowtie2out ${outdir}/${sample}.bowtie2.bz2 \
  -o ${outdir}/${sample}_profile.txt \
  --nproc ${SLURM_CPUS_PER_TASK} --bowtie2db /ref/gdlab/data/metaphlan_db_Oct2022 \
  -t rel_ab_w_read_stats
