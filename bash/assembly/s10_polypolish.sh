#!/bin/bash

#===============================================================================
# File Name    : s10_polypolish.sh
# Description  : s10_Polishes a long-read assembly with illumina reads
# Usage        : sbatch polypolish.sh
# Author       : Luke Diorio-Toth
# Modified     : Emily Benedict, ebenedict@wustl.edu
# Version      : 1.0
#===============================================================================

#SBATCH --job-name=polypolish
#SBATCH --array=1-96%12
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=slurm_out/polypolish/polypolish_%a_%A.out
#SBATCH --error=slurm_out/polypolish/polypolish_%a_%A.err

eval $( spack load --sh polypolish )

basedir="$PWD"
indir="${basedir}/d04_medaka_flye"
reads="/scratch/gdlab/j.ryou/analysis/Preterm_Plasmidome/d00_pilot_raw/pod5"
outdir="${basedir}/d05__polypolish"

#Using the short reads sample list because we're missing a few but have all the ont ones
short_sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/r1_p4_gtac.txt`

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/round1_barcodes.txt`

mkdir -p ${outdir}/${short_sample}

draft="${indir}/${sample}/assembly.fasta"

set -x
# align short reads


time bwa index ${draft}
# make sure -a option is enabled!
time bwa mem \
       -a \
       -t ${SLURM_CPUS_PER_TASK} \
       ${draft} ${reads}/${short_sample}_FW_clean.fastq.gz \
       > ${outdir}/${short_sample}/FW.sam
time bwa mem \
       -a \
       -t ${SLURM_CPUS_PER_TASK} \
       ${draft} ${reads}/${short_sample}_RV_clean.fastq.gz \
       > ${outdir}/${short_sample}/RV.sam
# exclude incorrect alignments based on insert size (optional, but recommended)
time polypolish_insert_filter.py \
       --in1 ${outdir}/${short_sample}/FW.sam \
       --in2 ${outdir}/${short_sample}/RV.sam \
       --out1 ${outdir}/${short_sample}/FW_filtered.sam \
       --out2 ${outdir}/${short_sample}/RV_filtered.sam
# polish the draft genome
time polypolish \
       ${draft} \
       ${outdir}/${short_sample}/FW_filtered.sam \
       ${outdir}/${short_sample}/RV_filtered.sam \
       > ${outdir}/${short_sample}/${short_sample}_polished.fasta
RC=$?
set +x
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
