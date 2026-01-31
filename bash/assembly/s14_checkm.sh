#!/bin/bash
#===============================================================================
#
# File Name    : s14_checkm.sh
# Description  : This script will assess the quality of genome assemblies.
# Usage        : sbatch s14_checkm.sh
# Author       : Jian Ryou
# Version      : 1.0
# Created On   : Wed Aug 17 17:37:20 CST 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkm
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --cpus-per-task=8
#SBATCH --mem=42G
#SBATCH --output=slurm_out/checkm/x_checkm_%a.out

eval $( spack load --sh py-checkm-genome )

basedir="$PWD"
indir="${basedir}/d04_flye"
tmpdir="${basedir}/d05_QC/tmp_flye_out"
outdir="${basedir}/d05_QC/checkm"

#make output directory
mkdir -p ${outdir}

#read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/barcodes.txt`

# IMPORTANT: the find command is meant to find the output dirs from flye.
mkdir -p ${tmpdir}
for dir in `find ${indir} -type d -maxdepth 1 -mindepth 1`; do cp $dir/assembly.fasta $tmpdir/$(basename $dir).fasta; done;

set -x

time checkm lineage_wf -f ${outdir}/231102_checkm_output.txt\
        -t ${SLURM_CPUS_PER_TASK} \
        -x fasta \
        --tab_table \
        ${tmpdir} ${outdir}

RC=$?

set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
