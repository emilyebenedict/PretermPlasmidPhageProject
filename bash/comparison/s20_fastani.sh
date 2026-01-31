#!/bin/bash
#===============================================================================
# File Name    : fastani.sh
# Description  : Fast all vs. all ANI calculation of many genomes
# Usage        : sbatch fastani.sh
# Author       : Luke Diorio-Toth, modified by Emily Benedict ebenedict@wustl.edu
# Version      : 1.2
# Created      : Wed Oct 12 13:30:06 CDT 2022
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=all-fastani
#SBATCH --mem=30G
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm_out/fastani/z_fastani_%A.out
#SBATCH --error=slurm_out/fastani/z_fastani_%A.out
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --mail-type=END

eval $( spack load --sh fastani )

basedir="$PWD"
indir="${basedir}/d08_pilot_r1_r3_r4_plasmids"
outdir="${basedir}/d12_all_plasmids_fastani"

# Make output directory
mkdir -p ${outdir}

# Find assemblies and write full path to a file
# for an all v all analysis, that file will be ref and query
find ${indir} -type f -name "*.fasta" >> ${outdir}/genome_paths.txt


set -x
time fastANI --ql ${outdir}/genome_paths.txt \
            --rl ${outdir}/genome_paths.txt \
            -o ${outdir}/fastani_out.txt \
            -t ${SLURM_CPUS_PER_TASK} --fragLen 250 --minFraction .75

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully!"
else
  echo "Error occurred!"
  exit $RC
fi
