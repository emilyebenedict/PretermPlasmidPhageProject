#!/bin/bash
#===============================================================================
#
# File Name    : s29_metaphlan.sh
# Description  : This script will concatenate metaphlan output
# Usage        : sbatch s29_metaphlan.sh
# Author       : Luke Diorio-Toth, modified by Jian Ryou j.ryou@wustl.edu
# Version      : 1.0
# Created On   : Wed Jan  6 11:49:47 CST 2021
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-0:00:00 # days-hh:mm:ss
#SBATCH --job-name=merge
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=slurm_out/metaphlan_merge/merge_%a_%A.out
#SBATCH --error=slurm_out/metaphlan_merge/merge_%a_%A.out

eval $( spack load --sh py-metaphlan@4.0.2 )

basedir="$PWD"
indir="${basedir}/metaphlan_all"
outdir="${basedir}/merged_outputs"
mkdir -p ${outdir}


set -x # Start debug mode (will send commands to outfile)

# Merge metaphlan profiles
merge_metaphlan_tables.py  ${indir}/*_profile.txt > ${outdir}/240923_P4_GM_dev_merged_metaphlan_n977.txt

RC=$? # Save error code for command
set +x # End debug mode (will stop sending commands to outfile)

# Output if job was successful
if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occured in ${sample}!"
    exit $RC
fi
