#!/bin/bash
#===============================================================================
#
# File Name    : s12c_coverm.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12c_coverm.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : 2022-06-01
# Last Modified: 2022-06-01
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=coverm
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1%50
#SBATCH --output=slurm_out/coverm_plasmid_count/x_coverm_%A_%a.out
#SBATCH --error=slurm_out/coverm_plasmid_count/y_coverm_%A_%a.out

#total is 3009 for phage, 894 for plasmid
eval $( spack load --sh miniconda3@4.10.3 )

CONDA_BASE=$(conda info --base)

source $CONDA_BASE/etc/profile.d/conda.sh

conda activate /home/gmark/conda_env/coverM # version 0.7.0


basedir="/scratch/gdlab/p.sri/250113_P4Pipeline/"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250113_P4Pipeline/250610_bowtieplasmidsandphages/250613_plasmidmapping_sampleIDpermute.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250113_P4Pipeline/250610_bowtieplasmidsandphages/250613_plasmidmapping_plasmidpermute.txt`
indir="/scratch/gdlab/p.sri/250113_P4Pipeline/250610_bowtieplasmidsandphages/d12_bowtie2/${sample}_${ref}"
outdir="/scratch/gdlab/p.sri/250113_P4Pipeline/250610_bowtieplasmidsandphages/d12_bowtie2/d12_coverm_count"

#make output directory
mkdir -p ${outdir}

time coverm genome --bam-files ${indir}/${sample}_sorted.bam --genome-fasta-files /scratch/gdlab/p.sri/250113_P4Pipeline/250610_bowtieplasmidsandphages/splitplasmids/${ref}_plasmids.fa  -m relative_abundance mean trimmed_mean covered_bases count --min-read-percent-identity 95 -o ${outdir}/${sample}_${ref}_coverm_output.tsv 
#time coverm genome --bam-files ${indir}/${sample}_sorted.bam --genome-fasta-files /scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/combinedviruses/splitviruses/${ref}.fa  -m relative_abundance mean trimmed_mean covered_bases count --min-read-percent-identity 95 -o ${outdir}/${sample}_${ref}_coverm_output.tsv 

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi

