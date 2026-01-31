#!/bin/bash
#===============================================================================
# File Name    : s03_CheckV.sh
# Description  : Assessing the quality of metagenome-assembled viral genomes
# Usage        : sbatch s03_CheckV.sh
# Author       : Kailun Zhang, modified by Sri Paruthiyil paruthiyil@wustl.edu
# Version      : 1.0
# Modified     : Nov 02 2021
# Created      : May 14 2022
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=checkv
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=slurm/checkv/checkv_%A.out
#SBATCH --error=slurm/checkv/checkv_%A.err

#eval $( spack load --sh miniconda3 )
eval $( spack load --sh /buvzc6u ) # miniconda3@4.10.3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /ref/gdlab/software/envs/CheckV
export CHECKVDB=/ref/gdlab/software/envs/CheckV/checkv-db-v1.5

#Requires user to take all fna files from CT3 and concatenate them to get the input file
set -x
checkv end_to_end /scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/combinedviruses/250114_combinedviruses.fna /scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/checkv/rawcombinedviruses -t 16

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
