#!/bin/bash
#===============================================================================
# File Name    : s06_vcontact2_binned.sh
# Description  : This script will run amrfinder in parallel
# Usage        : sbatch s06_vcontact2_binned.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.3
# Created On   : Dec 14 2021
# Last Modified: Feb 20 2023
#===============================================================================
#Submission script for HTCF
#SBATCH --job-name=vcontact2
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --output=slurm/vcontact/x_vcontact2_%A.out
#SBATCH --error=slurm/vcontact/y_vcontact2_%A.out


eval $( spack load --sh miniconda3@4.10.3 )
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

# activate env
conda activate /home/kailun/conda/vContact2 # version 0.9.19

#module load openjdk
eval $( spack load --sh openjdk@11.0.15_10 )

basedir="/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/"

set -x
time vcontact2 --raw-proteins /scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/checkv/rawcombinedviruses/tmp/proteins.faa --rel-mode 'Diamond' --proteins-fp /scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/checkv/rawcombinedviruses/tmp/250114_CT3_G2G.csv --db 'ProkaryoticViralRefSeq207-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/kailun/conda/MAVERICLab-vcontact2-34ae9c466982/bin/cluster_one-1.0.jar --output-dir ${basedir}/vConTACT_250114


RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured in vContact2"
  exit $RC
fi
