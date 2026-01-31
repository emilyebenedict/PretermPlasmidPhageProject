#!/bin/bash          
#===============================================================================      
# Name         : s07.1_CenoteTaker3.sh
# Description  : Discover and annoate virus sequences/genomes    
# Usage        : sbatch s07.1_CenoteTaker3.sh
# Author       : Kailun Zhang, modified by Sri Paruthiyil paruthiyil@wustl.edu       
# Version      : 1.0 
# Created On   : 2024-03-07
#===============================================================================      
#Submission script for HTCF
#SBATCH --job-name=CT3 
#SBATCH --array=2-201
#SBATCH --mail-type=END
#SBATCH --mail-user=paruthiyil@wustl.edu
#SBATCH --mem=48G      
#SBATCH --cpus-per-task=12
#SBATCH --output=/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/slurm/cenote/x_ct3_%A_%a.out  
#SBATCH --error=/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/slurm/cenote/y_ct3_%A_%a.err

eval $( spack load --sh /buvzc6u ) # miniconda3@4.10.3
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
source activate /ref/gdlab/software/envs/Mamba/Mambaforge-Linux-x86_64/envs/ct3_env

basedir="/scratch/gdlab/p.sri/250113_P4Pipeline"
indir="/lts/gdlab/users/current/ebenedict/p4/allrounds_polished_assemblies"
#outdir="${basedir}/d10a_megahit_ct3"
#mkdir -p ${outdir}

# Read in the slurm array task                               
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/mapping/250113_emilyfinalassemblies.txt`
runtitle=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/mapping/250113_emilyfinalassemblies_cenotenames.txt`

set -x

cenotetaker3 -c ${indir}/${sample}.fastq -r ${runtitle} -p T --lin_minimum_hallmark_genes 0 --circ_minimum_hallmark_genes 0 -t ${SLURM_CPUS_PER_TASK}

RC=$?

set +x

#output if job was successful
if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi
