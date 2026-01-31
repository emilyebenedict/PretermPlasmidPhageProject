#!/bin/bash
#===============================================================================
# me           : s08_fastqc.sh
# Description  : This script will run fastqc in parallel
# Usage        : sbatch s01_fastQC.sh
# Author       : Erin Newcomer, modified by Emily Benedict, ebenedict@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=fastqc
#SBATCH --array=1-15
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm_out/fastqc/z_fastqc_%A_%a.out
#SBATCH --error=slurm_out/fastqc/z_fastqc_%A_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $(spack load --sh fastqc@0.11.9)

basedir="$PWD"
rawin="/lts/gdlab/raw_data/seq5/"
cleanin="${basedir}/d01_trimmedreads"
rawout="${basedir}/d02_cleanreads/raw_reads"
cleanout="${basedir}/d02_cleanreads/clean_reads"
mkdir -p ${rawout}
mkdir -p ${cleanout}

export JAVA_ARGS="-Xmx8000M"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/r2_all_seq_list.txt`
gtac_name=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/r2_gtac_seq_list.txt`
set -x

time sh -c \
"fastqc ${rawin}/${gtac_name}_L003_R1_001.fastq.gz ${rawin}/${gtac_name}_L003_R2_001.fastq.gz -o ${rawout} -t ${SLURM_CPUS_PER_TASK}; \
fastqc ${cleanin}/${sample}_FW_clean.fastq ${cleanin}/${sample}_RV_clean.fastq ${cleanin}/${sample}_UP_clean.fastq -o ${cleanout} -t ${SLURM_CPUS_PER_TASK}"
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
