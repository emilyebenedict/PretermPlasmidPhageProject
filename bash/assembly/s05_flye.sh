#!/bin/bash
#===============================================================================
#
# File Name    : s05_flye.sh
# Description  : This script will assemble ONT long reads with Flye.
# Usage        : sbatch s06_flye.sh
# Author       : Erin Newcomer
# Modified     : Emily Benedict, ebenedict@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=flye
#SBATCH --time=3-00:00:00 # days-hh:mm:ss
#SBATCH --array=1-96%12
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=slurm_out/flye/flye-redo_%a.out
#SBATCH --error=slurm_out/flye/flye-redo_%a.err
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --mail-type=END

eval $( spack load --sh py-flye )

basedir="$PWD"
indir="${basedir}/d02_filtlong"
outdir="${basedir}/d03_flye"

mkdir -p ${outdir}

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/round3_barcodes.txt`

#read in the slurm array task
echo ${sample}
set -x

time flye \
	--nano-hq "${indir}/SQK-NBD114-96_barcode${sample}.fastq.gz" \
	--out-dir "${outdir}/${sample}" \
        --meta \
	--threads ${SLURM_CPUS_PER_TASK} \
        --iterations 5 \


RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred in ${sample}!"
  exit $RC
fi
