#!/bin/bash
#===============================================================================
# Name         : s08_multiqc.sh
# Description  : Consolidates output from s01_fastqc.sh
# Usage        : sbatch s08_multiqc.sh
# Author       : Erin Newcomer, modified by Emily Benedict, ebenedict@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=multiqc
#SBATCH --mem=8G
#SBATCH --output=slurm_out/multiqc/z_multiqc_%A.out
#SBATCH --error=slurm_out/multiqc/z_multiqc_%A.out
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu

eval $( spack load --sh py-multiqc )

basedir="$PWD"
indir="${basedir}/d03_read_qc_round2/subsampled"
rawin="${indir}/raw_reads/fastqc"
cleanin="${indir}/clean_reads/fastqc"
rawout="${indir}/raw_reads/multiqc"
cleanout="${indir}/clean_reads/multiqc"

mkdir -p ${outdir}

set -x
time sh -c \
"multiqc ${rawin} -o ${rawout} -n raw_multiqc; \
multiqc ${cleanin} -o ${cleanout} -n clean_multiqc"
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occurred!"
  exit $RC
fi
