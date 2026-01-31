#!/bin/bash

#===============================================================================
# Name         : s06_trimmomatic.sh
# Description  : Trims adapter sequences and low-quality bases from illumina
#              : reads, producing a "cleaned" set of paried reads
# Usage        : sbatch s06_trimmomatic.sh
#===============================================================================
#SBATCH --job-name=trimmomatic
#SBATCH --array=1-15
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --output=slurm_out/trim/x_trim_%a_%A.out
#SBATCH --error=slurm_out/trim/y_trim_%a_%A.err

eval $( spack load --sh trimmomatic )

basedir="$PWD"
indir="/lts/gdlab/raw_data/seq5/"
outdir="${basedir}/d01_trimmedreads"

mkdir -p ${outdir}

# need to declare memory explicitely
export JAVA_ARGS="-Xmx1000M"

# choose which adapters to use:
# NexteraPE-PE.fa,
# TruSeq2-PE.fa, TruSeq2-SE.fa,
# TruSeq3-PE-2.fa, TruSeq3-PE.fa, TruSeq3-SE.fa
adapt="/ref/gdlab/data/trimmomatic_adapters/0.39/NexteraPE-PE.fa"

# samplelist.txt contains a list of all file names
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/r2_all_seq_list.txt`
gtac_name=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/r2_gtac_seq_list.txt`

# R1 = fwd and R2 = rev
# P = paired, UP = unpaired
set -x
time trimmomatic \
    PE \
    -phred33 \
    -trimlog \
    ${outdir}/Paired_${sample}_trimlog.txt \
    ${indir}/${gtac_name}_L003_R1_001.fastq.gz \
    ${indir}/${gtac_name}_L003_R2_001.fastq.gz \
    ${outdir}/${sample}_FW_clean.fastq \
    ${outdir}/${sample}_FW_clean_UP.fastq.gz \
    ${outdir}/${sample}_RV_clean.fastq \
    ${outdir}/${sample}_RV_clean_UP.fastq.gz \
    ILLUMINACLIP:${adapt}:2:30:10:1:true \
    SLIDINGWINDOW:4:20 \
    LEADING:10 \
    TRAILING:10 \
    MINLEN:60
RC=$?
set +x

# combine unpaired reads into a single fastq and delete redundant files
cat ${outdir}/${sample}_FW_clean_UP.fastq.gz \
    ${outdir}/${sample}_RV_clean_UP.fastq.gz \
    > ${outdir}/${sample}_UP_clean.fastq.gz
rm ${outdir}/${sample}_FW_clean_UP.fastq.gz
rm ${outdir}/${sample}_RV_clean_UP.fastq.gz 

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
