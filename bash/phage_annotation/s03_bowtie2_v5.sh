#!/bin/bash
#===============================================================================
#
# File Name    : s12b_bowtie2.sh
# Description  : This script will align reads to a ref using bowtie2
# Usage        : sbatch s12b_bowtie2.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: May 31 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtie2
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=2-36
#SBATCH --output=slurm/bowtie2/x_bowtie2_%A_%a.out
#SBATCH --error=slurm/bowtie2/y_bowtie2_%A_%a.out

eval $( spack load --sh bowtie2@2.4.2 )
#load samtools/1.14
eval $( spack load --sh /6p5wlkk )

basedir="/scratch/gdlab/p.sri/250414_HVP_induction"
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250414_HVP_induction/bwamapping_permutereads.txt`
ref=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250414_HVP_induction/bwamapping_permutegenome.txt`
#indir="/lts/gdlab/users/current/kailun/Phageome/NEC_Olm_Reads/BBtools/d03_bbtools_single"
#index is where the output from bowtie2-build is stored
index="/scratch/gdlab/p.sri/250414_HVP_induction/bowtieBuild/${ref}/${ref}"
outdir="/scratch/gdlab/p.sri/250414_HVP_induction/d12_bowtie2/${sample}_${ref}"

#make output directory
mkdir -p ${outdir}

time bowtie2 -p ${SLURM_CPUS_PER_TASK} \
    -x ${index} \
    -S ${outdir}/${sample}_mappedreads.sam \
    -1 /scratch/gdlab/p.sri/250414_HVP_induction/d01a_trimmomatic/${sample}_FW_trim.fastq.gz \
    -2 /scratch/gdlab/p.sri/250414_HVP_induction/d01a_trimmomatic/${sample}_RV_trim.fastq.gz

# convert sam to bam, then generate sorted bam file and index it.
samtools view -bS ${outdir}/${sample}_mappedreads.sam > ${outdir}/${sample}_raw.bam
samtools sort ${outdir}/${sample}_raw.bam > ${outdir}/${sample}_sorted.bam
samtools index ${outdir}/${sample}_sorted.bam
# generate stats
samtools depth -aa ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_depth.txt
samtools stats ${outdir}/${sample}_sorted.bam > ${outdir}/${sample}_stats.txt

rm ${outdir}/${sample}_mappedreads.sam
rm ${outdir}/${sample}_raw.bam

RC=$?
set +x

if [ "$RC" -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi

