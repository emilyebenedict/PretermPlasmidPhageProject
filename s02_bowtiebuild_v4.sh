#!/bin/bash
#===============================================================================
#
# File Name    : s12a_bowtieBuild.sh
# Description  : This script will index a ref for Burrows-Wheeler alignment
# Usage        : sbatch s12a_bowtieBuild.sh
# Author       : Kailun Zhang, kailun@wustl.edu
# Version      : 1.0
# Created On   : Jul 09 2019 by Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Last Modified: May 31 2020
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=bowtieBuild
#SBATCH --array=1-3
#SBATCH --mem=4G
#SBATCH --output=slurm/bowtie/x_bowtieBuild_%A_%a.out
#SBATCH --error=slurm/bowtie/y_bowtieBuild_%A_%a.out

eval $( spack load --sh bowtie2 )


basedir="/scratch/gdlab/p.sri/250414_HVP_induction/"
indir="/scratch/gdlab/p.sri/250414_HVP_induction/genomes"
outdir="/scratch/gdlab/p.sri/250414_HVP_induction/bowtieBuild"
mkdir -p ${outdir}



ref=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/gdlab/p.sri/250414_HVP_induction/bwamapping.txt`
#ref=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/240821_coassembly_list_160.txt`
#sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/scaffold_list_new.txt`
mkdir -p ${outdir}/${ref}

set -x

time 
#grep "${ref}" /scratch/gdlab/p.sri/240429_phageome_P4/CT3/finalfna/240430_all_virus_sequence.fna | sed 's/>//' > ${outdir}/${ref}/241014_phage_contigs.txt
#cd ${outdir}/${ref}
#awk 'BEGIN{while((getline<"241014_phage_contigs.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f' /scratch/gdlab/p.sri/240429_phageome_P4/CT3/finalfna/240430_all_virus_sequence.fna >$


#grep -A 1 ">${ref}" ${indir}/241001_Olm_SRR_CT3_pass4.fasta > ${basedir}/d12_bowtie/${ref}_contigs.fasta
#sed -i '/--/d' ${ref}_contigs.fasta

bowtie2-build /scratch/gdlab/p.sri/250414_HVP_induction/genomes/${ref}.fa ${outdir}/${ref}/${ref}

RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error Occured!"
  exit $RC
fi

