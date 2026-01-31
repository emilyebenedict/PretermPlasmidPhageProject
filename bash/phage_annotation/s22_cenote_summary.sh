#!/bin/bash
#===============================================================================
#
# File Name    : s01_Cenote-Taker2_sum.sh
# Description  : Summarize busco data into a single file
# Usage        : s01_Cenote-Taker2_sum.sh
# Author       : Kailun Zhang, modified by Sri Paruthiyil paruthiyil@wustl.edu
# Version      : 1.0
# Created On   : Nov 16 2021
# Modified On  : Jan 13 2023
#===============================================================================
#
#Submission script for HTCF
#SBATCH --job-name=sum
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --output=/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/slurm/cenote/x_ct3sum_%A_%a.out  
#SBATCH --error=/scratch/gdlab/p.sri/250113_P4Pipeline/250113_CT3/slurm/cenote/y_ct3sum_%A_%a.err
basedir="/scratch/gdlab/p.sri/250113_P4Pipeline"
indir="${basedir}/250113_CT3"
#CT2dir="${indir}"
#samplelist="/scratch/gdlab/p.sri/230807_UPEC/ecoli_mapping_key_justnumbers.txt"
#outfile="${basedir}/240429_sample_list_new.txt"
samplelist="/scratch/gdlab/p.sri/250113_P4Pipeline/mapping/250113_emilyfinalassemblies_cenotenames.txt"
outfile="${basedir}/250113_CT3/250114_cenote_summary.tsv"
# create output file
touch ${outfile}
header="ORIGINAL_NAME\tCENOTE_NAME\tORGANISM_NAME\tEND_FEATURE\tLENGTH\tORF_CALLER\tNUM_HALLMARKS\tHALLMARK_NAMES\tBLASTN_INFO"
echo -e ${header} >> ${outfile}
# loop through samplelist and append stats to outfile
while read sample; do
        # load data
        CT_Data=`grep -P "${sample}" ${indir}/${sample}/${sample}_virus_summary.tsv`
        echo -e "${CT_Data}" >> ${outfile}
done < ${samplelist}
