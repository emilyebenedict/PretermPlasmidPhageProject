#!/bin/bash
#===============================================================================
# Name         : s17_mash_summary.sh
# Usage        : sbatch s17_mash_summary.sh
# Created On   : 2023-11-02
# Author       : Jian Ryou, j.ryou@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=0-6:00:00 # days-hh:mm:ss
#SBATCH --job-name=mash_cat
#SBATCH --mail-type=END
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --array=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/mash/z_combine_mash_%a.out
#SBATCH --error=slurm_out/mash/z_combine_mash_%a.out

#set up directories
basedir="$PWD"
mash_in="${basedir}/d05_mash/refseq_round3"
outdir="${basedir}/d05_mash"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/p4_round3_barcodes.txt`

mash_comb_out="${outdir}/p4_round3_mash_comb.txt"
header="Sample\tMash Identity\tShared Hashes\tP-value\tQuery-ID\tContaminant Hashes\n"
echo -n -e ${header} >> ${mash_comb_out}

for file in ${mash_in}/*.tab; do
	# Extract the base name without extension (contig name)
    sample_name=$(basename "$file" ".tab")
    # Extract only fields with a shared hash greater than or equal to 700
    while IFS= read -r line; do
    	echo -n -e "${sample_name}\t${line}\n" >> ${mash_comb_out}
    done < ${file}
done

echo "Done combining all mash output files"
