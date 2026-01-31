#!/bin/bash
#===============================================================================
# Name         : s12_cat_flye_info_files.sh
# Usage        : sbatch s12_cat_flye_info_files.sh (run after s_extract_flye_info.sh)
# Created On   : 2023-11-02
# Author       : Jian Ryou, j.ryou@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=0-6:00:00 # days-hh:mm:ss
#SBATCH --job-name=cat
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --array=1
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/grep/z_combine_flye_%a.out
#SBATCH --error=slurm_out/grep/z_combine_flye_%a.out

#set up directories
basedir="$PWD"
flye_in="${basedir}/d06_QC/flye-out"
outdir="${flye_in}"

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/polish_barcodes.txt`

flye_comb_out="${outdir}/240201_round3_flye_comb_out.txt"
header="Sample\tTotal length\tFragments\tFragments_N50\tLargest_Fragment\tScaffolds\tMean_Coverage\tLongest_contig (bp)\tLength (bp)\tCoverage\tCircular\tPlasmids\n"
echo -n -e ${header} >> ${flye_comb_out}

for file in ${flye_in}/*_comb_info.txt; do
	tmp=`head -2 $file | tail -1`
	echo -n -e "${tmp}\n" >> ${flye_comb_out}
done

echo "Done combining all flye output files"
