#!/bin/bash
#===============================================================================
# Name         : s11_extract_flye_info.sh
# Usage        : sbatch s11_extract_flye_info.sh
# Created On   : 2023-09-10
# Author       : Jian Ryou, j.ryou@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=0-6:00:00 # days-hh:mm:ss
#SBATCH --job-name=parse
#SBATCH --mail-type=END
#SBATCH --mail-user=j.ryou@wustl.edu
#SBATCH --array=1-45
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/grep/z_grep-3_%a.out
#SBATCH --error=slurm_out/grep/z_grep-3_%a.out

#set up directories
basedir="$PWD"
flye_in="/scratch/gdlab/ebenedict/preterm_plasmidome/round3/d03_flye/round3_barcodes"
outdir="${basedir}/d06_QC/flye-out"
mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/polish_barcodes.txt`

header="Sample\tTotal length\tFragments\tFragments_N50\tLargest_Fragment\tScaffolds\tMean_Coverage\tLongest_contig (bp)\tLength (bp)\tCoverage\tCircular\tPlasmids\n"

# Go into each flye barcode assembly directory and find assembly_info.txt and flye.log
## First, rename and move assembly fastas and assembly infos from flye 
### check if the directory with numbered with array number exists
if [ -d "${flye_in}/${sample}" ]; then
	# Define the source and target file paths
	assembly_log="${flye_in}/${sample}/flye.log"
	assembly_info="${flye_in}/${sample}/assembly_info.txt"
	assembly_fasta="${flye_in}/${sample}/assembly.fasta"
	# Check if the source file exists
	if [ -f "$assembly_log" ] && [ -f "$assembly_info" ]; then
		# Define the target file names with the barcode name
		target_assembly_log="${outdir}/${sample}_flye.log"
		target_assembly_info="${outdir}/${sample}_assembly_info.txt"
		tmp_info="${outdir}/tmp_${sample}_assembly_info.txt"
		target_assembly_fasta="${outdir}/${sample}_assembly.fasta"
		combined_info="${outdir}/${sample}_comb_info.txt"
		# Copy and rename the files to the target directory
		cp "$assembly_info" "$target_assembly_info"
		cp "$assembly_log" "$target_assembly_log"
		
		# Extract info from assembly log
		total_lengh=`awk '$1=$1' ${target_assembly_log} | grep -n "Total length:" | cut -d ":" -f3`
		fragments=`awk '$1=$1' ${target_assembly_log} | grep -n "Fragments:" | cut -d ":" -f3`
		fragments_n50=`awk '$1=$1' ${target_assembly_log} | grep -n "Fragments N50:" | cut -d ":" -f3`
		largest_frg=`awk '$1=$1' ${target_assembly_log} | grep -n "Largest frg:" | cut -d ":" -f3`
		scaffolds=`awk '$1=$1' ${target_assembly_log}| grep -n "Scaffolds:" | cut -d ":" -f3`
		mean_cov=`awk '$1=$1' ${target_assembly_log}| grep -n "Mean coverage:" | cut -d ":" -f3`


		# Write combined info file for each sample

		echo -n -e ${header} >> ${combined_info}
		echo -n -e "${sample}\t${total_lengh}\t${fragments}\t${fragments_n50}\t${largest_frg}\t${scaffolds}\t${mean_cov}" >> ${combined_info}
		
		# Extract info from assembly info txt file
		#longest_contig_info=`awk 'NR==2 { print $1, $2, $3, $4 }' ${target_assembly_info}`
		longest_contig_info=`head -2 ${target_assembly_info} | tail -1 | cut -f1-4`
		echo -n -e "\t${longest_contig_info}\t" >> ${combined_info}

		echo "Printed longest contig info into combined_info file"

		# Find putative plasmids
		cut -f 1-4 ${target_assembly_info} >> ${tmp_info} 
		if awk -F'\t' '$4=="Y" && $2<1000000' "${tmp_info}"; then
			# Store plasmid info
			plasmid_info=`awk -F'\t' '$4=="Y" && $2<1000000' "${tmp_info}"`

			tmp_plasmid_info="${outdir}/${sample}_tmp_plasmid.txt"
			echo -n -e "${plasmid_info}" >> ${tmp_plasmid_info}

			while IFS= read -ra line; do
				for i in ${line}; do
					echo -n -e "$i\t" >> ${combined_info}
				done
			done < 	"${tmp_plasmid_info}"
			echo "Completed writing of plasmids"

			#echo -n -e "${plasmid_info}" >> ${combined_info}

			# Extract the corresponding contigs from the fasta file
			# awk -F'\t' '$4=="Y" && $2<1000000 {print $1}' "${target_dir}/${SLURM_ARRAY_TASK_ID}_assembly_info.txt" | \
		    # while read -r contig; do
		    #     sed -n "/^>$contig$/{p;n;p}" "$target_dir/${SLURM_ARRAY_TASK_ID}_consensus.fasta"
		    # done > "${target_dir}/${SLURM_ARRAY_TASK_ID}_plasmid_contigs.fasta"
		fi

		echo "Copied files from $flye_in to $outdir with  names $target_assembly_info and $target_assembly_log"
	else 
		echo "Assembly log and info files not found for $sample"
	fi
else
	echo "Directory ${flye_in}/${sample} does not exist!"
fi



# Rename those two files with the barcode
# Copy them to new output directory
# In new directory, loop through each file and extract the following info: 
## "Assembly statistics:" in flye.log and the #seq_name (col1), length (col2), cov (col3),
## and circ. (col4). Extract this infor for first two rows (the two longest contigs) and 
## any circular plasmids 
