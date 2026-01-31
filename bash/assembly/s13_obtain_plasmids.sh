#!/bin/bash
#===============================================================================
#
# File Name    : s13_get_plasmid_fastas.sh
# Description  : This script will get plasmid fastas from flye assemblies
# Usage        : sbatch s13_get_plasmid_fastas.sh
# Author       : Emily Benedict & Sri Paruthiyil
# Version      : 1.0
# Created On   : 2023-09-16
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=plasmidfinder
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=slurm_out/plasmidfastas/z_plasmidfinder_%a.out
#SBATCH --error=slurm_out/plasmidfastas/z_plasmidfinder_%a.out
#SBATCH --mail-user=ebenedict@wustl.edu
#SBATCH --mail-type=END

source_dir="/scratch/gdlab/ebenedict/preterm_plasmidome/round3/d03_flye"
target_dir="/scratch/gdlab/ebenedict/preterm_plasmidome/d06_plasmids"

cat "/scratch/gdlab/ebenedict/preterm_plasmidome/round3/test_p4_r3_barcodes.txt" | while read x; do
    # Create the directory path for the current iteration
    current_dir="${source_dir}/${x}"

    # Check if the directory exists
    if [ -d "$current_dir" ]; then
        # Define the source and target file paths
        assembly_fasta="${current_dir}/assembly.fasta"
        assembly_info="${current_dir}/assembly_info.txt"
        
        # Check if the source files exist
        if [ -f "$assembly_fasta" ] && [ -f "$assembly_info" ]; then
            # Define the target file names with the value of x
            target_assembly_fasta="${target_dir}/${x}_assembly.fasta"
            target_assembly_info="${target_dir}/${x}_assembly_info.txt"
            
            # Move and rename the files to the target directory
            rsync -avr "$assembly_fasta" "$target_assembly_fasta"
            rsync -avr "$assembly_info" "$target_assembly_info"
            
            echo "Moved files from $current_dir to $target_dir with names ${x}_assembly.fasta and ${x}_assembly_info.txt"
        else
            echo "Assembly files not found in $current_dir"
        fi
    else
        echo "Directory $current_dir does not exist"
    fi
done

for txt_file in "$target_dir"/*_info.txt; do
    # Extract the base name without extension (contig name)
    contig_name=$(basename "$txt_file" "_info.txt")
    
    # Check if the 4th column contains 'Y' and the 2nd column is less than 1 million
    if awk -F'\t' '$4=="Y" && $2<1000000' "$txt_file"; then
        # Extract the corresponding contigs from the fasta file
        awk -F'\t' '$4=="Y" && $2<1000000 {print $1}' "$txt_file" | \
        while read -r contig; do
            sed -n "/^>$contig$/,/^>/ {/^>/!p}" "$target_dir/${contig_name}.fasta"
        done > "${target_dir}/${contig_name}_plasmids.fasta"
    fi
done
