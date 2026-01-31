# File name:      plasmid_tree_generate.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-08-15
# Description:    This script generates a phylogeny from FastANI output of P4 putative plasmids

# loading packages
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(phytools)



# setwd
setwd('~/Box Sync/EEB/Dantas/pppp/')

# read in fastani data
fastani = read.delim('./data/p4_plasmid_fastani_out.txt',header = F, stringsAsFactors = F)
colnames(fastani) <- c("source", "target", "ANI", "sourcelen", "targetlen")
fastani = fastani[grep('contig',fastani$source),]
fastani = fastani[grep('contig',fastani$target),]
fastani$source = gsub('/scratch/gdlab/ebenedict/preterm_plasmidome/250816_separate_plasmid_contigs/','', fastani$source)
fastani$source = gsub('plasmids.fastq_','', fastani$source)
fastani$source = gsub('.fa','', fastani$source)
fastani$target = gsub('/scratch/gdlab/ebenedict/preterm_plasmidome/250816_separate_plasmid_contigs/','', fastani$target)
fastani$target = gsub('plasmids.fastq_','', fastani$target)
fastani$target = gsub('.fa','', fastani$target)
fastani$source_sample = gsub('_contig.*','',fastani$source)
fastani$target_sample = gsub('_contig.*','',fastani$target)
fastani = fastani[fastani$target != '',]


# read in list of HQ genomes to keep in trees
mash_data = data.frame(read_excel('./data/PretermPlasmidPhageProject_polypolish_circular_mash_checkm_HQgenomes.xlsx'),  stringsAsFactors = F)
fastani = fastani[fastani$source_sample %in% mash_data$Sample,]
fastani = fastani[fastani$target_sample %in% mash_data$Sample,]

heatmap_matrix <- dcast(fastani, source ~ target, value.var = "ANI")
heatmap_matrix[is.na(heatmap_matrix)] = 50
rownames(heatmap_matrix) = heatmap_matrix$source
heatmap_matrix = heatmap_matrix[,-1]

distance_matrix = dist(heatmap_matrix)
h_clustered = hclust(distance_matrix)
plasmid_tree = as.phylo(h_clustered)

write.tree(plasmid_tree, './data/p4_plasmid_tree.treefile')
