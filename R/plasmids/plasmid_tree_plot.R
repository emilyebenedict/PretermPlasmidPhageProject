# File name:      plasmid_tree_plot.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-08-15
# Description:    This script annotates the p4 putative plasmid tree

# loading packages
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(NatParksPalettes)
library(RColorBrewer)
library(phytools)
library(ggtree)
library(pheatmap)
library(tidytree)
library(ggnewscale)
library(lubridate)
library(kableExtra)
library(ggpubr)
library(qgraph)


# setwd
setwd('~/Box Sync/EEB/Dantas/pppp/')

#### reading in and processing data ####

# read in tree
tree = read.tree('./data/p4_plasmid_tree.treefile') # tips: 436

# reading in mash data
mash_data = data.frame(read_excel('./data/PretermPlasmidPhageProject_polypolish_circular_mash_checkm_HQgenomes.xlsx'), stringsAsFactors = F)
mash_data$Sample = gsub('-','_', mash_data$Sample)
mash_data$Sample = gsub('Ent','ENT', mash_data$Sample)
mash_data$Sample = gsub('Se','SE',mash_data$Sample)
mash_data$Sample = gsub('Sa','SA',mash_data$Sample)
mash_data$Sample = gsub('Ko','KO',mash_data$Sample)
mash_data$Sample = gsub('Kp','KP',mash_data$Sample)
mash_data$Sample = gsub('Sw','SW',mash_data$Sample)
mash_data$genus = gsub(' .*','', mash_data$simplified_mash_call)

# reading in patient data
pts = data.frame(read_excel('./data/250220_P4_GM_Dev_Abx_groupings_updated.xlsx'), stringsAsFactors = F)
colnames(pts)[2] = 'Infant'
class(pts$Infant) = 'character'
pts$Site = gsub('St. Louis','MO',pts$Site)
pts$Site = gsub('Oklahoma','OK',pts$Site)
pts$Site = gsub('Louisville','KY',pts$Site)
pts = pts[pts$Infant %in% mash_data$Infant,]
pts$abx_group = gsub('A','a',pts$abx_group)
pts$abx_group = gsub('after','After',pts$abx_group)
class(mash_data$Infant) = 'character'

mash_pt_data = full_join(mash_data, pts, relationship = 'many-to-many')

# reading in cluster labels
ani_cluster_labels = read.csv('./data/251130_p4_plasmid_ANI_cluster_labels.csv', stringsAsFactors = F)

# reading in colors
colors = data.frame(read_excel('./notes/p4_colors.xlsx'), stringsAsFactors = F)

# ordering dataframes for plotting
ordered_df = data.frame(tree$tip.label, row.names = tree$tip.label)
colnames(ordered_df) = 'Tips'
ordered_df$Sample = gsub('_contig.*','', ordered_df$Tips)
ordered_df = ordered_df[ordered_df$Sample %in% mash_data$Sample,]
tree = drop.tip(tree, setdiff(tree$tip.label, ordered_df$Tips))
ani_cluster_labels = ani_cluster_labels[ani_cluster_labels$Contig %in% tree$tip.label,]
for(i in 1:nrow(ordered_df)){
  ordered_df$Mash[i] = mash_pt_data$simplified_mash_call[mash_pt_data$Sample == ordered_df$Sample[i]]
  ordered_df$NICU[i] = mash_pt_data$Site[mash_pt_data$Sample == ordered_df$Sample[i]]
  ordered_df$abx_group[i] = mash_pt_data$abx_group[mash_pt_data$Sample == ordered_df$Sample[i]]
  ordered_df$genus[i] = mash_pt_data$genus[mash_pt_data$Sample == ordered_df$Sample[i]]
  
  if(ordered_df$Tips[i] %in% ani_cluster_labels$Contig){
    ordered_df$ani_cluster[i] = ani_cluster_labels$cluster[ani_cluster_labels$Contig == ordered_df$Tips[i]]
  }else{
    ordered_df$ani_cluster[i] = 'none'
  }
  
} 
ordered_df$abx_group = gsub(' After',' after',ordered_df$abx_group)

# separating the dataframes
mash_df = ordered_df[,c('Tips','Mash')]
nicu_df = ordered_df[,'NICU', drop=F]
abx_df = ordered_df[,'abx_group',drop=F]
ani_cluster_df = ordered_df[,'ani_cluster',drop=F]
genus_df=ordered_df[,'genus',drop=F]

#### plotting the trees!! ####
mash_colors = colors$Hex[colors$Category == 'Simplified Mash']
names(mash_colors) = colors$Name[colors$Category == 'Simplified Mash']

treeplot = ggtree(tree) %<+% mash_df
circle_treeplot = ggtree(tree, layout = 'fan') %<+% mash_df

# plotting genus ring + abx group + plasmid cluster
genus_df$genus[genus_df$genus == 'Shigella'] = 'Escherichia'
genus_colors = colors$Hex[colors$Category == 'genus']
names(genus_colors) = colors$Name[colors$Category == 'genus']

genus_tree = gheatmap(circle_treeplot, genus_df, colnames_angle = 85, colnames_offset_y = -1,
         offset = 0, width = 0.08, font.size = 5)+
  geom_treescale(x = 30, y = 1)+
  scale_fill_manual(values = genus_colors, breaks = names(genus_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Host species genus')

genus_tree2 = genus_tree + new_scale_fill()

# adding abx group ring
abx_colors = colors$Hex[colors$Category == 'Abx cohort']
names(abx_colors) = colors$Name[colors$Category == 'Abx cohort']


genus_abxgroup_tree = gheatmap(genus_tree2, abx_df, colnames_angle = 85, colnames_offset_y = -1,
                               offset = 30, width = 0.08, font.size = 5)+
  geom_treescale(x = 30, y = 1)+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Infant antibiotic exposure')

genus_abxgroup_tree2 = genus_abxgroup_tree + new_scale_fill()

# adding ANI-based cluster
ani_cluster_colors = colors$Hex[colors$Category == 'plasmid_cluster']
names(ani_cluster_colors) = colors$Name[colors$Category == 'plasmid_cluster']

genus_abxgroup_cluster_tree = gheatmap(genus_abxgroup_tree2, ani_cluster_df, colnames_angle = 85, colnames_offset_y = -1,
                               offset = 60, width = 0.08, font.size = 5)+
  geom_treescale(x = 30, y = 1)+
  scale_fill_manual(values = ani_cluster_colors, breaks = names(ani_cluster_colors))+ggtree::vexpand(0.1,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.position = 'right')+
  labs(fill = 'Plasmid cluster')

#pdf(file = paste0('./figures/',Sys.Date(),'_p4_anicluster_abxgroup_genus_plasmid_tree.pdf'), width = 16, height = 16)
genus_abxgroup_cluster_tree
#dev.off()
