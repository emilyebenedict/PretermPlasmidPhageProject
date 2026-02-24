# File name:      250815_p4_plasmid_tree.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-08-15
# Description:    This script annotates the p4 plasmid tree

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
setwd('~/Library/CloudStorage/Box-Box/pppp/')

# read in tree
tree = read.tree('~/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/VIPtree/260128/all.bionj.asc.newick') # tips: 436

# reading in mash data
mash_data = data.frame(read_excel('/Users/gdlab/Library/CloudStorage/Box-Box/pppp/data/PretermPlasmidPhageProject_polypolish_circular_mash_checkm_HQgenomes.xlsx'), stringsAsFactors = F)
mash_data$Sample = gsub('-','_', mash_data$Sample)
mash_data$Sample = gsub('Ent','ENT', mash_data$Sample)
mash_data$Sample = gsub('Se','SE',mash_data$Sample)
mash_data$Sample = gsub('Sa','SA',mash_data$Sample)
mash_data$Sample = gsub('Ko','KO',mash_data$Sample)
mash_data$Sample = gsub('Kp','KP',mash_data$Sample)
mash_data$Sample = gsub('Sw','SW',mash_data$Sample)

# reading in patient data
pts = data.frame(read_excel('/Users/gdlab/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/260202_P4_GM_dev_meta_updated_groupings_v2 (1).xlsx'), stringsAsFactors = F)
#Removing all after only patients

pts <- pts %>% filter(Patient != "184.01")
pts <- pts %>% filter(Patient != "2031.01")
# pts <- pts %>% filter(Patient != "278.01")
# 
# pts <- pts %>% filter(Patient != "2038.01")
# 
# pts <- pts %>% filter(Patient != "2018.01")
# 
# pts <- pts %>% filter(Patient != "2007.01")
# 
# pts <- pts %>% filter(Patient != "1027.01")
colnames(pts)[6] = 'Infant'
class(pts$Infant) = 'character'
pts$Site = gsub('St. Louis','MO',pts$Site)
pts$Site = gsub('Oklahoma','OK',pts$Site)
pts$Site = gsub('Louisville','KY',pts$Site)
pts = pts[pts$Infant %in% mash_data$Infant,]
pts$abx_group_updated = gsub('A','a',pts$abx_group_updated)
pts$abx_group_updated = gsub('after','After',pts$abx_group_updated)
class(mash_data$Infant) = 'character'

mash_pt_data = full_join(mash_data, pts, relationship = 'many-to-many')
# mash_pt_data <- left_join(
#   mash_data,
#   pts %>% select(Infant, abx_group_updated, Site),
#   by = "Infant"
# )



# reading in plasmidfinder data
# plasmidfinder = read.table('./data/p4_plasmidfinder.tsv', sep = '\t')
# colnames(plasmidfinder) = plasmidfinder[1,]
# plasmidfinder = plasmidfinder[-1,]
# plasmidfinder = plasmidfinder[plasmidfinder$Database != 'Database',]
# plasmidfinder$Contig = gsub('plasmids.fastq_','',plasmidfinder$Contig)
# plasmidfinder$Sample = gsub('_contig.*','',plasmidfinder$Contig)
# 
# for(i in 1:nrow(plasmidfinder)){
#   if(plasmidfinder$Note[i] != ''){
#     plasmidfinder$Plasmid[i] = plasmidfinder$Note[i]
#   }else{
#     plasmidfinder$Plasmid[i] = plasmidfinder$Plasmid[i]
#   }
# }

# getting high-level plasmid groupings
# plasmidfinder$Group = gsub('rep.*','rep',plasmidfinder$Plasmid)
# plasmidfinder$Group = gsub('prgW.*','prgW',plasmidfinder$Group)
# plasmidfinder$Group = gsub('Inc.*','Inc',plasmidfinder$Group)
# plasmidfinder$Group = gsub('Col.*','Col',plasmidfinder$Group)
# plasmidfinder$Group = gsub('orf','ORF',plasmidfinder$Group)
# plasmidfinder$Plasmid = gsub('orf','ORF',plasmidfinder$Plasmid)
# groups_to_keep = c('rep','prgW','Inc','Col')
# plasmidfinder$Group[!plasmidfinder$Group %in% groups_to_keep] = 'Other'
# 
# # reading in AMRFinder results
# amrfinder = data.frame(read_excel('./data/p4_genome_amrfinder.xlsx'), stringsAsFactors = F)
# amrfinder = amrfinder[amrfinder$Name != 'Name',]
# amrfinder$Contig = paste0(amrfinder$Name,'_',amrfinder$Contig.id)
# amrfinder = amrfinder[amrfinder$Element.type == 'AMR',]
# 
# # reading in cluster labels
# cluster_labels = read.csv('./data/250818_p4_plasmid_blast_cluster_labels.csv', stringsAsFactors = F)
# ani_cluster_labels = read.csv('./data/250819_p4_plasmid_ANI_cluster_labels.csv', stringsAsFactors = F)

# reading in colors
colors = data.frame(read_excel('./notes/p4_colors.xlsx'), stringsAsFactors = F)

# ordering dataframes for plotting
tree$tip.label <- sub("^(([^_]*_){2}[^_]*)_", "\\1_contig_", tree$tip.label)  # Insert "contig" before third underscore
tree$tip.label <- sub("^([^_]*)_([^_]*)", "\\1.\\2", tree$tip.label)  # Replace first underscore with a period
#tree$tip.label <- sub("^(([^_]*_){2}[^_]*)_", "\\1_contig_", tree$tip.label)  # Insert "contig" before third underscore

ordered_df = data.frame(tree$tip.label, row.names = tree$tip.label)
colnames(ordered_df) = 'Tips'
ordered_df$Sample = gsub('_contig.*','', ordered_df$Tips)
ordered_df = ordered_df[ordered_df$Sample %in% mash_data$Sample,]
tree = drop.tip(tree, setdiff(tree$tip.label, ordered_df$Tips))


#write_lines(tree$tip.label, "/Users/gdlab/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/VIPtree/260128_tipsprefilter.txt" )


#Dropping phages that didn't pass pharokka and mapping

#tips_keep <- readLines("~/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/VIPtree/250919_pharokkapassedphages.txt")
tips_keep <- readLines("/Users/gdlab/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/VIPtree/260202_succesfullymappedphages_V2.txt")




# Transformations:
# 1. Replace first underscore with a period
# 2. Insert the word "chunk" between a double underscore ("__")
# 3. Replace @ with underscore

tips_keep_fmt <- tips_keep |>
  str_replace("_", ".") |>                    # only first underscore
  str_replace("__", "_contig_") |>           # insert chunk between double underscores
  str_replace_all("@", "_")                   # replace all @ with _

common <- intersect(tree$tip.label, tips_keep_fmt)
if (length(common) == 0) {
  stop("None of the formatted tip labels match the tree tips. Check transformations.")
}

# Keep only matching tips
tree <- keep.tip(tree, common)


#write_lines(tree$tip.label, "/Users/gdlab/Library/CloudStorage/Box-Box/DantasLabSriParuthiyil/DantasLab/Phage/plasmidome/250113_P4phagepipelinefinalizing/VIPtree/260128_tipspostfilter.txt" )



not_in_tree <- setdiff(tips_keep_fmt, tree$tip.label)
# 
# # Things in the tree that are NOT in tips_keep_fmt
not_in_tips <- setdiff(tree$tip.label, tips_keep_fmt)


#cluster_labels = cluster_labels[cluster_labels$Contig %in% tree$tip.label,]
for(i in 1:nrow(ordered_df)){
  ordered_df$Mash[i] = mash_pt_data$simplified_mash_call[mash_pt_data$Sample == ordered_df$Sample[i]]
  ordered_df$NICU[i] = mash_pt_data$Site[mash_pt_data$Sample == ordered_df$Sample[i]]
  ordered_df$abx_group_updated[i] = mash_pt_data$abx_group_updated[mash_pt_data$Sample == ordered_df$Sample[i]]
  
  # if(ordered_df$Sample[i] %in% plasmidfinder$Sample){
  # ordered_df$plasmid_group[i] = plasmidfinder$Group[plasmidfinder$Sample == ordered_df$Sample[i]]
  # ordered_df$plasmid[i] = plasmidfinder$Plasmid[plasmidfinder$Sample == ordered_df$Sample[i]]
  # }else{
  #   ordered_df$plasmid[i] = 'Unidentified'
  #   ordered_df$plasmid_group[i] = 'Unidentified'
  # }

  # if(ordered_df$Tips[i] %in% amrfinder$Contig){
  #   ordered_df$ARG[i] = amrfinder$Class[amrfinder$Contig == ordered_df$Tips[i]][1]
  # }else{
  #   ordered_df$ARG[i] = 'None'
  # }
  
  # if(ordered_df$Tips[i] %in% cluster_labels$Contig){
  #   ordered_df$Cluster[i] = cluster_labels$cluster[cluster_labels$Contig == ordered_df$Tips[i]]
  # }else{
  #   ordered_df$Cluster[i] = 'none'
  # }
  
  # if(ordered_df$Tips[i] %in% ani_cluster_labels$Contig){
  #   ordered_df$ani_cluster[i] = ani_cluster_labels$cluster[ani_cluster_labels$Contig == ordered_df$Tips[i]]
  # }else{
  #   ordered_df$ani_cluster[i] = 'none'
  # }
  
} 


# separating the dataframes
mash_df = ordered_df[,c('Tips','Mash')]
nicu_df = ordered_df[,'NICU', drop=F]
abx_df = ordered_df[,'abx_group_updated',drop=F]
mash_df$Mash <- sub("^(\\S+)\\s+.*$", "\\1 sp.", mash_df$Mash)
mash_df$Mash <- ifelse(
  grepl("^(Escherichia|Shigella)", mash_df$Mash),
  "Escherichia coli",
  mash_df$Mash
)

abx_df$abx_group_updated <- gsub("After", "after", abx_df$abx_group_updated)
#simple_plasmid_df = ordered_df[,'plasmid_group',drop=F]
#plasmid_df = ordered_df[,'plasmid',drop=F]
#arg_df = ordered_df[,c('Tips','ARG'),drop=F]
#cluster_df = ordered_df[,'Cluster',drop=F]
#ani_cluster_df = ordered_df[,'ani_cluster',drop=F]

# plotting the trees!!
mash_colors = colors$Hex[colors$Category == 'Simplified Mash']
names(mash_colors) = colors$Name[colors$Category == 'Simplified Mash']
# adding AMR class ring to tree
#arg_colors = c(natparks.pals('Torres',4),'0xFF')
#names(arg_colors) = c('AMINOGLYCOSIDE','BETA-LACTAM','LINCOSAMIDE','MACROLIDE','None')  


treeplot = ggtree(tree) %<+% mash_df
circle_treeplot = ggtree(tree, layout = 'fan') %<+% mash_df
# generating mash-tip treeplot
tip_treeplot = treeplot + geom_tiplab(align = TRUE, size = 0)+ geom_tippoint(aes(color = Mash), size = 3) +
  scale_color_manual(values = mash_colors, breaks = names(mash_colors))+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

circle_tip_treeplot = circle_treeplot + geom_tiplab(align = TRUE, size = 0)+ geom_tippoint(aes(color = Mash), size = 1.5) +
  scale_color_manual(values = mash_colors, breaks = names(mash_colors))+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

#arg_circle_tip = ggtree(tree, layout = 'fan') %<+% arg_df

# arg_circle_tip_treeplot = arg_circle_tip + geom_tiplab(align = TRUE, size = 0)+ geom_tippoint(aes(color = ARG), size = 3) +
#   scale_color_manual(values = arg_colors, breaks = names(arg_colors))+
#   theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

# plotting Mash as a ring instead
mash_df = mash_df[,2, drop =F]
# mash_tree = gheatmap(treeplot, mash_df, colnames_angle = 85, colnames_offset_y = -1,
#            offset = 0.0000001, width = 0.05, font.size = 5)+
#   geom_treescale(x = 30)+
#   scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+ggtree::vexpand(0.1,-1)+
#   theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'none')+
#   labs(fill = 'Bacterial host')
# 
# mash_tree2 = mash_tree + new_scale_fill()

circle_mash_tree <- gheatmap(circle_treeplot, mash_df,
                             colnames_angle = 85,
                             colnames_offset_y = -0.1,  # keep labels close
                             offset = 0.01,             # very small gap
                             width = 0.1,               # thinner ring
                             font.size = 5) +
  #geom_treescale(x = .3) +
  scale_fill_manual(values = mash_colors,
                    breaks = names(mash_colors)) +
  theme(legend.position = "right",
        plot.margin = margin(0,0,0,0)) +  # remove extra white space
  labs(fill = "Bacterial host")


circle_mash_tree2 = circle_mash_tree + new_scale_fill()

# adding abx ring
abx_colors = colors$Hex[colors$Category == 'Abx cohort']
names(abx_colors) = colors$Name[colors$Category == 'Abx cohort']

# mash_nicu_tree = gheatmap(mash_tree2, nicu_df, colnames_angle = 85, colnames_offset_y = -1,
#                           offset = 20, width = 0.05, font.size = 5)+
#   scale_fill_manual(values = nicu_colors, breaks = names(nicu_colors))+ggtree::vexpand(0.1,-1)+
#   theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
#         legend.position = 'right')+
#   labs(fill = 'NICU Location')
# 
# mash_nicu_tree2 = mash_nicu_tree + new_scale_fill()

circle_mash_nicu_tree <- gheatmap(circle_mash_tree2, abx_df,
                                  colnames_angle = 85,
                                  colnames_offset_y = -0.2,   # keep labels closer
                                  offset = 0.1,              # small spacing from previous ring
                                  width = 0.08,               # slimmer ring
                                  font.size = 5) +
  geom_treescale() +
  scale_fill_manual(values = abx_colors,
                    breaks = names(abx_colors)) +
  theme(legend.position = "right",        # hide legend completely
        plot.margin = margin(0,0,0,0)) + # remove extra whitespace
  labs(fill = "Infant Antibiotic Exposure")

circle_mash_nicu_tree2 = circle_mash_nicu_tree + new_scale_fill()

plot_file <- file.path('./figures/', paste0(Sys.Date(), "abx_mash_phage.pdf"))
#ggsave(filename = plot_file, plot = circle_mash_nicu_tree2, width = 18, height = 18)














