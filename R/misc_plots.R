# File name:      misc_plots.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-08-18
# Description:    This script generates non-tree and non-plasmid figures

# loading packages
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(NatParksPalettes)
library(RColorBrewer)
library(lubridate)
library(kableExtra)
library(ggpubr)
library(qgraph)
library(KEGGREST)


# setwd
setwd('~/Box Sync/EEB/Dantas/pppp/')

#### reading in data ####

# reading in mash data
mash_data = data.frame(read_excel('./data/PretermPlasmidPhageProject_polypolish_circular_mash_checkm_HQgenomes.xlsx'), stringsAsFactors = F)

# reading in patient data
pts = data.frame(read_excel('./data/250220_P4_GM_Dev_Abx_groupings_updated.xlsx'), stringsAsFactors = F)
colnames(pts)[2] = 'Infant'
class(pts$Infant) = 'character'
pts = pts[pts$Cohort == 'P4_cross_sectional',]
pts$Site = gsub('St. Louis','MO',pts$Site)
pts$Site = gsub('Oklahoma','OK',pts$Site)
pts$Site = gsub('Louisville','KY',pts$Site)
pts = pts[pts$Infant %in% mash_data$Infant,]
pts$abx_group = gsub('A','a',pts$abx_group)
class(mash_data$Infant) = 'character'

mash_pt_data = full_join(mash_data, pts)
mash_pt_data = mash_pt_data[!is.na(mash_pt_data$Sample),]

plasmid_contigs = read.csv('./data/plasmid_contigs.csv', stringsAsFactors = F)
plasmid_contigs$Sample = gsub('_contig.*','',plasmid_contigs$Contig)
plasmid_contigs$Sample = gsub('[>]','',plasmid_contigs$Sample)
plasmid_contigs = plasmid_contigs[plasmid_contigs$Sample %in% mash_data$Sample,]

mash_pt_data$plasmid = 'No plasmids'
mash_pt_data$plasmid[mash_pt_data$Sample %in% plasmid_contigs$Sample] = 'Plasmid'

# reading in AMRFinder results
amrfinder = data.frame(read_excel('./data/p4_genome_amrfinder.xlsx'), stringsAsFactors = F)
amrfinder = amrfinder[amrfinder$Name != 'Name',]
amrfinder$Contig = paste0(amrfinder$Name,'_',amrfinder$Contig.id)
amrfinder = amrfinder[amrfinder$Element.type == 'AMR',]
amrfinder$sampgene = paste0(amrfinder$Name, amrfinder$Gene.symbol)
amrfinder$contig = paste0(amrfinder$Name, '_',amrfinder$Contig.id)
amrfinder$chromplasmid = 'chromosome'
plasmid_contigs$Contig = gsub('>','',plasmid_contigs$Contig)
amrfinder$chromplasmid[amrfinder$contig %in% plasmid_contigs$Contig] = 'plasmid'
unique_amrfinder = amrfinder[!duplicated(amrfinder$sampgene),]
for(i in 1:nrow(mash_pt_data)){
  mash_pt_data$args_per_isolate[i] = sum(unique_amrfinder$Name == mash_pt_data$Sample[i])
  mash_pt_data$all_args[i] = sum(amrfinder$Name == mash_pt_data$Sample[i])
  mash_pt_data$chrom_args[i] = sum(unique_amrfinder$chromplasmid == 'chromosome' & unique_amrfinder$Name == mash_pt_data$Sample[i])
  mash_pt_data$all_chrom_args[i] = sum(amrfinder$chromplasmid == 'chromosome' & amrfinder$Name == mash_pt_data$Sample[i])
  mash_pt_data$plasmid_args[i] = sum(unique_amrfinder$chromplasmid == 'plasmid' & unique_amrfinder$Name == mash_pt_data$Sample[i])
  mash_pt_data$all_plasmid_args[i] = sum(amrfinder$chromplasmid == 'plasmid' & amrfinder$Name == mash_pt_data$Sample[i])
}

# reading in cluster labels
cluster_labels = read.csv('./data/251130_p4_plasmid_ANI_cluster_labels.csv', stringsAsFactors = F)
cluster_labels$Sample = gsub('_contig.*','',cluster_labels$Contig)
cluster_labels = cluster_labels[cluster_labels$Sample %in% mash_pt_data$Sample,]

# reading in colors
colors = data.frame(read_excel('./notes/p4_colors.xlsx'), stringsAsFactors = F)

# reading in metaphlan data
met_data = read.delim('./data/251013_P4_crosssectional_merged_metaphlan3.txt', sep = '\t')
met_data$clade_name = gsub('.*s__','',met_data$clade_name)
met_data = met_data[!grepl('k__',met_data$clade_name),]
met_data$clade_name = gsub('_',' ',met_data$clade_name)
colnames(met_data) = gsub('X','',colnames(met_data))
colnames(met_data) = gsub('_paired_metaphlan','',colnames(met_data))
t_metdata = data.frame(t(met_data))
colnames(t_metdata) = t_metdata[1,]
t_metdata = t_metdata[-1,]
t_metdata$Infant = gsub('_','.',rownames(t_metdata))
# long
long_t_metdata = melt(t_metdata, id = 'Infant')
# subsetting to remove infants who are not included
long_t_metdata = long_t_metdata[long_t_metdata$Infant %in% mash_data$Infant,]
# adding on abx group and site
abxgroup_site = mash_pt_data[,c('Infant','abx_group','Site')]
abxgroup_site = abxgroup_site[!duplicated(abxgroup_site),]
long_t_metdata = left_join(long_t_metdata, abxgroup_site)
long_t_metdata$genus = gsub(' .*','', long_t_metdata$variable)

genus_pal = c(colors$Hex[colors$Category == 'genus'],'grey85')
names(genus_pal) = c(colors$Name[colors$Category == 'genus'],'Non-target taxa')

mash_colors = c(colors$Hex[colors$Category == 'Simplified Mash'], 'grey85')
names(mash_colors) = c(colors$Name[colors$Category == 'Simplified Mash'],'Non-target taxa')

mash_data$genus = gsub(' .*','',mash_data$simplified_mash_call)
long_t_metdata$variable = gsub('Enterobacter .*','Enterobacter sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' oxytoca',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' michiganensis',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' aerogenes',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' quasipneumoniae',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' variicola',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' avium',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' casseliflavus',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' gallinarum',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' hominis',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' warneri',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' durans',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' gilvus',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' hirae',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' sp ESNIH1',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' argenteus',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' caprae',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' schweitzeri',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' haemolyticus',' sp.',long_t_metdata$variable)
long_t_metdata$variable = gsub(' lugdunensis',' sp.',long_t_metdata$variable)

long_t_metdata$variable[!(long_t_metdata$variable %in% mash_data$simplified_mash_call)] = 'Non-target taxa'
long_t_metdata$genus[!(long_t_metdata$genus %in% mash_data$genus)] = 'Non-target taxa'

class(long_t_metdata$value) = 'numeric'

# reading in antibiotic exposure data for large cohort
abx_data = data.frame(read_excel('./data/260131_P4_GM_dev_meta_updated_groupings_v1.xlsx'), stringsAsFactors = F)

#### making cultured isolates barplot - Figure 3B ####
#pdf(paste0('./figures/',Sys.Date(),'_p4_cultured_species_barplot.pdf'), width = 14, height = 9)
ggplot(mash_pt_data, aes(x = simplified_mash_call, fill = simplified_mash_call))+
  geom_bar(stat = 'count')+
  facet_wrap(~Site, scales = 'free')+
  scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+
  theme_bw()+
  labs(x = 'Taxa', y = 'Count of isolates', fill = 'Taxa')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('Cross-sectional cultured isolates')
#dev.off()


#pdf(paste0('./figures/',Sys.Date(),'_p4_cultured_species_barplot.pdf'), width = 14, height = 9)
ggplot(mash_pt_data, aes(x = simplified_mash_call, fill = simplified_mash_call))+
  geom_bar(stat = 'count')+
  facet_wrap(~plasmid, scales = 'free')+
  scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+
  theme_bw()+
  labs(x = 'Taxa', y = 'Count of isolates', fill = 'Taxa')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('Cross-sectional cultured isolates')
#dev.off()

mash_pt_data$genus = gsub(' .*','',mash_pt_data$simplified_mash_call)
mash_pt_data$genus[mash_pt_data$genus == 'Shigella'] = 'Escherichia'
#pdf(paste0('./figures/',Sys.Date(),'_p4_cultured_species_barplot.pdf'), width = 14, height = 9)
ggplot(mash_pt_data, aes(x = factor(genus, levels = c('Enterococcus','Klebsiella','Escherichia','Staphylococcus','Enterobacter')), 
                         fill = simplified_mash_call))+
  geom_bar(stat = 'count')+
  geom_text(
    aes(label = after_stat(count)),
    stat = "count",
    position = position_stack(vjust = 0.5),  # Place text in the middle of each segment
    color = "grey85", size = 5  # Adjust label size and color
  )+
  scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+
  theme_bw()+
  labs(x = 'Genus', y = 'Count of isolates', fill = 'Taxa')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('Cross-sectional cultured isolates')
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_p4_cultured_species_barplot_abxgroup.pdf'), width = 14, height = 9)
ggplot(mash_pt_data, aes(x = factor(genus, levels = c('Enterococcus','Klebsiella','Escherichia','Staphylococcus','Enterobacter','Unidentified')), 
                         fill = simplified_mash_call))+
  geom_bar(stat = 'count')+
  facet_wrap(~abx_group, scales = 'free')+
  geom_text(
    aes(label = after_stat(count)),  # Label by MLST
    stat = "count",
    position = position_stack(vjust = 0.5),  # Place text in the middle of each segment
    color = "grey85", size = 5  # Adjust label size and color
  )+
  scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+
  theme_bw()+
  labs(x = 'Genus', y = 'Count of isolates', fill = 'Taxa')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('Cross-sectional cultured isolates')
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_p4_cultured_species_barplot_NICU.pdf'), width = 14, height = 9)
ggplot(mash_pt_data, aes(x = factor(genus, levels = c('Enterococcus','Klebsiella','Escherichia','Staphylococcus','Enterobacter','Unidentified')), 
                         fill = simplified_mash_call))+
  geom_bar(stat = 'count')+
  facet_wrap(~Site, scales = 'free')+
  geom_text(
    aes(label = after_stat(count)),  # Label by MLST
    stat = "count",
    position = position_stack(vjust = 0.5),  # Place text in the middle of each segment
    color = "grey85", size = 5  # Adjust label size and color
  )+
  scale_fill_manual(values = mash_colors, breaks = names(mash_colors))+
  theme_bw()+
  labs(x = 'Genus', y = 'Count of isolates', fill = 'Taxa')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('Cross-sectional cultured isolates')
#dev.off()

#### cross-sectional cohort MetaPhlAn plot - Figure 3C ####
## doing an average for each abx group
for(i in 1:nrow(long_t_metdata)){
  long_t_metdata$group_avg[i] = 100*(sum(long_t_metdata$value[long_t_metdata$variable == long_t_metdata$variable[i] & long_t_metdata$abx_group == long_t_metdata$abx_group[i]]) / 
                                       sum(long_t_metdata$value[long_t_metdata$abx_group == long_t_metdata$abx_group[i]]))
}

avg_metaphlan = long_t_metdata
avg_metaphlan$specgroup = paste0(avg_metaphlan$variable, avg_metaphlan$abx_group)
avg_metaphlan = avg_metaphlan[!duplicated(avg_metaphlan$specgroup),]

#pdf(paste0('./figures/', Sys.Date(),'_p4_crossectional_metaphlan3_avgs.pdf'), width = 9, height = 7)
ggplot(avg_metaphlan, aes(x = abx_group, y = group_avg, fill = 
                            factor(variable, levels = c('Enterobacter sp.','Enterococcus faecalis','Enterococcus faecium','Enterococcus sp.',
                                                        'Escherichia coli','Klebsiella pneumoniae','Klebsiella sp.','Staphylococcus aureus',
                                                        'Staphylococcus epidermidis','Staphylococcus sp.','Non-target taxa'))))+
  geom_bar(position = 'stack', stat = 'identity')+
  labs(fill = 'Species', y = 'MetaPhlAn3 abundance results', x = 'Antibiotic exposure group')+
  theme_bw()+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100.05))+
  scale_fill_manual(values = mash_colors)+
  theme(axis.text.x = element_text(size = 14),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'bold'))+
  ggtitle('MetaPhlAn3 results for p4 cohort')
#dev.off()


#### AMRFinder ARGs per isolate plots - Supplementary Figure 5A-E ####
abx_colors = colors$Hex[colors$Category == 'Abx cohort']
names(abx_colors) = colors$Name[colors$Category == 'Abx cohort']

#pdf(paste0('./figures/', Sys.Date(),'_uniqueARGs_per_isolate.pdf'), width = 9, height = 6)
ggplot(mash_pt_data, aes(x = abx_group, y = args_per_isolate, fill = abx_group))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of unique ARGs per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 14), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of unique ARGs per isolate')
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_uniqueARGs_per_isolate_chromosome.pdf'), width = 9, height = 6)
ggplot(mash_pt_data, aes(x = abx_group, y = chrom_args, fill = abx_group))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of unique ARGs per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 14), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of unique ARGs per isolate on chromosome')
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_uniqueARGs_per_isolate_plasmid.pdf'), width = 9, height = 6)
ggplot(mash_pt_data, aes(x = abx_group, y = plasmid_args, fill = abx_group))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of unique ARGs per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 14), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of unique ARGs per isolate on plasmids')
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_allARGs_per_isolate.pdf'), width = 9, height = 6)
ggplot(mash_pt_data, aes(x = abx_group, y = all_args, fill = abx_group))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of ARGs per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 14), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of ARGs per isolate')
#dev.off()


#### Antibiotic exposure by day of life - Supplementary Figure 1D ####
#pdf(paste0('./figures/', Sys.Date(),'_largecohort_acuteabxexposure_by_dol.pdf'), width = 15)
ggplot(abx_data[abx_data$acute_abx_exposure_count > 0,], aes(x = DOL, fill = abx_group_updated, color = abx_group_updated))+
  geom_histogram(stat = 'count', alpha = 0.6)+
  labs(fill = 'Antibiotic exposure group',color = 'Antibiotic exposure group', y = 'Number of samples', x = 'Day of life')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  scale_color_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 12),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_largecohort_acuteabxexposure_by_dol_earlyonly.pdf'), width = 15)
ggplot(abx_data[abx_data$acute_abx_exposure_count > 0 & abx_data$abx_group == 'Early only',], aes(x = DOL, fill = abx_group, color = abx_group))+
  geom_histogram(stat = 'count', alpha = 0.6)+
  labs(fill = 'Antibiotic exposure group',color = 'Antibiotic exposure group', y = 'Number of samples', x = 'Day of life')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  scale_color_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 12),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_largecohort_acuteabxexposure_by_dol_earlyandafter.pdf'), width = 15)
ggplot(abx_data[abx_data$acute_abx_exposure_count > 0 & abx_data$abx_group == 'Early and after',], aes(x = DOL, fill = abx_group, color = abx_group))+
  geom_histogram(stat - Supp= 'count', alpha = 0.6)+
  labs(fill = 'Antibiotic exposure group',color = 'Antibiotic exposure group', y = 'Number of samples', x = 'Day of life')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  scale_color_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 12),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))
#dev.off()

