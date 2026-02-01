# File name:      plasmid_nontree_plots.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-08-18
# Description:    This script generates non-tree plasmid figures for final p4 plasmids

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

#### reading in and prepping data ####
# reading in mash data
mash_data = data.frame(read_excel('./data/PretermPlasmidPhageProject_polypolish_circular_mash_checkm_HQgenomes.xlsx'), stringsAsFactors = F)

# reading in patient data
pts = data.frame(read_excel('./data/260131_P4_GM_dev_meta_updated_groupings_v1.xlsx'), stringsAsFactors = F)
colnames(pts)[colnames(pts) == 'Patient'] = 'Infant'
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

# reading in plasmidfinder data
plasmidfinder = read.table('./data/p4_plasmidfinder.tsv', sep = '\t')
colnames(plasmidfinder) = plasmidfinder[1,]
plasmidfinder = plasmidfinder[-1,]
plasmidfinder = plasmidfinder[plasmidfinder$Database != 'Database',]
plasmidfinder$Contig = gsub('plasmids.fastq_','',plasmidfinder$Contig)
plasmidfinder$Sample = gsub('_contig.*','',plasmidfinder$Contig)

for(i in 1:nrow(plasmidfinder)){
  if(plasmidfinder$Note[i] != ''){
    plasmidfinder$Plasmid[i] = plasmidfinder$Note[i]
  }else{
    plasmidfinder$Plasmid[i] = plasmidfinder$Plasmid[i]
  }
}

# getting high-level plasmid groupings
plasmidfinder$Group = gsub('rep.*','rep',plasmidfinder$Plasmid)
plasmidfinder$Group = gsub('prgW.*','prgW',plasmidfinder$Group)
plasmidfinder$Group = gsub('Inc.*','Inc',plasmidfinder$Group)
plasmidfinder$Group = gsub('Col.*','Col',plasmidfinder$Group)
plasmidfinder$Group = gsub('orf','ORF',plasmidfinder$Group)
plasmidfinder$Plasmid = gsub('orf','ORF',plasmidfinder$Plasmid)
groups_to_keep = c('rep','prgW','Inc','Col')
plasmidfinder$Group[!plasmidfinder$Group %in% groups_to_keep] = 'Other'

# reading in AMRFinder results
amrfinder = data.frame(read_excel('./data/p4_genome_amrfinder.xlsx'), stringsAsFactors = F)
amrfinder = amrfinder[amrfinder$Name != 'Name',]
amrfinder$Contig = paste0(amrfinder$Name,'_',amrfinder$Contig.id)
amrfinder = amrfinder[amrfinder$Element.type == 'AMR',]

# reading in cluster labels
cluster_labels = read.csv('./data/251130_p4_plasmid_ANI_cluster_labels.csv', stringsAsFactors = F)
cluster_labels$Sample = gsub('_contig.*','',cluster_labels$Contig)
cluster_labels = cluster_labels[cluster_labels$Sample %in% mash_pt_data$Sample,]
# reading in bed data
bed_data = data.frame(read_excel('./data/220712_PP_Hospital_bed_room_stay.xlsx'), stringsAsFactors = F)
colnames(bed_data)[3] = 'Infant'
bed_data$Infant = as.character(bed_data$Infant)

# putting data together for plotting
bed_data = bed_data[bed_data$Infant %in% mash_pt_data$Infant,]
mash_pt_bed = full_join(mash_pt_data, bed_data, relationship = 'many-to-many')

# reading in colors
colors = data.frame(read_excel('./notes/p4_colors.xlsx'), stringsAsFactors = F)

# dealing with bed data
bed_date = mash_pt_bed[,c('Infant','start','end')]
long_bed_date = melt(bed_date, id = 'Infant')
sample_bed_date = mash_pt_bed[,c('Sample','start','end')]
long_sample_beddate = melt(sample_bed_date, id = 'Sample')
colnames(long_bed_date)[1] = 'Sample'
colnames(long_sample_beddate)[1] = 'Sample'
# adding back on the rest of the data
for(i in 1:nrow(long_bed_date)){
  if(long_bed_date$variable[i] == 'start'){
    long_bed_date$bed[i] = mash_pt_bed$bed[mash_pt_bed$Infant == long_bed_date$Sample[i] & mash_pt_bed$start == long_bed_date$value[i]]
    long_bed_date$room[i] = mash_pt_bed$room[mash_pt_bed$Infant == long_bed_date$Sample[i] & mash_pt_bed$start == long_bed_date$value[i]]
  }else if(long_bed_date$variable[i] == 'end'){
    long_bed_date$bed[i] = mash_pt_bed$bed[mash_pt_bed$Infant == long_bed_date$Sample[i] & mash_pt_bed$end == long_bed_date$value[i]]
    long_bed_date$room[i] = mash_pt_bed$room[mash_pt_bed$Infant == long_bed_date$Sample[i] & mash_pt_bed$end == long_bed_date$value[i]]
  }
}

for(i in 1:nrow(long_sample_beddate)){
  if(long_sample_beddate$variable[i] == 'start'){
    long_sample_beddate$bed[i] = mash_pt_bed$bed[mash_pt_bed$Sample == long_sample_beddate$Sample[i] & mash_pt_bed$start == long_sample_beddate$value[i]]
    long_sample_beddate$room[i] = mash_pt_bed$room[mash_pt_bed$Sample == long_sample_beddate$Sample[i] & mash_pt_bed$start == long_sample_beddate$value[i]]
  }else if(long_sample_beddate$variable[i] == 'end'){
    long_sample_beddate$bed[i] = mash_pt_bed$bed[mash_pt_bed$Sample == long_sample_beddate$Sample[i] & mash_pt_bed$end == long_sample_beddate$value[i]]
    long_sample_beddate$room[i] = mash_pt_bed$room[mash_pt_bed$Sample == long_sample_beddate$Sample[i] & mash_pt_bed$end == long_sample_beddate$value[i]]
  }
}

infant_colors = sample(natparks.pals('Acadia',length(table(long_bed_date$Sample))))
names(infant_colors) = names(table(long_bed_date$Sample))

shared_plasmid_samples = cluster_labels$Sample[cluster_labels$cluster != 'none']
shared_plasmid_samples = gsub('_.*','',shared_plasmid_samples)
nonshared_plasmid_samples = cluster_labels$Sample[cluster_labels$cluster == 'none']
nonshared_plasmid_samples = gsub('_.*','',nonshared_plasmid_samples)

long_bed_date$shared = 'no'
long_bed_date$shared[long_bed_date$Sample %in% shared_plasmid_samples] = 'yes'

# getting a column of whether or not that cluster is shared between infants
cluster_labels$infant = gsub('_.*','',cluster_labels$Contig)
cluster_labels$Multi_infant = '1 infant'
for(i in 1:nrow(cluster_labels)){
  if(length(unique(cluster_labels$infant[cluster_labels$cluster == cluster_labels$cluster[i]])) > 1){
    cluster_labels$Multi_infant[i] = '>1 infant'
  }
}
cluster_labels$Multi_infant[cluster_labels$cluster == 'none'] = '1 infant'


rm2 = long_bed_date[long_bed_date$room == 'Room 2',]
rm4 = long_bed_date[long_bed_date$room == 'Room 4',]

# getting plasmid contigs
plasmid_contigs = read.csv('./data/plasmid_contigs.csv', stringsAsFactors = F)
plasmid_contigs$Sample = gsub('_contig.*','',plasmid_contigs$Contig)
plasmid_contigs$Sample = gsub('[>]','',plasmid_contigs$Sample)
plasmid_contigs = plasmid_contigs[plasmid_contigs$Sample %in% mash_data$Sample,]

mash_pt_data$plasmid = 'No plasmids'
mash_pt_data$plasmid[mash_pt_data$Sample %in% plasmid_contigs$Sample] = 'Plasmid'

#### plotting bed trace with plasmid cluster ####
cluster_beddate = full_join(long_sample_beddate, cluster_labels, relationship = 'many-to-many')
cluster_beddate$cluster[is.na(cluster_beddate$cluster)] = 'none'
cluster_beddate = cluster_beddate[!is.na(cluster_beddate$value),]

cluster_colors = colors$Hex[colors$Category == 'plasmid_cluster']
names(cluster_colors) = colors$Name[colors$Category == 'plasmid_cluster']

sharing_shapes = c(19,17)
names(sharing_shapes) = c('1 infant','>1 infant')

#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_bedtrace.pdf'), width = 14, height = 12)
ggplot(cluster_beddate[cluster_beddate$bed != 'KY' & cluster_beddate$bed != 'OK' & cluster_beddate$cluster != 'none',], 
       aes(x = value, y = bed, color = as.factor(cluster), shape = Multi_infant))+
  geom_point(size = 8, alpha = 0.95)+
  scale_shape_manual(values = sharing_shapes, breaks = names(sharing_shapes))+
  facet_wrap(~room, scales = 'free_y')+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Bed', color = 'Plasmid cluster',shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time and beds')
#dev.off()

cluster_beddate$infant = gsub('_.*','',cluster_beddate$Sample)
cluster_beddate$shared = 'no'
cluster_beddate$shared[cluster_beddate$infant %in% shared_plasmid_samples] = 'yes'
cluster_beddate$Multi_infant[cluster_beddate$cluster == 'none'] = '1 infant'
#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_time_OKKY.pdf'), width = 8, height = 4)
ggplot(cluster_beddate[cluster_beddate$room == 'KY' & cluster_beddate$cluster != 'none' | cluster_beddate$room == 'OK' & cluster_beddate$cluster != 'none',], aes(x = value, y = cluster, color = as.factor(cluster), 
                                                                                          shape = Multi_infant))+
  geom_point(size = 8, alpha = 0.95)+
  scale_shape_manual(values = sharing_shapes, breaks = names(sharing_shapes))+
  facet_wrap(~room, scales = 'free_y')+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Cluster', color = 'Plasmid cluster', shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time in OK and KY')
#dev.off()

cluster_rm2 = cluster_beddate[cluster_beddate$room == 'Room 2',]
cluster_rm4 = cluster_beddate[cluster_beddate$room == 'Room 4',]

#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_bedtrace_room2.pdf'), width = 6, height = 5)
ggplot(cluster_rm2[cluster_rm2$cluster != 'none',], aes(x = value, y = bed, color = as.factor(cluster),
                                                        shape = Multi_infant))+
  geom_point(size = 16, alpha = 0.95)+
  facet_wrap(~room, scales = 'free_y')+
  scale_shape_manual(values = sharing_shapes, breaks = names(sharing_shapes))+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Bed', color = 'Plasmid cluster',shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time and beds')
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_bedtrace_room4.pdf'), width = 6, height = 5)
ggplot(cluster_rm4[cluster_rm4$cluster != 'none',], aes(x = value, y = bed, color = as.factor(cluster),
                                                        shape = Multi_infant))+
  geom_point(size = 16, alpha = 0.95)+
  facet_wrap(~room, scales = 'free_y')+
  scale_shape_manual(values = sharing_shapes, breaks = names(sharing_shapes))+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Bed', color = 'Plasmid cluster', shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time and beds')
#dev.off()


#### plasmids per isolate per genus per site plots ####

nicu_colors = colors$Hex[colors$Category == 'NICU location']
names(nicu_colors) = colors$Name[colors$Category == 'NICU location']



# plotting number of plasmids per isolate by genus
plasmid_contigs = plasmid_contigs[!duplicated(plasmid_contigs$Contig),]
plasmid_counts = data.frame(table(plasmid_contigs$Sample), stringsAsFactors = F)
colnames(plasmid_counts) = c('Sample','Plasmid_contigs')
for(i in 1:nrow(mash_pt_data)){
  if(mash_pt_data$Sample[i] %in% plasmid_counts$Sample){
    mash_pt_data$Plasmid_counts[i] = plasmid_counts$Plasmid_contigs[plasmid_counts$Sample == mash_pt_data$Sample[i]]
  }else{
    mash_pt_data$Plasmid_counts[i] = 0
  }
}

mash_pt_data$genus = gsub(' .*','',mash_pt_data$simplified_mash_call)
mash_pt_data$genus = gsub('Shigella','Escherichia',mash_pt_data$genus)
mash_pt_data$abx_group = gsub(' After',' after',mash_pt_data$abx_group)

abx_colors = colors$Hex[colors$Category == 'Abx cohort']
names(abx_colors) = colors$Name[colors$Category == 'Abx cohort']

#pdf(paste0('./figures/', Sys.Date(),'_plasmids_per_isolate_by_genus.pdf'), width = 16, height = 8)
ggplot(mash_pt_data, aes(x = factor(abx_group, levels = c('None','Early only','Early and after')), 
                         y = Plasmid_counts, fill = abx_group))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  facet_grid(~genus, scales = 'free')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of plasmids per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(axis.text.x = element_text(size = 12), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of plasmids per isolate')
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_plasmids_per_isolate_by_genus.pdf'), width = 16, height = 12)
ggplot(mash_pt_data, aes(x = Site, 
                         y = Plasmid_counts, fill = Site))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = Site), size = 2, color = 'grey10')+
  facet_grid(~genus, scales = 'free')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of plasmids per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of plasmids per isolate')
#dev.off()

genus_colors = colors$Hex[colors$Category == 'genus']
names(genus_colors) = colors$Name[colors$Category == 'genus']

mash_colors = colors$Hex[colors$Category == 'Simplified Mash']
names(mash_colors) = colors$Name[colors$Category == 'Simplified Mash']

#pdf(paste0('./figures/', Sys.Date(),'_plasmids_per_isolate_by_genus.pdf'), width = 16, height = 12)
ggplot(mash_pt_data, aes(x = genus, y = Plasmid_counts, fill = genus, color = simplified_mash_call))+
  geom_boxplot(alpha = 0.7,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = genus, color = simplified_mash_call), size = 3)+
  labs(fill = 'Genus',color = 'Species', y = 'Number of plasmids per isolate', x = 'Genus')+
  theme_bw()+
  scale_fill_manual(values = genus_colors, breaks = names(genus_colors))+
  scale_color_manual(values=mash_colors, breaks = names(mash_colors))+
  theme(axis.text.x = element_text(size = 12), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of plasmids per isolate')
#dev.off()

#pdf(paste0('./figures/', Sys.Date(),'_plasmids_per_isolate.pdf'), width = 16, height = 12)
ggplot(mash_pt_data, aes(x = factor(abx_group, levels = c('Early and after','Early only','None')), 
                         y = Plasmid_counts, fill = abx_group, color = genus))+
  geom_boxplot(alpha = 0.8,outliers = F, color = 'grey10')+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), aes(group = abx_group), size = 2, color = 'grey10')+
  labs(fill = 'Antibiotic exposure',color = 'Antibiotic exposure', y = 'Number of plasmids per isolate', x = 'Antibiotic exposure')+
  theme_bw()+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  scale_color_manual(values=genus_colors, breaks = names(genus_colors))+
  theme(axis.text.x = element_text(size = 12), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 14, face = 'italic'))+
  ggtitle('Number of plasmids per isolate')
#dev.off()


#### Preparing data for for plasmid cluster visualizations in Gephi ####
fastani = read.delim('./data/p4_plasmid_fastani_out.txt', header = F)
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
fastani = fastani[fastani$source_sample %in% mash_data$Sample,]
fastani = fastani[fastani$target_sample %in% mash_data$Sample,]

fastani = fastani[fastani$ANI > 99.9,]
fastani = fastani[fastani$source != fastani$target,]
fastani = fastani[abs(fastani$sourcelen - fastani$targetlen) < 0.1*pmax(fastani$sourcelen, fastani$targetlen),]

cluster_labels = read.csv('./data/251130_p4_plasmid_ANI_cluster_labels.csv', stringsAsFactors = F)
cluster_labels$Sample = gsub('_contig.*','',cluster_labels$Contig)
#write_tsv(fastani,'./251011_plasmid_gephi/251130_fastani_sameplasmids.tsv')
gephi_nodes = mash_pt_data
gephi_nodes$genus = gsub(' .*','',gephi_nodes$simplified_mash_call)
gephi_nodes$genus[gephi_nodes$genus == 'Shigella'] = 'Escherichia'
colnames(gephi_nodes)[1] = 'id'
gephi_nodes = gephi_nodes[!duplicated(gephi_nodes$id),]
gephi_nodes = gephi_nodes[gephi_nodes$id %in% fastani$source_sample | gephi_nodes$id %in% fastani$target_sample,]
#write.csv(gephi_nodes, './251011_plasmid_gephi/251130_fastani_nodes_table.csv', quote = F, row.names = F)

#### genome fastANI results ####
genome_fastani = read.table('./data/251117_p4_genome_fastani_out.txt', sep = '\t')
colnames(genome_fastani) = c('source','target','ANI','sourcelen','targetlen')
genome_fastani$source = gsub('.*polished_assemblies/','',genome_fastani$source)
genome_fastani$source = gsub('_polished.*','',genome_fastani$source)
genome_fastani$target = gsub('.*polished_assemblies/','',genome_fastani$target)
genome_fastani$target = gsub('_polished.*','',genome_fastani$target)
genome_fastani = genome_fastani[genome_fastani$source != genome_fastani$target,]

genome_fastani = genome_fastani[abs(genome_fastani$sourcelen - genome_fastani$targetlen) < 0.1*pmax(genome_fastani$sourcelen, genome_fastani$targetlen),]

clustered_plasmids_genome_fastani = genome_fastani[genome_fastani$source %in% fastani$source_sample |
                                                     genome_fastani$target %in% fastani$source_sample |
                                                     genome_fastani$source %in% fastani$target_sample |
                                                     genome_fastani$target %in% fastani$target_sample,]
annotated_ani = read.csv('./data/251117_plasmids_fastani_clusters_genomesani.csv', stringsAsFactors = F)
fastani = full_join(fastani, annotated_ani)

for(i in 1:nrow(cluster_beddate)){
  if('no' %in% fastani$highly_similar_sources[fastani$source == cluster_beddate$Contig[i] | fastani$target == cluster_beddate$Contig[i]]){
    cluster_beddate$same_isolate_plasmid[i] = 'distinct isolates' 
  }else{
    cluster_beddate$same_isolate_plasmid[i] = 'highly similar isolates'
  }
}

cluster_beddate$Multi_infant_Multi_isolate = paste0(cluster_beddate$Multi_infant, ', ', cluster_beddate$same_isolate_plasmid)


# plotting bedtrace with this information
multi_sharing_shapes = c(1,19,2,17)
names(multi_sharing_shapes) = c('1 infant, highly similar isolates','1 infant, distinct isolates','>1 infant, highly similar isolates',
                                '>1 infant, distinct isolates')

#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_bedtrace.pdf'), width = 12, height = 9)
ggplot(cluster_beddate[cluster_beddate$cluster != 'none' & cluster_beddate$bed != 'KY' & cluster_beddate$bed != 'OK',], aes(x = value, y = bed, color = as.factor(cluster),
                                                        shape = Multi_infant_Multi_isolate))+
  geom_point(size = 8, alpha = 0.95)+
  facet_wrap(~room, scales = 'free_y')+
  scale_shape_manual(values = multi_sharing_shapes, breaks = names(multi_sharing_shapes))+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Bed', color = 'Plasmid cluster', shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time and beds')
#dev.off()

nonone_clusterbed = cluster_beddate[cluster_beddate$cluster != 'none',]

#pdf(paste0('./figures/',Sys.Date(),'_p4_plasmidclusters_bedtrace_OKKY.pdf'), width = 10, height = 4)
ggplot(nonone_clusterbed[nonone_clusterbed$bed == 'KY' | nonone_clusterbed$bed == 'OK',], aes(x = value, y = cluster, color = as.factor(cluster),
                                                                                                                            shape = Multi_infant_Multi_isolate))+
  geom_point(size = 8, alpha = 0.95)+
  facet_wrap(~room, scales = 'free_y')+
  scale_shape_manual(values = multi_sharing_shapes, breaks = names(multi_sharing_shapes))+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Date', y = 'Bed', color = 'Plasmid cluster', shape = 'Infants with plasmid cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 12, face = 'bold'),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'),
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid clusters across time and beds')
#dev.off()

# getting N plasmids found in multiple rooms
rooms = cluster_beddate
rooms$clusterroom = paste0(rooms$cluster, rooms$room)
rooms = rooms[!duplicated(rooms$clusterroom),]
table(rooms$cluster[rooms$bed != 'KY' & rooms$bed != 'OK'])
site = rooms
site$site = 'MO'
site$site[site$bed == 'KY'] = 'KY'
site$site[site$bed == 'OK'] = 'OK'
site$clustersite = paste0(site$cluster, site$site)
site = site[!duplicated(site$clustersite),]
table(site$cluster)

#### bakta for shared plasmids ####
# reading in bakta results
plasmid_bakta = read.delim('./data/250914_plasmid_bakta.tsv')
plasmid_bakta$Contig = gsub('[.]tsv.*','',plasmid_bakta$Contig)
plasmid_bakta = full_join(plasmid_bakta, cluster_labels)
plasmid_bakta = plasmid_bakta[!is.na(plasmid_bakta$Product),]
plasmid_bakta = plasmid_bakta[plasmid_bakta$Product != 'hypothetical protein',]
plasmid_bakta$Sample = gsub('_contig.*','',plasmid_bakta$Contig)
plasmid_bakta$UniRef50 = gsub('.*UniRef50_','UniRef50_',plasmid_bakta$DbRef)
plasmid_bakta$UniRef50 = gsub(',.*','',plasmid_bakta$UniRef50)
plasmid_bakta$UniRef90 = gsub('.*UniRef90_','UniRef90_',plasmid_bakta$DbRef)
plasmid_bakta$UniRef90 = gsub(',.*','',plasmid_bakta$UniRef90)
plasmid_bakta = plasmid_bakta[plasmid_bakta$Sample %in% mash_pt_data$Sample,]
plasmid_bakta$cluster[is.na(plasmid_bakta$cluster)] = 'none'
plasmid_bakta = left_join(plasmid_bakta, mash_pt_bed)
plasmid_bakta = plasmid_bakta[!is.na(plasmid_bakta$Contig),]

plasmid_bakta$Product = gsub(',','--',plasmid_bakta$Product)
#write.csv(plasmid_bakta[,c(1:8,10)], './data/251208_plasmid_bakta.csv', row.names = F, quote = F)


# getting functions in shared plasmids vs not in shared plasmids
shared_only_product = plasmid_bakta[plasmid_bakta$Product %in% 
                      setdiff(plasmid_bakta$Product[plasmid_bakta$cluster != 'none'],plasmid_bakta$Product[plasmid_bakta$cluster == 'none']),]
nonshared_only_product = plasmid_bakta[plasmid_bakta$Product %in% 
                         setdiff(plasmid_bakta$Product[plasmid_bakta$cluster == 'none'],plasmid_bakta$Product[plasmid_bakta$cluster != 'none']),]

shared_only_gene = plasmid_bakta[plasmid_bakta$Gene %in% 
                   setdiff(plasmid_bakta$Gene[plasmid_bakta$cluster != 'none'],plasmid_bakta$Gene[plasmid_bakta$cluster == 'none']),]
nonshared_only_gene = plasmid_bakta[plasmid_bakta$Gene %in% 
                      setdiff(plasmid_bakta$Gene[plasmid_bakta$cluster == 'none'],plasmid_bakta$Gene[plasmid_bakta$cluster != 'none']),]

#pdf(paste0('./figures/',Sys.Date(),'_p4_cluster8_plasmid_bakta.pdf'), width = 10, height  = 18)
ggplot(plasmid_bakta[plasmid_bakta$cluster == '8' & plasmid_bakta$Contig != '253.03_EF1_contig_7',], aes(x = Contig, y = Product, color = as.factor(cluster), fill = as.factor(cluster)))+
  geom_tile(size = 5)+
  #facet_wrap(~room, scales = 'free_y')+
  scale_color_manual(values = cluster_colors, breaks = names(cluster_colors))+
  scale_fill_manual(values = cluster_colors, breaks = names(cluster_colors))+
  theme_bw()+
  labs(x = 'Plasmid', y = 'Bakta annotation', color = 'Cluster')+
  theme(strip.text = element_text(size = 14, face = 'bold'), axis.text = element_text(size = 9),
        axis.title = element_text(size = 14, face = 'bold'),
        title = element_text(size = 14, face = 'bold'), legend.position = 'none',
        plot.margin = margin(19,19,19))+
  ggtitle('P4 plasmid cluster 8')
#dev.off()


