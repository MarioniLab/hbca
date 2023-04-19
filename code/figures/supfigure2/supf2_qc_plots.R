# Make QC supplementary plots.
# conda env preprocessing


### Load the required packages
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DropletUtils))

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(tidyverse))




### Inputs
mastertable_pw = '/nfs/research/marioni/areed/projects/hbca/metadata/2021-08-10_master_table.csv'
pD_pwd = '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/pD'
random_seed = 42



### Load data
mst_tbl <- as.data.frame(read.csv(file = mastertable_pw))
samples_list <- mst_tbl$sample_id

### Colour dictionaies
sample_type_colour_dictionary = c('Organoid unsorted' = '#fddd9d', 
                                  'Organoid LP sorted' = '#ff8d1f', 
                                  'Supernatant unsorted' = '#76e1f8', 
                                  'Supernatant live-sorted' = '#b3cd8f')

### Join pD together for all samples
pD_list <- list()
pD_joint <- data.frame()
for (smpl in samples_list){
  if (file.exists(paste0(pD_pwd, '/', smpl, '.pD.csv'))){
    pD_temp <- as.data.frame(read.csv(file = paste0(pD_pwd, '/', smpl, '.pD.csv')))
    pD_joint <- rbind(pD_joint, pD_temp)
  }
}
#make a clean pD
pD_joint_clean <- pD_joint[pD_joint$pass_all, ]

print(head(pD_joint)) 
print(dim(pD_joint)) #1651763, 15
print(dim(pD_joint_clean)) #880362, 15

### Make plots

## Directories
supfig2_dir = '/nfs/research/marioni/areed/projects/hbca/figures/supfigure2/'
dir.create(supfig2_dir, showWarnings = FALSE, recursive = TRUE)

## ggplots
umi_violin <- ggplot(pD_joint, aes(x=sample_id, y=sum, fill=sample_type)) +
  geom_violin() +
  scale_y_log10() +
  ylab("umi") +
  ggtitle("Total UMI Counts by Sample") +
  scale_fill_manual(values = sample_type_colour_dictionary) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, hjust=1))

detected_violin <- ggplot(pD_joint, aes(x=sample_id, y=detected, fill=sample_type)) +
  geom_violin() +
  scale_y_log10() +
  ylab("detected") +
  ggtitle("Total Genes by Sample") +
  scale_fill_manual(values = sample_type_colour_dictionary) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

mt_violin <- ggplot(pD_joint, aes(x=sample_id, y=subsets_Mito_percent, fill=sample_type)) +
  geom_violin() +
  scale_y_log10() +
  ylab("subsets_Mito_percent") +
  ggtitle("Mito Percent by Sample") +
  scale_fill_manual(values = sample_type_colour_dictionary) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))


## save plots
pdf(file = paste0(supfig2_dir, 'supfig2_preclean_umi_violins.pdf'), width = 35, height= 7)
print(umi_violin)
dev.off()

pdf(file = paste0(supfig2_dir, 'supfig2_preclean_detected_violins.pdf'), width = 35, height= 7)
print(detected_violin)
dev.off()

pdf(file = paste0(supfig2_dir, 'supfig2_preclean_mito_violins.pdf'), width = 35, height= 7)
print(mt_violin)
dev.off()

pdf(file = paste0(supfig2_dir, 'supfig2_clean_umi_violins.pdf'), width = 35, height= 7)
print(umi_violin %+% pD_joint_clean)
dev.off()

pdf(file = paste0(supfig2_dir, 'supfig2_clean_detected_violins.pdf'), width = 35, height= 7)
print(detected_violin %+% pD_joint_clean)
dev.off()

pdf(file = paste0(supfig2_dir, 'supfig2_clean_mito_violins.pdf'), width = 35, height= 7)
print(mt_violin %+% pD_joint_clean)
dev.off()




