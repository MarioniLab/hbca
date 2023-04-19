#! Make plots of gene expression (from dge ARvsHR) across celltypes.
#scVI new
#C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/figures/celltype_markers/dge_gene_violins.R


suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(Matrix))

suppressMessages(library(optparse))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(ggrastr))
suppressMessages(library(ggbreak))
suppressMessages(library(reshape2))








input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi_2/output'
sce <- readRDS(paste0(input_path, '/input/input_sce.rds'))
rownames(sce) <- rowData(sce)$X



###PDL1

#violin plots over epithelial clusters coloured by tissue_condition
sce_sub <- sce[, sce$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4', 'LP5',
                                   'HS1', 'HS2', 'HS3', 'HS4', 
                                   'BSL1', 'BSL2')]

colour_set_tissue_condition <- c('Mammoplasty WT' = '#f1ce63', 'Mastectomy BRCA1' = '#4e79a7',
                                 'Mastectomy BRCA2' = '#a0cbe8', 'HR-unk' = '#bab0ac')
colour_set_level1 <- c('Luminal progenitor' = '#DDA0DD', 'Luminal hormone sensing' = '#EE8298', 'Basal' = '#E6554A')
level2_colour_dictionary = c('LP1' = "#DA80DA", 'LP2' = "#815481", 'LP3' = "#C040C0", 'LP4' = "#E1AFE1", 'LP5' = "#3F0034",
                             'HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028",
                             'BSL1' = "#EB675E", 'BSL2' = "#A23E36",
                             'DDC1' = "#540F54", 'DDC2' = "#53407F",
                             'FB1' = "#DFA38A", 'FB2' = "#8C3612", 'FB3' = "#623623", 'FB4' = "#916350", 'FB5' = "#DAC3C3",
                             'VM1' = "#F8770B", 'VM2' = "#E09E3A", 'VM3' = "#CD7225", 'VM4' = "#FFC990", 'VM5' = "#AC5812",
                             'EC venous' = "#FEE083", 'EC capillary' = "#897538", 'EC arterial' = "#E7B419", 'EC angiogenic tip' = "#BCA048",
                             'LEC1' = "#6F8BE2", 'LEC2' = "#3053BC",
                             'CD8T 1' = "#6D9F58", 'CD8T 2' = "#9EB766", 'CD8T 3' = "#BDCB10", 'CD4T' = "#3A6527", 'IFNG+ T' = "#9EA743",
                             'NK1' = "#E2E8A7", 'NK2' = "#5A6209", 'NK3' = "#8FE36B",
                             'ILC3' = "#818A31",
                             'B cell' = "#9FC5E8", 'Plasma cell' = "#23D9F1",
                             'Macrophage' = "#64C6A6") #ignore doublets
pD <- colData(sce_sub)
pD$CD274 <- logcounts(sce_sub)['CD274',]
#pD$logCD274 <- log2(pD$CD274+1)
pD <- as.data.frame(pD)

# # gg <- plotExpression(sce_sub,
# #                      features='CD274',
# #                      x='level2',
# #                      colour_by = 'tissue_condition') +
# #   scale_colour_manual(values = colour_set_tissue_condition)
# # gg2 <- plotExpression(sce_sub,
# #                      features='CD274',
# #                      x='level2',
# #                      colour_by = 'tissue_condition',
# #                      log2_values=TRUE) +
# #   scale_colour_manual(values = colour_set_tissue_condition)
# gg <- ggplot(pD) +
#   geom_violin(mapping = aes_string(x='level2', y='CD274', fill='level1')) +
#   geom_beeswarm_rast(mapping = aes_string(x='level2', y='CD274', colour='tissue_condition')) +
#   scale_colour_manual(values = colour_set_tissue_condition) +
  # scale_fill_manual(values = colour_set_level1) +
  # theme_classic()
# gg2 <- ggplot(pD) +
#   geom_violin(mapping = aes_string(x='level2', y='logCD274', fill='level1')) +
#   geom_beeswarm_rast(mapping = aes_string(x='level2', y='CD274', colour='tissue_condition')) +
#   scale_colour_manual(values = colour_set_tissue_condition) +
#   scale_fill_manual(values = colour_set_level1) +
#   theme_classic()
# gg3 <- ggplot(pD) +
#   geom_violin(mapping = aes_string(x='level2', y='CD274', fill='level1')) +
#   # geom_beeswarm_rast(mapping = aes_string(x='level2', y='CD274', colour='tissue_condition')) +
#   scale_colour_manual(values = colour_set_tissue_condition) +
#   scale_fill_manual(values = colour_set_level1) +
#   theme_classic()
# gg <- rasterize(gg, layers='GeomPoint', dpi=500)
# gg2 <- rasterize(gg2, layers='GeomPoint', dpi=500)
# gg3 <- rasterize(gg3, layers='GeomPoint', dpi=500)
# pdf(paste0(fig_dir, 'violin/ggplot_PDL1_epi_level2.pdf'), width = 12, height = )
# gg
# dev.off()
# pdf(paste0(fig_dir, 'violin/ggplot_PDL1_epi_level2_log2.pdf'), width = 12, height = 7)
# gg2
# dev.off()
# pdf(paste0(fig_dir, 'violin/ggplot_PDL1_epi_level2_no_dots.pdf'), width = 12, height = 7)
# gg3
# dev.off()

#save directories
fig_dir <- '/nfs/research/marioni/areed/projects/hbca/figures/src/marker_plots/output/'
dir.create(paste0(fig_dir, 'violin/'), showWarnings = F, recursive = T)
dir.create(paste0(fig_dir, 'box/'), showWarnings = F, recursive = T)
dir.create(paste0(fig_dir, 'box2/'), showWarnings = F, recursive = T)
dir.create(paste0(fig_dir, 'box3/'), showWarnings = F, recursive = T)
dir.create(paste0(fig_dir, 'dot_age/'), showWarnings = F, recursive = T)

#plot violins
pD_bsl <- pD[(pD$level2 %in% c('BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]

pdf(paste0(fig_dir, 'violin/bsl_PDL1_tissue_condition.pdf'), width = 12, height = 8)
ggplot(pD_bsl) +
  geom_violin(mapping = aes_string(x='tissue_condition', y='CD274', fill='tissue_condition')) +
  geom_beeswarm_rast(mapping = aes_string(x='tissue_condition', y='CD274')) + #, colour='tissue_condition')) +
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic()
dev.off()

print('done violins')

#Violin not working as this is very lowly expressed do same as for milk genes

pD_bsl$ratio1 <- pD_bsl$CD274 > 0

metadata_sub <- pD_bsl[pD_bsl$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]

metadata_sum1 <- metadata_sub %>%
  group_by(patientID, level1, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum2 <- metadata_sub %>%
  group_by(patientID, level1, tissue_condition, patient_age, .drop=FALSE) %>%
  summarise(n=n(), avg_gene=mean(logcounts))
metadata_sum1.1 <- metadata_sum1 %>%
  group_by(tissue_condition, ratio1, .drop=FALSE) %>%
  summarise(mean=mean(freq), sd=sd(freq)) 

dir.create(paste0(fig_dir, 'barplot/'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(fig_dir, 'box/'), recursive = TRUE, showWarnings = FALSE)

pdf(paste0(fig_dir, 'barplot/', 'bsl_tissue_condition_PDL1.pdf'))
print(ggplot(data = metadata_sum1.1[metadata_sum1.1$ratio1,], aes(x=tissue_condition, y=mean, fill=tissue_condition)) + 
        geom_bar(stat="identity") +  
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.4, alpha=0.7, size=1) +
        scale_fill_manual(values = colour_set_tissue_condition) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()
pdf(paste0(fig_dir, 'barplot/', 'bsl_tissue_condition_PDL1_no_errbar.pdf'))
print(ggplot(data = metadata_sum1.1[metadata_sum1.1$ratio1,], aes(x=tissue_condition, y=mean, fill=tissue_condition)) + 
        geom_bar(stat="identity") +  
        #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.4, alpha=0.7, size=1) +
        scale_fill_manual(values = colour_set_tissue_condition) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
dev.off()

print('done bar plots')

#even better using box plots
dir.create(paste0(fig_dir, 'box/'), showWarnings = F, recursive = T)

gene_list <- c('PDCD1LG2', 'CD80', 'CD86', 'PVR', 'LGALS9', 'CD274')
for (gene in gene_list) {
  
  pD[, gene] <- logcounts(sce_sub)[gene,]
  pD$logcounts <- pD[, gene]
  pD <- as.data.frame(pD)
  pD$tissue_condition[pD$tissue_condition %in% c('Mastectomy WT', 'Mastectomy unknown')] <- 'HR-unk'
  pD$tissue_condition <- factor(pD$tissue_condition, levels=c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')) #, 'HR-unk'))
  
  ##BSL
  pD_bsl <- pD[(pD$level2 %in% c('BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  pD_bsl$ratio1 <- factor(pD_bsl[, gene] > 0)
  
  metadata_sub <- pD_bsl[pD_bsl$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1,tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box/', 'bsl_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box/', 'bsl_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=tissue_condition, y=avg_gene, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##LP
  pD_lp <- pD[(pD$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4', 'LP5')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_lp$ratio1 <- factor(pD_lp[, gene] > 0)
  
  metadata_sub <- pD_lp[pD_lp$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box/', 'lp_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box/', 'lp_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=tissue_condition, y=avg_gene, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##HS
  
  pD_hs <- pD[(pD$level2 %in% c('HS1', 'HS2', 'HS3', 'HS4')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_hs$ratio1 <- factor(pD_hs[, gene] > 0)
  
  metadata_sub <- pD_hs[pD_hs$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts)) 
  
  pdf(paste0(fig_dir, 'box/', 'hs_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box/', 'hs_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=tissue_condition, y=avg_gene, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(aes()) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
}

##Try accounting for/visualising celltype subclusters

gene_list <- c('PDCD1LG2', 'CD80', 'CD86', 'PVR', 'LGALS9', 'CD274')
for (gene in gene_list) {
  
  pD[, gene] <- logcounts(sce_sub)[gene,]
  pD$logcounts <- pD[, gene]
  pD <- as.data.frame(pD)
  
  ##Epi
  pD_epi <- pD[(pD$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4', 'HS1', 'HS2', 'HS3', 'HS4', 'BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  pD_epi$ratio1 <- factor(pD_epi[, gene] > 0)
  
  metadata_sub <- pD_epi[pD_epi$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  #check have greater than 20 cells of this subcluster for this patient
  metadata_sum1_check <- metadata_sum1 %>%
    group_by(patientID, level2) %>%
    summarise(n_total=sum(n))
  metadata_sum1_check$keep <- metadata_sum1_check$n_total > 20
  metadata_sum1_check2 <- dcast(metadata_sum1_check, patientID ~ level2, value.var='keep')
  metadata_sum1_check2[is.na(metadata_sum1_check2)] <- FALSE
  rownames(metadata_sum1_check2) <- metadata_sum1_check2$patientID
  metadata_sum1_check2$patientID <- NULL
  for (patientID_use in rownames(metadata_sum1_check2)) {
    for (level2_use in colnames(metadata_sum1_check2)) {
      #print(patientID_use)
      #print(level2_use)
      #print((metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use))
      if (is.na(metadata_sum1_check2[patientID_use, level2_use])){
        #print('Use datapoint.')
      } else if (metadata_sum1_check2[patientID_use, level2_use]) {
        #print('Use datapoint.')
      } else {
        #print('Do not use datapoint.')
        metadata_sum1[(metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use), c('n', 'freq')] <- NA 
        metadata_sum2[(metadata_sum2$patientID == patientID_use) & (metadata_sum2$level2 == level2_use), c('avg_gene')] <- NA 
      }
    }
  }
  
  pdf(paste0(fig_dir, 'box2/', 'epi_tissue_condition_', gene, '.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black', 
                                        'black', 'black', 'black', 'black', 
                                        'black', 'black')) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'epi_tissue_condition_', gene, '_2.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'epi_tissue_condition_', gene, '_logcounts.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##BSL
  pD_bsl <- pD[(pD$level2 %in% c('BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  pD_bsl$ratio1 <- factor(pD_bsl[, gene] > 0)
  
  metadata_sub <- pD_bsl[pD_bsl$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  #check have greater than 20 cells of this subcluster for this patient
  metadata_sum1_check <- metadata_sum1 %>%
    group_by(patientID, level2) %>%
    summarise(n_total=sum(n))
  metadata_sum1_check$keep <- metadata_sum1_check$n_total > 20
  metadata_sum1_check2 <- dcast(metadata_sum1_check, patientID ~ level2, value.var='keep')
  metadata_sum1_check3 <- dcast(metadata_sum1, patientID ~ level2, value.var='n', fun.aggregate=sum)
  metadata_sum1_check2[is.na(metadata_sum1_check2)] <- FALSE
  rownames(metadata_sum1_check2) <- metadata_sum1_check2$patientID
  metadata_sum1_check2$patientID <- NULL
  rownames(metadata_sum1_check3) <- metadata_sum1_check3$patientID
  metadata_sum1_check3$patientID <- NULL
  metadata_sum1_check3 <- metadata_sum1_check3 > 20   
  for (patientID_use in rownames(metadata_sum1_check2)) {
    for (level2_use in colnames(metadata_sum1_check2)) {
      #print(patientID_use)
      #print(level2_use)
      #print((metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use))
      if (is.na(metadata_sum1_check2[patientID_use, level2_use])){
        #print('Use datapoint.')
      } else if (metadata_sum1_check2[patientID_use, level2_use]) {
        #print('Use datapoint.')
      } else {
        #print('Do not use datapoint.')
        metadata_sum1[(metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use), c('n', 'freq')] <- NA 
        metadata_sum2[(metadata_sum2$patientID == patientID_use) & (metadata_sum2$level2 == level2_use), c('avg_gene')] <- NA 
      }
    }
  }
  
  pdf(paste0(fig_dir, 'box2/', 'bsl_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'bsl_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'bsl_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##LP
  pD_lp <- pD[(pD$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_lp$ratio1 <- factor(pD_lp[, gene] > 0)
  
  metadata_sub <- pD_lp[pD_lp$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  #check have greater than 20 cells of this subcluster for this patient
  metadata_sum1_check <- metadata_sum1 %>%
    group_by(patientID, level2) %>%
    summarise(n_total=sum(n))
  metadata_sum1_check$keep <- metadata_sum1_check$n_total > 20
  metadata_sum1_check2 <- dcast(metadata_sum1_check, patientID ~ level2, value.var='keep')
  metadata_sum1_check2[is.na(metadata_sum1_check2)] <- FALSE
  rownames(metadata_sum1_check2) <- metadata_sum1_check2$patientID
  metadata_sum1_check2$patientID <- NULL
  for (patientID_use in rownames(metadata_sum1_check2)) {
    for (level2_use in colnames(metadata_sum1_check2)) {
      #print(patientID_use)
      #print(level2_use)
      #print((metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use))
      if (is.na(metadata_sum1_check2[patientID_use, level2_use])){
        #print('Use datapoint.')
      } else if (metadata_sum1_check2[patientID_use, level2_use]) {
        #print('Use datapoint.')
      } else {
        #print('Do not use datapoint.')
        metadata_sum1[(metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use), c('n', 'freq')] <- NA 
        metadata_sum2[(metadata_sum2$patientID == patientID_use) & (metadata_sum2$level2 == level2_use), c('avg_gene')] <- NA 
      }
    }
  }
  
  pdf(paste0(fig_dir, 'box2/', 'lp_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'lp_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'lp_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##HS
  
  pD_hs <- pD[(pD$level2 %in% c('HS1', 'HS2', 'HS3', 'HS4')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_hs$ratio1 <- factor(pD_hs[, gene] > 0)
  
  metadata_sub <- pD_hs[pD_hs$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  #check have greater than 20 cells of this subcluster for this patient
  metadata_sum1_check <- metadata_sum1 %>%
    group_by(patientID, level2) %>%
    summarise(n_total=sum(n))
  metadata_sum1_check$keep <- metadata_sum1_check$n_total > 20
  metadata_sum1_check2 <- dcast(metadata_sum1_check, patientID ~ level2, value.var='keep')
  metadata_sum1_check2[is.na(metadata_sum1_check2)] <- FALSE
  rownames(metadata_sum1_check2) <- metadata_sum1_check2$patientID
  metadata_sum1_check2$patientID <- NULL
  for (patientID_use in rownames(metadata_sum1_check2)) {
    for (level2_use in colnames(metadata_sum1_check2)) {
      #print(patientID_use)
      #print(level2_use)
      #print((metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use))
      if (is.na(metadata_sum1_check2[patientID_use, level2_use])){
        #print('Use datapoint.')
      } else if (metadata_sum1_check2[patientID_use, level2_use]) {
        #print('Use datapoint.')
      } else {
        #print('Do not use datapoint.')
        metadata_sum1[(metadata_sum1$patientID == patientID_use) & (metadata_sum1$level2 == level2_use), c('n', 'freq')] <- NA 
        metadata_sum2[(metadata_sum2$patientID == patientID_use) & (metadata_sum2$level2 == level2_use), c('avg_gene')] <- NA 
      }
    }
  }
  
  pdf(paste0(fig_dir, 'box2/', 'hs_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'hs_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'hs_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}


#Make some plots without the cutoff
gene_list <- c('PDCD1LG2', 'CD80', 'CD86', 'PVR', 'LGALS9', 'CD274')
for (gene in gene_list) {
  
  pD[, gene] <- logcounts(sce_sub)[gene,]
  pD$logcounts <- pD[, gene]
  pD <- as.data.frame(pD)
  
  ##Epi
  pD_epi <- pD[(pD$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4', 'HS1', 'HS2', 'HS3', 'HS4', 'BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  pD_epi$ratio1 <- factor(pD_epi[, gene] > 0)
  
  metadata_sub <- pD_epi[pD_epi$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_epi_tissue_condition_', gene, '.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black', 
                                        'black', 'black', 'black', 'black', 
                                        'black', 'black')) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_epi_tissue_condition_', gene, '_2.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_epi_tissue_condition_', gene, '_logcounts.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_epi_tissue_condition_', gene, '_3.pdf'), height = 8, width = 16)
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##BSL
  pD_bsl <- pD[(pD$level2 %in% c('BSL1', 'BSL2')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  pD_bsl$ratio1 <- factor(pD_bsl[, gene] > 0)
  
  metadata_sub <- pD_bsl[pD_bsl$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_bsl_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_bsl_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_bsl_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_bsl_tissue_condition_', gene, '_3.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##LP
  pD_lp <- pD[(pD$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_lp$ratio1 <- factor(pD_lp[, gene] > 0)
  
  metadata_sub <- pD_lp[pD_lp$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_lp_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_lp_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_lp_tissue_condition_', gene, '_3.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_lp_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  #This doesn't show much at all
  pdf(paste0(fig_dir, 'dot_age/', 'nocutoff_lp_tissue_condition_', gene, '_3.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=patient_age, y=freq, color=tissue_condition, shape = level2)) + 
          geom_point() +
          scale_color_manual(values = colour_set_tissue_condition) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  
  ##HS
  
  pD_hs <- pD[(pD$level2 %in% c('HS1', 'HS2', 'HS3', 'HS4')) & (pD$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')), ]
  
  pD_hs$ratio1 <- factor(pD_hs[, gene] > 0)
  
  metadata_sub <- pD_hs[pD_hs$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
  
  metadata_sum1 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, ratio1, .drop=FALSE) %>%
    summarise(n=n()) %>%
    mutate(freq = n / sum(n))
  metadata_sum2 <- metadata_sub %>%
    group_by(patientID, level1, level2, tissue_condition, patient_age, .drop=FALSE) %>%
    summarise(n=n(), avg_gene=mean(logcounts))
  
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_hs_tissue_condition_', gene, '.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=tissue_condition, y=freq, color=level2, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black', 'black')) + #c('HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028")) +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_hs_tissue_condition_', gene, '_2.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_point(position=position_dodge(width=0.75)) +
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_hs_tissue_condition_', gene, '_3.pdf'))
  print(ggplot(data = metadata_sum1[metadata_sum1$ratio1==TRUE,], aes(x=level2, y=freq, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
  pdf(paste0(fig_dir, 'box2/', 'nocutoff_hs_tissue_condition_', gene, '_logcounts.pdf'))
  print(ggplot(data = metadata_sum2, aes(x=level2, y=avg_gene, color=tissue_condition, fill=tissue_condition)) + 
          geom_boxplot() +  
          geom_count(position=position_dodge(width=0.75)) +
          scale_size_area() + 
          scale_fill_manual(values = colour_set_tissue_condition) +
          scale_color_manual(values = c('black', 'black', 'black')) + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
}
















###CSN2 etc

sce_lp <- sce[, sce$level2 %in% c('LP1', 'LP2', 'LP3', 'LP4')]
pD_lp <- as.data.frame(colData(sce_lp))

genes0 = 'CSN2'
genes1 = c('CSN2', 'CSN1S1', 'LALBA')
genes2 = c(genes1, 'CSN3')

#make srt object
srt <- Seurat::CreateSeuratObject(counts = logcounts(sce_lp), min.cells = 3, project = "lp")
srt <- Seurat::AddMetaData(object = srt, metadata = pD_lp$sampleID, col.name = "sampleID")
srt <- Seurat::AddMetaData(object = srt, metadata = pD_lp$patientID, col.name = "patientID")
srt <- Seurat::AddMetaData(object = srt, metadata = pD_lp$tissue_condition, col.name = "tissue_condition")
srt <- Seurat::AddMetaData(object = srt, metadata = pD_lp$level2, col.name = "level2")

#Calculate ModuleScore's
srt <- Seurat::AddModuleScore(srt, features = list(genes1), ctrl = 20, name= "genes1")
srt <- Seurat::AddModuleScore(srt, features = list(genes2), ctrl = 20, name= "genes2")

#convert scores back into sce
sce_lp$CSN2 <- as.numeric(logcounts(sce_lp)['CSN2',])
sce_lp$CSN1S1 <- as.numeric(logcounts(sce_lp)['CSN1S1',])
sce_lp$CSN3 <- as.numeric(logcounts(sce_lp)['CSN3',])
sce_lp$LALBA <- as.numeric(logcounts(sce_lp)['LALBA',])
sce_lp$BCL11A <- as.numeric(logcounts(sce_lp)['BCL11A',])
sce_lp$genes11 <- srt@meta.data$genes11
sce_lp$genes21 <- srt@meta.data$genes21
pD_lp <- as.data.frame(colData(sce_lp))
pD_lp$parous <- 'nul'
pD_lp$parous[pD_lp$parity %in% c('1','2','3','4')] <- 'par'
pD_lp$parous[pD_lp$parity %in% c('Unknown')] <- 'unk'
pD_lp$tc_par <- paste0(pD_lp$tissue_condition, ' ', pD_lp$parous)

#plot as violin
for (gene_set in c('CSN2', 'genes11', 'genes21')) {
  gg_violin1 <- ggplot(pD_lp,
                mapping = aes_string(x='tc_par', y=gene_set, fill='tissue_condition')) +
    geom_violin() +
    #geom_beeswarm_rast() +
    #scale_colour_manual(values = colour_set_tissue_condition) +
    scale_fill_manual(values = colour_set_tissue_condition)
  gg_violin1 <- rasterize(gg_violin1, layers='GeomPoint', dpi=500)
  
  gg_box1 <- ggplot(pD_lp,
                       mapping = aes_string(x='tc_par', y=gene_set, 
                                            #colour='tissue_condition',
                                            fill = 'tissue_condition')) +
    geom_boxplot() +
    #geom_beeswarm_rast() +
    #scale_colour_manual(values = colour_set_tissue_condition) +
    scale_fill_manual(values = colour_set_tissue_condition)
  gg_box1 <- rasterize(gg_violin1, layers='GeomPoint', dpi=500)
  
  dir.create(paste0(fig_dir, 'violin/'), showWarnings = F, recursive = T)
  dir.create(paste0(fig_dir, 'box/'), showWarnings = F, recursive = T)
  
  pdf(paste0(fig_dir, 'violin/ggplot_', gene_set, '_lp_tc_par.pdf'), width = 12, height = 7)
  print(gg_violin1)
  dev.off()
  
  pdf(paste0(fig_dir, 'box/ggplot_', gene_set, '_lp_tc_par.pdf'), width = 12, height = 7)
  print(gg_box1)
  dev.off()
  }  


#These are not very instructive so lets look at ratios of cells with non-zero counts

pD_lp$ratio1 <- pD_lp$CSN2 > 0
pD_lp$ratio2 <- pD_lp$CSN1S1 > 0
pD_lp$ratio3 <- pD_lp$CSN3 > 0
pD_lp$ratio4 <- pD_lp$LALBA > 0
pD_lp$ratio5 <- pD_lp$BCL11A > 0
pD_lp$ratio_all <- pD_lp$ratio1 | pD_lp$ratio2 | pD_lp$ratio3 | pD_lp$ratio4

metadata_sub <- pD_lp[pD_lp$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2'), ]
metadata_sub$patientID <- factor(metadata_sub$patientID)
#metadata_sub$tissue_condition <- factor(metadata_sub$tissue_condition)
#metadata_sub$tc_par <- factor(metadata_sub$tc_par)
metadata_sub$ratio1 <- factor(metadata_sub$ratio1)
metadata_sub$ratio2 <- factor(metadata_sub$ratio2)
metadata_sub$ratio3 <- factor(metadata_sub$ratio3)
metadata_sub$ratio4 <- factor(metadata_sub$ratio4)
metadata_sub$ratio5 <- factor(metadata_sub$ratio5)
metadata_sub$ratio_all <- factor(metadata_sub$ratio_all)


metadata_sum1 <- metadata_sub %>%
  group_by(patientID, ratio1, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum1 <- left_join(metadata_sum1, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')

metadata_sum2 <- metadata_sub %>%
  group_by(patientID, ratio2, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum2 <- left_join(metadata_sum2, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')

metadata_sum3 <- metadata_sub %>%
  group_by(patientID, ratio3, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum3 <- left_join(metadata_sum3, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')

metadata_sum4 <- metadata_sub %>%
  group_by(patientID, ratio4, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum4 <- left_join(metadata_sum4, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')

metadata_sum5 <- metadata_sub %>%
  group_by(patientID, ratio5, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum5 <- left_join(metadata_sum5, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')

metadata_sum_all <- metadata_sub %>%
  group_by(patientID, ratio_all, .drop=F) %>% #, .drop = c(T, T, F)) %>%
  summarise(n=n()) %>%
  mutate(freq = n / sum(n))
metadata_sum_all <- left_join(metadata_sum_all, unique(metadata_sub[, c('patientID', 'tc_par', 'tissue_condition')]), by='patientID')



#patientID columns with celltype fill
pdf(paste0(fig_dir, 'box/', 'lp_tc_par_milk_gene.pdf'))
ggplot(data = metadata_sum_all[metadata_sum_all$ratio_all == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN2.pdf'))
ggplot(data = metadata_sum1[metadata_sum1$ratio1 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  scale_fill_manual(values = colour_set_tissue_condition) +
  scale_y_break(c(0.0035,0.005), scales=0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN2_no_cut.pdf'))
ggplot(data = metadata_sum1[metadata_sum1$ratio1 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  #scale_y_break(c(0.0035,0.005), scales=0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN1S1.pdf'))
ggplot(data = metadata_sum2[metadata_sum2$ratio2 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  
  scale_fill_manual(values = colour_set_tissue_condition) +
  scale_y_break(c(0.04,0.06), scales=0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN1S1_no_cut.pdf'))
ggplot(data = metadata_sum2[metadata_sum2$ratio2 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN3.pdf'))
ggplot(data = metadata_sum3[metadata_sum3$ratio3 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) + 
  scale_fill_manual(values = colour_set_tissue_condition) +
  scale_y_break(c(0.03,0.08), scales=0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
pdf(paste0(fig_dir, 'box/', 'lp_tc_par_CSN3_no_cut.pdf'))
ggplot(data = metadata_sum3[metadata_sum3$ratio3 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) + 
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(fig_dir, 'box/', 'lp_tc_par_LALBA_no_cut.pdf'))
ggplot(data = metadata_sum4[metadata_sum4$ratio4 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(fig_dir, 'box/', 'lp_tc_par_BCL11A.pdf'))
ggplot(data = metadata_sum5[metadata_sum5$ratio5 == 'TRUE',], aes(x=tc_par, y=freq, fill=tissue_condition)) + 
  geom_boxplot() + 
  geom_point(aes()) +  
  scale_fill_manual(values = colour_set_tissue_condition) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

