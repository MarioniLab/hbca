#Make patient bar plots over the celltypes.

library(plyr)
library(dplyr)
library(ggplot2)




### Load data

metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_sub <- metadata_all[,c('patientID', 'level1', 'level2')]

### Colours choices 13
colour_palette <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                    "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

tissue_condition_colour_dictionary = c('Mammoplasty WT' = '#f1ce63', 'Mastectomy BRCA1' = '#4e79a7',
                                       'Mastectomy BRCA2' = '#a0cbe8', 'Mastectomy WT' = '#f28e2b',
                                       'Mastectomy unknown' = '#bab0ac', 'Contralateral BRCA1' =  '#ff9d9a')

level1_colour_dictionary = c('Luminal progenitor' = '#DDA0DD', 'Luminal hormone sensing' = '#EE3A8C', 'Basal' = '#FF6347', 
                             'Fibroblast' = '#804E1A', 'Vascular mural' = '#F89440', 'Endothelial' = '#FFF08D',
                             "Lymphoid" = '#9FC5E8', "Myeloid" = '#AAB256') #ignore doublets
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

#make colour mapping for patientIDs to tissue_condition colours.
patientID_tissue_condition <- unique(data.frame(patientID = metadata_all$patientID, 'tissue_condition' = metadata_all$tissue_condition))
patientID_tissue_condition$tissue_condition <- factor(patientID_tissue_condition$tissue_condition, levels = names(tissue_condition_colour_dictionary))
patientID_tissue_condition <- patientID_tissue_condition[order(patientID_tissue_condition$tissue_condition, patientID_tissue_condition$patientID), ]
patientID_tissue_condition$colour <- tissue_condition_colour_dictionary[patientID_tissue_condition$tissue_condition]
patientID_tissue_condition_colours <- setNames(as.character(patientID_tissue_condition$colour), patientID_tissue_condition$patientID)

### Preprocessing

#remove Doublet cells
metadata_sub <- metadata_sub[metadata_sub$level2 != 'Doublet', ]

#aggregate
metadata_sum1 <- metadata_sub %>%
  group_by(patientID, level1) %>%
  summarise(n=n())

metadata_sum2 <- metadata_sub %>%
  group_by(patientID, level1, level2) %>%
  summarise(n=n())

celltype_list1 <- c("Luminal progenitor", "Luminal hormone sensing", "Basal",
                    "Fibroblast", "Endothelial", "Vascular mural",
                    "Lymphoid", "Myeloid")
celltype_list2 <- c("LP1", "LP2", "LP3", "LP4", "LP5",
                    "HS1", "HS2", "HS3", "HS4", 
                    "BSL1", "BSL2",
                    "DDC1", "DDC2",
                    "FB1", "FB2", "FB3", "FB4", "FB5",
                    "VM1", "VM2", "VM3", "VM4", "VM5",
                    "EC venous", "EC capillary", "EC arterial", 
                    "EC angiogenic tip", "LEC1", "LEC2",
                    'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                    'IFNG+ T', 'NK1', 'NK2', 'NK3',
                    'ILC3', 'B cell', 'Plasma cell', 'Macrophage')
metadata_sum1$level1 <- factor(metadata_sum1$level1,
                               levels = celltype_list1)
metadata_sum2$level2 <- factor(metadata_sum2$level2,
                               levels = rev(celltype_list2))
metadata_sum1$patientID <- factor(metadata_sum1$patientID,
                               levels = rev(patientID_tissue_condition$patientID)) #rev
metadata_sum2$patientID <- factor(metadata_sum2$patientID,
                               levels = (patientID_tissue_condition$patientID)) #rev
### Plots

fig_dir <- '/nfs/research/marioni/areed/projects/hbca/figures/src/patientcelltype_barplot/barplot/'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#patientID columns with celltype fill
pdf(paste0(fig_dir, 'patientID_x_level1_fill.pdf'))
ggplot(data = metadata_sum1, aes(x=patientID, y=n, fill=level1)) + 
  geom_bar(position='fill', stat='identity')  + 
  scale_fill_manual(values = level1_colour_dictionary) +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = rev(patientID_tissue_condition_colours))) #rev
dev.off()

pdf(paste0(fig_dir, 'patientID_x_level2_fill.pdf'), width = 20, height=6)
ggplot(data = metadata_sum2, aes(x=patientID, y=n, fill=level2)) + 
  geom_bar(position='fill', stat='identity')  + 
  scale_fill_manual(values = level2_colour_dictionary) +
  #coord_flip() + 
  theme_classic() +
  theme(axis.text.x = element_text(colour = (patientID_tissue_condition_colours),
                                   angle = 90, vjust = 0.5, hjust=1)) #rev
dev.off()

#Celltype with patientID fill.

#fix patientID order
metadata_sum2$patientID <- factor(metadata_sum2$patientID,
                                  levels = patientID_tissue_condition$patientID)

#Epithelia
pdf(paste0(fig_dir, 'level2_x_patientID_fill_epi.pdf'), width = 15, height = 7)
ggplot(data = metadata_sum2[metadata_sum2$level1 %in% c('Luminal progenitor', 'Luminal hormone sensing', 'Basal'), ], aes(x=level2, y=n, fill=patientID)) + 
  geom_bar(position='fill', stat='identity') + 
  scale_fill_manual(values = colour_palette[1:length(unique(metadata_sum2$patientID))]) +
  theme_classic() +
  coord_flip() 
dev.off()

#Stroma
pdf(paste0(fig_dir, 'level2_x_patientID_fill_str.pdf'), width = 15, height = 7)
ggplot(data = metadata_sum2[metadata_sum2$level1 %in% c("Fibroblast", "Endothelial", "Vascular mural"), ], aes(x=level2, y=n, fill=patientID)) + 
  geom_bar(position='fill', stat='identity') + 
  scale_fill_manual(values = colour_palette[1:length(unique(metadata_sum2$patientID))]) +
  theme_classic() +
  coord_flip()
dev.off()

#Immune
pdf(paste0(fig_dir, 'level2_x_patientID_fill_imm.pdf'), width = 15, height = 7)
ggplot(data = metadata_sum2[metadata_sum2$level1 %in% c("Lymphoid", "Myeloid"), ], aes(x=level2, y=n, fill=patientID)) + 
  geom_bar(position='fill', stat='identity') + 
  scale_fill_manual(values = colour_palette[1:length(unique(metadata_sum2$patientID))]) +
  theme_classic() +
  coord_flip()
dev.off()


















