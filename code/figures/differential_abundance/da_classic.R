#Regular differential abundance testing - AR vs HR-unk


####A few issues here 

### I haven't removed the facs sorted sample types
### There is no kind of normalisation per sample - not sure if the model can account for this...

#librarys
library(plyr)
library(dplyr)
library(scran)
library(edgeR)
library(scater)
library(ggplot2)
library(EnhancedVolcano)

#load data
metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')

metadata_sub <- metadata_all[,c('sampleID', 'patientID', 'before', 'parity', 'patient_age', 'tissue_condition', 'level2')]
metadata_sub_test <- metadata_sub[(metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT', 'Mastectomy unknown')) & 
                                    metadata_sub$level2 != 'Doublet' & 
                                    metadata_sub$before %in% c('Organoid unsorted', 'Supernatant unsorted')
                                  , ] ####previously did not account for this!
metadata_sub_test$test <- 'AR'
metadata_sub_test$test[metadata_sub_test$tissue_condition != 'Mammoplasty WT'] <- 'HR-unk'
metadata_sub_test$parous <- 'Nulliparous'
metadata_sub_test$parous[metadata_sub_test$parity != '0'] <- 'Parous'

abundances <- table(metadata_sub_test$level2, metadata_sub_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

#colours
level2_colour_dictionary = c('LP1' = "#DA80DA", 'LP2' = "#815481", 'LP3' = "#C040C0", 'LP4' = "#E1AFE1", 
                             'HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028",
                             'BSL1' = "#EB675E", 'BSL2' = "#A23E36",
                             'FB1' = "#DFA38A", 'FB2' = "#8C3612", 'FB3' = "#623623", 'FB4' = "#916350", 'FB5' = "#DAC3C3",
                             'VM1' = "#F8770B", 'VM2' = "#E09E3A", 'VM3' = "#CD7225", 'VM4' = "#FFC990", 'VM5' = "#AC5812",
                             'EC venous' = "#FEE083", 'EC capillary' = "#897538", 'EC arterial' = "#E7B419", 'EC angiogenic tip' = "#BCA048",
                             'LEC1' = "#6F8BE2", 'LEC2' = "#3053BC",
                             'CD8T 1' = "#6D9F58", 'CD8T 2' = "#9EB766", 'CD8T 3' = "#BDCB10", 'CD4T' = "#3A6527", 'IFNG+ T' = "#9EA743",
                             'NK1' = "#E2E8A7", 'NK2' = "#5A6209", 'NK3' = "#8FE36B",
                             'ILC3' = "#818A31",
                             'B cell' = "#9FC5E8", 'Plasma cell' = "#23D9F1",
                             'Macrophage' = "#64C6A6") 


extra.info <- metadata_sub_test[match(colnames(abundances), metadata_sub_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$test)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

#Do not normalise like this for DA testing.
#Normalize if not already done
# y.ab <- calcNormFactors(y.ab)
# y.ab

design <- model.matrix(~ before + parous + patient_age + test, y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

#plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

res$table$FDR <- p.adjust(res$table$PValue, method="BH")
DAIs <- res$table[order(res$table$FDR),]
dim(DAIs) #41x5
save_path = '/nfs/research/marioni/areed/projects/hbca/figures/src/differential_abundance/output/classical/'
dir.create(save_path, showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'DA_level2.csv'))

Volcano <- EnhancedVolcano(DAIs,
                           lab = rownames(DAIs),
                           x = 'logFC',
                           y = 'FDR',
                           title = 'DA for cell types AR vs HR-unk',
                           subtitle = '',
                           subtitleLabSize = 2,
                           legendPosition = "bottom",
                           pointSize = 6.0,
                           labSize = 2.0,
                           FCcutoff = 1,
                           pCutoff = 0.05,
                           xlim=c(-5,7),
                           ylim=c(0,4),
                           col = c("grey", "forestgreen", "steelblue", "red"),
                           #legendVisible = FALSE,
                           drawConnectors = TRUE,
                           typeConnectors = 'open',
                           max.overlaps=Inf)
DAIs$celltype <- rownames(DAIs)
volcano2 <- ggplot(DAIs, aes(x=logFC, y=-log10(FDR),fill=celltype)) +
  geom_point(pch=21, size=6, alpha=0.8) +
  geom_hline(yintercept=-log10(0.05),size=1.1, lty="dashed") +
  geom_text_repel(data=DAIs[DAIs$FDR<0.1,], label=DAIs[DAIs$FDR<0.1,"celltype"],aes(fill=NULL),
                  size=6, color="grey30",force=9) +
  scale_fill_manual(values=level2_colour_dictionary) +
  #theme_pub() +
  theme(legend.position="none") +
  xlab("Log2(FC)")

pdf(paste0(save_path, 'DA_volcano_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
pdf(paste0(save_path, 'DA_volcano_level2_2.pdf'), width=12,height = 8)
volcano2
dev.off()



#####try splitting supernatant and organoids

#Organoid

metadata_sub <- metadata_all[,c('sampleID', 'patientID', 'before', 'parity', 'patient_age', 'tissue_condition', 'level2')]
metadata_sub_test <- metadata_sub[(metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT', 'Mastectomy unknown')) & 
                                    metadata_sub$level2 != 'Doublet' & 
                                    metadata_sub$before %in% c('Organoid unsorted')
                                  , ] ####previously did not account for this!
metadata_sub_test$test <- 'AR'
metadata_sub_test$test[metadata_sub_test$tissue_condition != 'Mammoplasty WT'] <- 'HR-unk'
metadata_sub_test$parous <- 'Nulliparous'
metadata_sub_test$parous[metadata_sub_test$parity != '0'] <- 'Parous'

abundances <- table(metadata_sub_test$level2, metadata_sub_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)


extra.info <- metadata_sub_test[match(colnames(abundances), metadata_sub_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$test)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

#Do not normalise like this for DA testing.
#Normalize if not already done
# y.ab <- calcNormFactors(y.ab)
# y.ab

design <- model.matrix(~ parous + patient_age + test, y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

#plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

res$table$FDR <- p.adjust(res$table$PValue, method="BH")
DAIs <- res$table[order(res$table$FDR),]
dim(DAIs) #41x5
save_path = '/nfs/research/marioni/areed/projects/hbca/figures/src/differential_abundance/output/classical/'
dir.create(save_path, showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'DA_level2_organoid.csv'))

Volcano <- EnhancedVolcano(DAIs,
                           lab = rownames(DAIs),
                           x = 'logFC',
                           y = 'FDR',
                           title = 'DA for cell types AR vs HR-unk',
                           subtitle = '',
                           subtitleLabSize = 2,
                           legendPosition = "bottom",
                           pointSize = 6.0,
                           labSize = 2.0,
                           FCcutoff = 1,
                           pCutoff = 0.05,
                           xlim=c(-5,7),
                           ylim=c(0,4),
                           col = c("grey", "forestgreen", "steelblue", "red"),
                           #legendVisible = FALSE,
                           drawConnectors = TRUE,
                           typeConnectors = 'open',
                           max.overlaps=Inf)

pdf(paste0(save_path, 'DA_volcano_level2_organoid.pdf'), width=12,height = 8)
Volcano
dev.off()

#Supernatant

metadata_sub <- metadata_all[,c('sampleID', 'patientID', 'before', 'parity', 'patient_age', 'tissue_condition', 'level2')]
metadata_sub_test <- metadata_sub[(metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT', 'Mastectomy unknown')) & 
                                    metadata_sub$level2 != 'Doublet' & 
                                    metadata_sub$before %in% c('Supernatant unsorted')
                                  , ] ####previously did not account for this!
metadata_sub_test$test <- 'AR'
metadata_sub_test$test[metadata_sub_test$tissue_condition != 'Mammoplasty WT'] <- 'HR-unk'
metadata_sub_test$parous <- 'Nulliparous'
metadata_sub_test$parous[metadata_sub_test$parity != '0'] <- 'Parous'

abundances <- table(metadata_sub_test$level2, metadata_sub_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)


extra.info <- metadata_sub_test[match(colnames(abundances), metadata_sub_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$test)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

#Do not normalise like this for DA testing.
#Normalize if not already done
# y.ab <- calcNormFactors(y.ab)
# y.ab

design <- model.matrix(~ parous + patient_age + test, y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

#plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

res$table$FDR <- p.adjust(res$table$PValue, method="BH")
DAIs <- res$table[order(res$table$FDR),]
dim(DAIs) #41x5
save_path = '/nfs/research/marioni/areed/projects/hbca/figures/src/differential_abundance/output/classical/'
dir.create(save_path, showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'DA_level2_supernatant.csv'))

Volcano <- EnhancedVolcano(DAIs,
                           lab = rownames(DAIs),
                           x = 'logFC',
                           y = 'FDR',
                           title = 'DA for cell types AR vs HR-unk',
                           subtitle = '',
                           subtitleLabSize = 2,
                           legendPosition = "bottom",
                           pointSize = 6.0,
                           labSize = 2.0,
                           FCcutoff = 1,
                           pCutoff = 0.05,
                           xlim=c(-5,7),
                           ylim=c(0,4),
                           col = c("grey", "forestgreen", "steelblue", "red"),
                           #legendVisible = FALSE,
                           drawConnectors = TRUE,
                           typeConnectors = 'open',
                           max.overlaps=Inf)

pdf(paste0(save_path, 'DA_volcano_level2_supernatant.pdf'), width=12,height = 8)
Volcano
dev.off()


#### No age parity blocking

metadata_sub <- metadata_all[,c('sampleID', 'patientID', 'before', 'parity', 'patient_age', 'tissue_condition', 'level2')]
metadata_sub_test <- metadata_sub[(metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT', 'Mastectomy unknown')) & 
                                    metadata_sub$level2 != 'Doublet' & 
                                    metadata_sub$before %in% c('Organoid unsorted', 'Supernatant unsorted')
                                  , ] ####previously did not account for this!
metadata_sub_test$test <- 'AR'
metadata_sub_test$test[metadata_sub_test$tissue_condition != 'Mammoplasty WT'] <- 'HR-unk'
metadata_sub_test$parous <- 'Nulliparous'
metadata_sub_test$parous[metadata_sub_test$parity != '0'] <- 'Parous'

abundances <- table(metadata_sub_test$level2, metadata_sub_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)


extra.info <- metadata_sub_test[match(colnames(abundances), metadata_sub_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$test)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

#Do not normalise like this for DA testing.
#Normalize if not already done
# y.ab <- calcNormFactors(y.ab)
# y.ab

design <- model.matrix(~ before + test, y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)

#plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

res$table$FDR <- p.adjust(res$table$PValue, method="BH")
DAIs <- res$table[order(res$table$FDR),]
dim(DAIs) #41x5
save_path = '/nfs/research/marioni/areed/projects/hbca/figures/src/differential_abundance/output/classical/'
dir.create(save_path, showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'DA_level2_no_block.csv'))

Volcano <- EnhancedVolcano(DAIs,
                           lab = rownames(DAIs),
                           x = 'logFC',
                           y = 'FDR',
                           title = 'DA for cell types AR vs HR-unk',
                           subtitle = '',
                           subtitleLabSize = 2,
                           legendPosition = "bottom",
                           pointSize = 6.0,
                           labSize = 2.0,
                           FCcutoff = 1,
                           pCutoff = 0.05,
                           xlim=c(-5,7),
                           ylim=c(0,4),
                           col = c("grey", "forestgreen", "steelblue", "red"),
                           #legendVisible = FALSE,
                           drawConnectors = TRUE,
                           typeConnectors = 'open',
                           max.overlaps=Inf)

pdf(paste0(save_path, 'DA_volcano_level2_no_block.pdf'), width=12,height = 8)
Volcano
dev.off()
