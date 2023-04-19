#Regular differential abundance testing - AR vs HR-unk
#This time try the 'more robust' rlm testing (this is what Karsten actually used).
#based off karstens code

#librarys
library(plyr)
library(dplyr)
library(scran)
library(edgeR)
library(scater)
library(ggplot2)
library(EnhancedVolcano)
library(MASS)
library(sfsmisc)


#load data
metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_sub <- metadata_all[,c('sampleID', 'patientID', 'before', 'parity', 'patient_age', 'tissue_condition', 'level2')]
metadata_sub_test <- metadata_sub[(metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT', 'Mastectomy unknown')) & 
                                    metadata_sub$level2 != 'Doublet' & 
                                    metadata_sub$before %in% c('Organoid unsorted', 'Supernatant unsorted'), ]
metadata_sub_test$test <- 'AR'
metadata_sub_test$test[metadata_sub_test$tissue_condition != 'Mammoplasty WT'] <- 'HR-unk'
metadata_sub_test$parous <- 'Nulliparous'
metadata_sub_test$parous[metadata_sub_test$parity != '0'] <- 'Parous'

abundances <- table(metadata_sub_test$level2, metadata_sub_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

#set up DGE
extra.info <- metadata_sub_test[match(colnames(abundances), metadata_sub_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info)

tab <- data.frame()
for (cluster_use in rownames(abundances)) {
  y.ab$samples$Counts <- cpm(y.ab, log=TRUE)[cluster_use,]
  fit <- rlm(formula=Counts ~  before + parous + patient_age + test, data=y.ab$samples, maxit=100)
  ftest <- f.robftest(fit,var="testHR-unk")
  tmp <- data.frame("CellTypes"=cluster_use,
                    "lfc"=fit$coefficients[length(fit$coefficient)],
                    "intercept"=fit$coefficients[1],
                    "PVal"=ftest$p.value)
  tab <- rbind(tab,tmp)
}
tab$FDR <- as.numeric(p.adjust(tab$PVal, method="BH"))
tab[order(tab$FDR,decreasing=FALSE),]

##save results

save_path = '/nfs/research/marioni/areed/projects/hbca/figures/src/differential_abundance/output/classical_robust/'
dir.create(save_path, showWarnings = F, recursive = T)
write.csv(tab, file= paste0(save_path, 'DA_level2.csv'))

Volcano <- EnhancedVolcano(tab,
                           lab = rownames(tab),
                           x = 'lfc',
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

pdf(paste0(save_path, 'DA_volcano_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
