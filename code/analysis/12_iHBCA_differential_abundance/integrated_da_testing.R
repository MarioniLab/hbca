#Complete da testing on ageing and parity of integrated data


#librarys
library(plyr)
library(dplyr)
library(scran)
library(edgeR)
library(scater)
library(ggplot2)
library(EnhancedVolcano)



metadata <- read.csv('/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/colData/integrated_datasets_with_extra_meta.csv')
metadata$sampleID <- paste0(metadata$patientID, '_', metadata$sample_type)

metadata_sub <- metadata[, c('sampleID', 'patientID', 'dataset', 'risk_status', 'live_sorted_boolean', 'sample_type', 'parity', 'patient_age', 'map_celltype_map2reed')]

#colours
level2_colour_dictionary = c('LASP1' = "#DA80DA", 'LASP2' = "#815481", 'LASP3' = "#C040C0", 'LASP4' = "#E1AFE1", 'LASP5' = "#3F0034",
                             'LHS1' = "#EDABB9", 'LHS2' = "#EB5C79", 'LHS3' = "#A06A75",
                             'BMYO1' = "#EB675E", 'BMYO2' = "#A23E36",
                             'FB1' = "#DFA38A", 'FB2' = "#8C3612", 'FB3' = "#623623", 'FB4' = "#916350", 'FB5' = "#DAC3C3",
                             'PV1' = "#F8770B", 'PV2' = "#E09E3A", 'PV3' = "#CD7225", 'PV4' = "#FFC990", 'PV5' = "#AC5812",
                             'VEV' = "#FEE083", 'VEC' = "#897538", 'VEA' = "#E7B419", 'VEAT' = "#BCA048",
                             'LE1' = "#6F8BE2", 'LE2' = "#3053BC",
                             'CD4_naive' = "#5D9047", 'CD4_Th' = "#3A6527",
                             'CD8_Tem' = "#8EDD6D", 'CD8_Trm' = "#9EB766",
                             'CD8_Tc1' = "#9EA743",
                             'NKT' = "#E2E8A7", 'NK' = "#5A6209",
                             'ILC' = "#818A31",
                             'B_naive' = "#1A64AA", 'B_mem_unswitched' = "#158AF6", 'B_mem_switched' = "#9FC5E8",
                             'Plasma_cell' = "#23D9F1",
                             'Macro' = "#64C6A6", 'Macro-lipo' = '#1F9D74',
                             'DC' = '#023E69') 

metadata_sub$map_celltype_map2reed <- gsub('HBCA-', '', metadata_sub$map_celltype_map2reed)
metadata_sub$map_celltype_map2reed <- gsub('True-', '', metadata_sub$map_celltype_map2reed)

metadata_sub <- metadata_sub[!(metadata_sub$map_celltype_map2reed %in% c('', 'Doublet', 'stripped_nuclei', 'DDC1', 'DDC2')),]
  
#Clean up for samples we can test
age_given_integer <- !(metadata_sub$patient_age %in% c('O','Y','Unknown', 'n/a')) & !is.na(metadata_sub$patient_age)
parity_known <- !(metadata_sub$parity %in% c('unknown', 'Unknown', 'n/a')) & !is.na(metadata_sub$parity)
not_celltype_sorted <- !(metadata_sub$live_sorted_boolean %in% c('remove', 'Unknown')) 
not_bad_sampletype <- !(metadata_sub$sample_type %in% c('Organoid LP', 'Unknown'))
AR_samples <- metadata_sub$risk_status %in% 'AR'
not_milk_samples <- !(metadata_sub$patientID %in% c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 'HMC4', 'HMC5', 
                                                    'HMC6', 'HMC7', 'HMC8', 'HMC9'))

good_cells <- age_given_integer & parity_known & not_celltype_sorted & not_bad_sampletype & AR_samples & not_milk_samples
reed_dataset_cells <- metadata_sub$dataset=='HBCA'
tbl <- table(good_cells, reed_dataset_cells)
print('Cells with required metadata (and AR):')
print(tbl)
print('Percentage from Reed dataset:')
print(tbl[2,2] / (tbl[2,1] + tbl[2,2])) #....


metadata_test <- metadata_sub[good_cells, ]
metadata_test$patient_age <- strtoi(metadata_test$patient_age)
metadata_test$parous <- 'Nulliparous'
metadata_test$parous[metadata_test$parity != '0'] <- 'Parous'

#issue with design matrix rank. HBCA dataset column is combination of sample_type columns! So simplify and do not block for sample_type
metadata_test$dataset <- paste0(metadata_test$dataset, '_', metadata_test$sample_type)

abundances <- table(metadata_test$map_celltype_map2reed, metadata_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)



#TEST AGEING IMPACT:

# Need to discretise ages to include Kumar data (Old/Young). Otherwise we would have to exclude Kumar data.

#Clean up for samples we can test
age_given <- !(metadata_sub$patient_age %in% c('Unknown', 'n/a')) & !is.na(metadata_sub$patient_age)
parity_known <- !(metadata_sub$parity %in% c('unknown', 'Unknown', 'n/a')) & !is.na(metadata_sub$parity)
not_celltype_sorted <- !(metadata_sub$live_sorted_boolean %in% c('remove', 'Unknown')) 
not_bad_sampletype <- !(metadata_sub$sample_type %in% c('Organoid LP', 'Unknown'))
AR_samples <- metadata_sub$risk_status %in% 'AR'
not_milk_samples <- !(metadata_sub$patientID %in% c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 'HMC4', 'HMC5', 
                                                    'HMC6', 'HMC7', 'HMC8', 'HMC9'))

good_cells <- age_given & parity_known & not_celltype_sorted & not_bad_sampletype & AR_samples & not_milk_samples
reed_dataset_cells <- metadata_sub$dataset=='HBCA'
tbl <- table(good_cells, reed_dataset_cells)
print('Cells with required metadata (and AR):')
print(tbl)
print('Percentage from Reed dataset:')
print(tbl[2,2] / (tbl[2,1] + tbl[2,2])) #0.2904089

metadata_test <- metadata_sub[good_cells, ]
metadata_test$patient_age_est <- mapvalues(metadata_test$patient_age, 
                                           from=c('O', 'Y'),
                                           to=c('51', '49'))
metadata_test$old_young <- strtoi(metadata_test$patient_age_est) >= 50  #TRUE is equivalent to Old (O), False = Young (Y)
metadata_test$parous <- 'Nulliparous'
metadata_test$parous[metadata_test$parity != '0'] <- 'Parous'

#issue with design matrix rank. HBCA dataset column is combination of sample_type columns! So simplify and do not block for sample_type
metadata_test$dataset <- paste0(metadata_test$dataset, '_', metadata_test$sample_type)

abundances <- table(metadata_test$map_celltype_map2reed, metadata_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

extra.info <- metadata_test[match(colnames(abundances), metadata_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$old_young)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

design <- model.matrix(~ dataset + live_sorted_boolean + parous + old_young, y.ab$samples)
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

save_path = '/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/'
dir.create(paste0(save_path, 'raw/'), showWarnings = F, recursive = T)
dir.create(paste0(save_path, 'volcano/'), showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'raw/DA_age_discrete_level2.csv'))

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
volcano2 <- ggplot(DAIs, aes(x=logFC, y=-log10(FDR), fill=celltype)) +
  geom_point(pch=21, size=6, alpha=0.8) +
  geom_hline(yintercept=-log10(0.05),size=1.1, lty="dashed") +
  geom_text_repel(data=DAIs[DAIs$FDR<0.1,], label=DAIs[DAIs$FDR<0.1,"celltype"],aes(fill=NULL),
                  size=6, color="grey30",force=9) +
  scale_fill_manual(values=level2_colour_dictionary) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Log2(FC)")

pdf(paste0(save_path, 'volcano/DA_age_discrete_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
pdf(paste0(save_path, 'volcano/DA_age_discrete_level2_2.pdf'), width=12,height = 8)
volcano2
dev.off()





#TEST PARITY IMPACT:

#much the same setup just different design matrix
abundances <- table(metadata_test$map_celltype_map2reed, metadata_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

extra.info <- metadata_test[match(colnames(abundances), metadata_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$old_young)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

design <- model.matrix(~ dataset + live_sorted_boolean + old_young + parous, y.ab$samples)
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

save_path = '/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/'
dir.create(paste0(save_path, 'raw/'), showWarnings = F, recursive = T)
dir.create(paste0(save_path, 'volcano/'), showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'raw/DA_parity_b_age_discrete_level2.csv'))

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
volcano2 <- ggplot(DAIs, aes(x=logFC, y=-log10(FDR), fill=celltype)) +
  geom_point(pch=21, size=6, alpha=0.8) +
  geom_hline(yintercept=-log10(0.05),size=1.1, lty="dashed") +
  geom_text_repel(data=DAIs[DAIs$FDR<0.1,], label=DAIs[DAIs$FDR<0.1,"celltype"],aes(fill=NULL),
                  size=6, color="grey30",force=9) +
  scale_fill_manual(values=level2_colour_dictionary) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Log2(FC)")

pdf(paste0(save_path, 'volcano/DA_parity_b_age_discrete_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
pdf(paste0(save_path, 'volcano/DA_parity_b_age_discrete_level2_2.pdf'), width=12,height = 8)
volcano2
dev.off()




#TEST BRCA1 IMPACT:

#Clean up for samples we can test
age_given <- !(metadata_sub$patient_age %in% c('Unknown', 'n/a')) & !is.na(metadata_sub$patient_age)
parity_known <- !(metadata_sub$parity %in% c('unknown', 'Unknown', 'n/a')) & !is.na(metadata_sub$parity)
not_celltype_sorted <- !(metadata_sub$live_sorted_boolean %in% c('remove', 'Unknown')) 
not_bad_sampletype <- !(metadata_sub$sample_type %in% c('Organoid LP', 'Unknown'))
AR_BR1_samples <- metadata_sub$risk_status %in% c('AR', 'HR-BR1')
not_milk_samples <- !(metadata_sub$patientID %in% c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 'HMC4', 'HMC5', 
                                                    'HMC6', 'HMC7', 'HMC8', 'HMC9'))

good_cells <- age_given & parity_known & not_celltype_sorted & not_bad_sampletype & AR_BR1_samples & not_milk_samples
reed_dataset_cells <- metadata_sub$dataset=='HBCA'
tbl <- table(good_cells, reed_dataset_cells)
print('Cells with required metadata (and AR):')
print(tbl)
print('Percentage from Reed dataset:')
print(tbl[2,2] / (tbl[2,1] + tbl[2,2])) #0.3904928

metadata_test <- metadata_sub[good_cells, ]
metadata_test$patient_age_est <- mapvalues(metadata_test$patient_age, 
                                           from=c('O', 'Y'),
                                           to=c('51', '49'))
metadata_test$old_young <- strtoi(metadata_test$patient_age_est) >= 50  #TRUE is equivalent to Old (O), False = Young (Y)
metadata_test$parous <- 'Nulliparous'
metadata_test$parous[metadata_test$parity != '0'] <- 'Parous'

#issue with design matrix rank. HBCA dataset column is combination of sample_type columns! So simplify and do not block for sample_type
metadata_test$dataset <- paste0(metadata_test$dataset, '_', metadata_test$sample_type)

abundances <- table(metadata_test$map_celltype_map2reed, metadata_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

extra.info <- metadata_test[match(colnames(abundances), metadata_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$old_young)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

design <- model.matrix(~ dataset + live_sorted_boolean + parous + old_young + risk_status, y.ab$samples)
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

save_path = '/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/'
dir.create(paste0(save_path, 'raw/'), showWarnings = F, recursive = T)
dir.create(paste0(save_path, 'volcano/'), showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'raw/DA_BR1_b_age_discrete_level2.csv'))

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
volcano2 <- ggplot(DAIs, aes(x=logFC, y=-log10(FDR), fill=celltype)) +
  geom_point(pch=21, size=6, alpha=0.8) +
  geom_hline(yintercept=-log10(0.05),size=1.1, lty="dashed") +
  geom_text_repel(data=DAIs[DAIs$FDR<0.1,], label=DAIs[DAIs$FDR<0.1,"celltype"],aes(fill=NULL),
                  size=6, color="grey30",force=9) +
  scale_fill_manual(values=level2_colour_dictionary) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Log2(FC)")

pdf(paste0(save_path, 'volcano/DA_BR1_b_age_discrete_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
pdf(paste0(save_path, 'volcano/DA_BR1_b_age_discrete_level2_2.pdf'), width=12,height = 8)
volcano2
dev.off()





#TEST BRCA2 IMPACT:

#Clean up for samples we can test
age_given <- !(metadata_sub$patient_age %in% c('Unknown', 'n/a')) & !is.na(metadata_sub$patient_age)
parity_known <- !(metadata_sub$parity %in% c('unknown', 'Unknown', 'n/a')) & !is.na(metadata_sub$parity)
not_celltype_sorted <- !(metadata_sub$live_sorted_boolean %in% c('remove', 'Unknown')) 
not_bad_sampletype <- !(metadata_sub$sample_type %in% c('Organoid LP', 'Unknown'))
AR_BR2_samples <- metadata_sub$risk_status %in% c('AR', 'HR-BR2')
not_milk_samples <- !(metadata_sub$patientID %in% c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 'HMC4', 'HMC5', 
                                                    'HMC6', 'HMC7', 'HMC8', 'HMC9'))

good_cells <- age_given & parity_known & not_celltype_sorted & not_bad_sampletype & AR_BR2_samples & not_milk_samples
reed_dataset_cells <- metadata_sub$dataset=='HBCA'
tbl <- table(good_cells, reed_dataset_cells)
print('Cells with required metadata (and AR):')
print(tbl)
print('Percentage from Reed dataset:')
print(tbl[2,2] / (tbl[2,1] + tbl[2,2])) #0.4007739

metadata_test <- metadata_sub[good_cells, ]
metadata_test$patient_age_est <- mapvalues(metadata_test$patient_age, 
                                           from=c('O', 'Y'),
                                           to=c('51', '49'))
metadata_test$old_young <- strtoi(metadata_test$patient_age_est) >= 50  #TRUE is equivalent to Old (O), False = Young (Y)
metadata_test$parous <- 'Nulliparous'
metadata_test$parous[metadata_test$parity != '0'] <- 'Parous'

#issue with design matrix rank. HBCA dataset column is combination of sample_type columns! So simplify and do not block for sample_type
metadata_test$dataset <- paste0(metadata_test$dataset, '_', metadata_test$sample_type)

abundances <- table(metadata_test$map_celltype_map2reed, metadata_test$sampleID) 
abundances <- unclass(abundances) 
head(abundances)

extra.info <- metadata_test[match(colnames(abundances), metadata_test$sampleID),]
y.ab <- DGEList(abundances, samples=extra.info[extra.info$sampleID %in% colnames(abundances),])
y.ab

keep <- filterByExpr(y.ab, group=y.ab$samples$old_young)
y.ab <- y.ab[keep,] #old: #for now try keep in BSL2 cluster which this would remove
print(keep)

design <- model.matrix(~ dataset + live_sorted_boolean + parous + old_young + risk_status, y.ab$samples)
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

save_path = '/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/'
dir.create(paste0(save_path, 'raw/'), showWarnings = F, recursive = T)
dir.create(paste0(save_path, 'volcano/'), showWarnings = F, recursive = T)
write.csv(DAIs, file= paste0(save_path, 'raw/DA_BR2_b_age_discrete_level2.csv'))

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
volcano2 <- ggplot(DAIs, aes(x=logFC, y=-log10(FDR), fill=celltype)) +
  geom_point(pch=21, size=6, alpha=0.8) +
  geom_hline(yintercept=-log10(0.05),size=1.1, lty="dashed") +
  geom_text_repel(data=DAIs[DAIs$FDR<0.1,], label=DAIs[DAIs$FDR<0.1,"celltype"],aes(fill=NULL),
                  size=6, color="grey30",force=9) +
  scale_fill_manual(values=level2_colour_dictionary) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Log2(FC)")

pdf(paste0(save_path, 'volcano/DA_BR2_b_age_discrete_level2.pdf'), width=12,height = 8)
Volcano
dev.off()
pdf(paste0(save_path, 'volcano/DA_BR2_b_age_discrete_level2_2.pdf'), width=12,height = 8)
volcano2
dev.off()



