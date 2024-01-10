#R

#Prepare the Pal et al dataset for joining and integration in iHBCA.

#conda env milo


#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)

library(plyr)
library(dplyr)
library(Matrix)


##Read in data

#counts matrixes (per lane/batch)
mat_list <- list()
mat_list[['PM0092_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909253_N-PM0092-Total-matrix.mtx.gz')
mat_list[['PM0019_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909254_N-PM0019-Total-matrix.mtx.gz')
mat_list[['N280_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909255_N-N280-Epi-matrix.mtx.gz')
mat_list[['PM0095_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909256_N-PM0095-Epi-matrix.mtx.gz')
mat_list[['PM0095_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909257_N-PM0095-Total-matrix.mtx.gz')
mat_list[['NF_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909258_N-NF-Epi-matrix.mtx.gz')
mat_list[['NE_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909259_N-NE-Epi-matrix.mtx.gz')
mat_list[['N1105_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909260_N-N1105-Epi-matrix.mtx.gz')
mat_list[['PM0230_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909261_N-PM0230-Total-matrix.mtx.gz')
mat_list[['MH0064_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909262_N-MH0064-Epi-matrix.mtx.gz')
mat_list[['MH0064_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909263_N-MH0064-Total-matrix.mtx.gz')
mat_list[['N1B_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909264_N-N1B-Epi-matrix.mtx.gz')
mat_list[['PM0233_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909265_N-PM0233-Total-matrix.mtx.gz')
mat_list[['MH0169_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909266_N-MH0169-Total-matrix.mtx.gz')
mat_list[['MH0023_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909267_N-MH0023-Epi-matrix.mtx.gz')
mat_list[['MH0023_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909268_N-MH0023-Total-matrix.mtx.gz')
mat_list[['PM0342_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909269_N-PM0342-Epi-matrix.mtx.gz')
mat_list[['PM0342_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909270_N-PM0342-Total-matrix.mtx.gz')
mat_list[['MH288_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909271_N-MH288-Total-matrix.mtx.gz')
mat_list[['MH0021_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909272_N-MH0021-Total-matrix.mtx.gz')
mat_list[['MH275_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909273_N-MH275-Epi-matrix.mtx.gz')
mat_list[['MH275_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909274_N-MH275-Total-matrix.mtx.gz')
mat_list[['PM0372_epi']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909275_N-PM0372-Epi-matrix.mtx.gz')
mat_list[['PM0372_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909276_N-PM0372-Total-matrix.mtx.gz')
mat_list[['KCF0894_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909277_B1-KCF0894-matrix.mtx.gz')
mat_list[['MH0033_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909278_B1-MH0033-matrix.mtx.gz')
mat_list[['MH0023_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909279_B1-MH0023-matrix.mtx.gz')
mat_list[['MH0090_mix']] <- readMM('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909280_B1-MH0090-matrix.mtx.gz')

#do same for barcodes:
bc_list <- list()
bc_list[['PM0092_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909253_N-PM0092-Total-barcodes.tsv.gz', sep='\t')
bc_list[['PM0019_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909254_N-PM0019-Total-barcodes.tsv.gz', sep='\t')
bc_list[['N280_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909255_N-N280-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0095_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909256_N-PM0095-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0095_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909257_N-PM0095-Total-barcodes.tsv.gz', sep='\t')
bc_list[['NF_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909258_N-NF-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['NE_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909259_N-NE-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['N1105_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909260_N-N1105-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0230_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909261_N-PM0230-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH0064_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909262_N-MH0064-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['MH0064_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909263_N-MH0064-Total-barcodes.tsv.gz', sep='\t')
bc_list[['N1B_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909264_N-N1B-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0233_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909265_N-PM0233-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH0169_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909266_N-MH0169-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH0023_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909267_N-MH0023-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['MH0023_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909268_N-MH0023-Total-barcodes.tsv.gz', sep='\t')
bc_list[['PM0342_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909269_N-PM0342-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0342_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909270_N-PM0342-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH288_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909271_N-MH288-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH0021_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909272_N-MH0021-Total-barcodes.tsv.gz', sep='\t')
bc_list[['MH275_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909273_N-MH275-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['MH275_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909274_N-MH275-Total-barcodes.tsv.gz', sep='\t')
bc_list[['PM0372_epi']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909275_N-PM0372-Epi-barcodes.tsv.gz', sep='\t')
bc_list[['PM0372_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909276_N-PM0372-Total-barcodes.tsv.gz', sep='\t')
bc_list[['KCF0894_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909277_B1-KCF0894-barcodes.tsv.gz', sep='\t')
bc_list[['MH0033_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909278_B1-MH0033-barcodes.tsv.gz', sep='\t')
bc_list[['MH0023_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909279_B1-MH0023-barcodes.tsv.gz', sep='\t')
bc_list[['MH0090_mix']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSM4909280_B1-MH0090-barcodes.tsv.gz', sep='\t')

#features
features <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/GSE161529_features.tsv.gz', sep='\t')


#Seurat objects
srt_normb1total <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/SeuratObject_NormB1Total.rds')
srt_normb1total@meta.data$cell_type <- mapvalues(srt_normb1total$seurat_clusters, 
                                               from = c('0','1','2','3','4','5','6','7', '8', '9'), 
                                               to = c('normBr1_Stroma', 'normBr1_Epithelial', 'normBr1_Stroma', 'normBr1_Epithelial', 
                                                      'normBr1_Epithelial', 'normBr1_Stroma', 'normBr1_Stroma', 'normBr1_Epithelial', 
                                                      'normBr1_Stroma', 'normBr1_Stroma'))

srt_normepi <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/SeuratObject_NormEpi.rds')
srt_normepi@meta.data$cell_type <- mapvalues(srt_normepi$seurat_clusters, 
                                               from = c('0','1','2','3'), 
                                               to = c('normEpi_LP', 'normEpi_ML', 'normEpi_Basal', 'normEpi_Fb'))

srt_normtotal <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/SeuratObject_NormTotal.rds')
srt_normtotal@meta.data$cell_type <- mapvalues(srt_normtotal$seurat_clusters, 
                                               from = c('0','1','2','3','4','5','6','7'), 
                                               to = c('normTotal_LP', 'normTotal_Stroma', 'normTotal_Basal', 'normTotal_ML', 'normTotal_Stroma', 
                                                      'normTotal_Stroma', 'normTotal_Stroma', 'normTotal_Stroma'))

#get rid of overlap normal cells in the b1 dataset
srt_normb1total <- srt_normb1total[, srt_normb1total@meta.data$orig.ident=='B1']

#merge
srt_combined <- merge(srt_normtotal, srt_normepi, add.cell.ids = c("", ""), project='hbca')
srt_combined <- merge(srt_combined, srt_normb1total, add.cell.ids = c("", ""), project='hbca')

#saveRDS of combined seurat
saveRDS(srt_combined, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/seuratobject_combined.rds')
# srt_combined <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/seuratobject_combined.rds')

#mapping sample names across srt and GEO raw datasets
srt_batch <- paste0(sapply(strsplit(srt_combined@meta.data$cell_type, '_'), function(x) x[1]),
                    '_',
                    sapply(strsplit(srt_combined@meta.data$group, '_'), function(x) paste0(x[1], '_', x[2])))
srt_combined@meta.data$srt_batch <- srt_batch
srt_combined@meta.data$cellID <- gsub('^__|^_', '', rownames(srt_combined@meta.data))
srt_combined@meta.data$barcode_10x <- sapply(strsplit(srt_combined@meta.data$cellID, '_'), function(x) x[4])
srt_combined@meta.data$barcode_10x[is.na(srt_combined@meta.data$barcode_10x)] <- sapply(strsplit(srt_combined@meta.data$cellID, '_'), function(x) x[3])[is.na(srt_combined@meta.data$barcode_10x)]

GEO_batch_unq <- names(bc_list)
srt_batch_unq <- unique(srt_batch)

srt_bc_list <- list()
for (batch in srt_batch_unq){
  srt_bc_list[[batch]] <- srt_combined@meta.data$barcode_10x[srt_combined@meta.data$srt_batch == batch]
}

df <- data.frame(matrix(data=0, ncol = length(GEO_batch_unq), nrow = length(srt_batch_unq)))
rownames(df) <- srt_batch_unq
colnames(df) <- GEO_batch_unq
df2 <- data.frame(row.names=GEO_batch_unq, 'srt_mapping' = rep(NA, length(GEO_batch_unq)))

for (rw in rownames(df)) {
  for (cl in colnames(df)) {
    test_var <- srt_bc_list[[rw]] %in% bc_list[[cl]]$V1
    #df[rw, cl] <- sum(test_var)
    if (sum(test_var) == length(srt_bc_list[[rw]])) {
      df[rw, cl] <- sum(test_var) 
      df2[cl,1] <- rw
    }
  }
}

write.csv(df2, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_GEO_to_srt_sampleID_mapping.csv')

srt_combined@meta.data$GEO_sampleID <- mapvalues(srt_combined@meta.data$srt_batch, from = df2$srt_mapping, to=rownames(df2))
srt_combined@meta.data$GEO_cellID <- paste0(srt_combined@meta.data$GEO_sampleID, '_', srt_combined@meta.data$barcode_10x)

saveRDS(srt_combined, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/seuratobject_combined2.rds')




## Make seperate counts, rowdata and coldata of set formatting

dir.create('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted', showWarnings = FALSE, recursive = TRUE)


#counts
for (batch in GEO_batch_unq) {
  rownames(mat_list[[batch]]) <- features$V2
  colnames(mat_list[[batch]]) <- paste0(batch, '_', bc_list[[batch]]$V1)
}
pal_counts <- do.call(cbind, mat_list)

#save after cleaning genes and cells up:

#rowdata
pal_rowData <- data.frame(#row.names = rownames(pal_counts),
                          'symbol' = rownames(pal_counts),
                          'gene_id' = features$V1)

#identify the duplicate gene symbols
n_occur <- data.frame(table(pal_rowData$symbol))
# pal_rowData[pal_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1],]
dupped_symbols <- unique(pal_rowData[pal_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1], 'symbol'])
write.csv(dupped_symbols, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/duplicated_symbols.csv')

temp_counts <- pal_counts
#temp_sub <- temp_counts[rownames(temp_counts) %in% dupped_symbols,]
temp_sums <- data.frame('symbols' = names(rowSums(temp_counts)),
                        'sums' = rowSums(temp_counts),
                        'keep' = TRUE)
for (symbol in dupped_symbols){
  temp <- temp_sums[temp_sums$symbols %in% symbol,]
  # print(symbol)
  # print(temp)
  max_val <- max(temp$sums)
  
  temp$keep <- (temp$sums == max_val) & (max_val > 0)
  temp_sums[temp_sums$symbols %in% symbol,] <- temp
  
}

#use the keep column to remove duplicated genes.
pal_rowData <- pal_rowData[temp_sums$keep,]
rownames(pal_rowData) <- pal_rowData$symbol

print(head(pal_rowData))
write.csv(pal_rowData, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_features.csv')


#coldata

#Grab and format the coldata info
pal_colData <- srt_combined@meta.data

#save a raw colData
write.csv(pal_colData, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_full_phenodata.csv')

#subset and save common colData
pal_colData_short <- data.frame(row.names = pal_colData$GEO_cellID,
                                'cellID' = pal_colData$GEO_cellID,
                                'patientID' = sapply(strsplit(pal_colData$GEO_sampleID, '_'), function(x) x[1]),
                                'batch' = pal_colData$GEO_sampleID,
                                'level0' = mapvalues(pal_colData$cell_type,
                                                     from = c("normBr1_Epithelial", "normBr1_Stroma",
                                                              "normEpi_Basal", "normEpi_Fb", "normEpi_LP", "normEpi_ML",
                                                              "normTotal_Basal", "normTotal_LP", "normTotal_ML", "normTotal_Stroma"),
                                                     to = c("Epithelial", "Stroma", 
                                                            "Epithelial", "Stroma", "Epithelial", "Epithelial", 
                                                            "Epithelial", "Epithelial", "Epithelial", "Stroma")),
                                'level1' = mapvalues(pal_colData$cell_type,
                                                     from = c("normBr1_Epithelial", "normBr1_Stroma",
                                                              "normEpi_Basal", "normEpi_Fb", "normEpi_LP", "normEpi_ML",
                                                              "normTotal_Basal", "normTotal_LP", "normTotal_ML", "normTotal_Stroma"),
                                                     to = c("Epithelial", 'Stroma', 
                                                            "Basal", "Fibroblast", "Luminal Progenitor", "Mature Luminal", 
                                                            "Basal", "Luminal Progenitor", "Mature Luminal", "Stroma")),
                                'level2' = NA)

pal_colData_short <- pal_colData_short[colnames(pal_counts), ]
pal_colData_short <- pal_colData_short[!is.na(pal_colData_short$batch),] #remove those without counts.

print(head(pal_colData_short))
write.csv(pal_colData_short, '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_phenodata.csv')


#save counts
pal_counts <- pal_counts[temp_sums$keep, pal_colData_short$cellID]
writeMM(pal_counts, file = '/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_counts.mtx')
