#R

#Prepare the Twigger et al dataset for joining and integration in iHBCA.

#conda env milo




#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(Seurat)

library(plyr)
library(dplyr)
library(Matrix)

#load data
sce <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/original/sce_all_nospike_2.rds')


#make seperate counts, rowdata and coldata of set formatting
dir.create('/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted')

#rowdata
twigger_rowData <- data.frame(#row.names = rownames(sce),
                              'symbol' = rownames(sce))
twigger_rowData$gene_id <- mapIds(org.Hs.eg.db,
                              keys=twigger_rowData$symbol,
                              column="ENSEMBL",
                              keytype="SYMBOL",
                              multiVals="first")

#identify the duplicate gene symbols
n_occur <- data.frame(table(twigger_rowData$symbol))
# twigger_rowData[twigger_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1],]
dupped_symbols <- unique(twigger_rowData[twigger_rowData$symbol %in% n_occur$Var1[n_occur$Freq > 1], 'symbol'])
write.csv(dupped_symbols, '/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/duplicated_symbols.csv')

temp_counts <- counts(sce)
#temp_sub <- temp_counts[rownames(temp_counts) %in% dupped_symbols,]
temp_sums <- data.frame('symbols' = names(rowSums(temp_counts)),
                            'sums' = rowSums(temp_counts),
                            'keep' = TRUE)
for (symbol in dupped_symbols){
  temp <- temp_sums[temp_sums$symbols %in% symbol,]
  max_val <- max(temp$sums)
  temp$keep <- temp$sums == max_val
  temp_sums[temp_sums$symbols %in% symbol,] <- temp
}

#use the keep column to remove duplicated genes.
twigger_rowData <- twigger_rowData[temp_sums$keep,]
rownames(twigger_rowData) <- twigger_rowData$symbol

write.csv(twigger_rowData, '/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/twigger_features.csv')


#coldata
twigger_meta <- colData(sce)
twigger_colData <- data.frame(row.names = twigger_meta$Barcode,
                              'cellID' = twigger_meta$Barcode,
                              'patientID' = twigger_meta$Sample,
                              'batch' = twigger_meta$Batches,
                              'level0' = mapvalues(twigger_meta$Identity, 
                                                   from = c("FB", "LP", "LC1", "LC2", "IM", 
                                                            "EN", "BA", "VA", "HR"),
                                                   to = c("Stroma", "Epithelial", "Epithelial", "Epithelial", "Immune",
                                                          "Stroma", "Epithelial", "Stroma", "Epithelial")),
                              'level1' = twigger_meta$Identity,
                              'level2' = NA)
print(head(twigger_colData))
write.csv(twigger_colData, '/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/twigger_phenodata.csv')

#counts
twigger_counts <- counts(sce)
twigger_counts <- twigger_counts[temp_sums$keep, ] #remove dupped genes
writeMM(twigger_counts, file = '/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/twigger_counts.mtx')
