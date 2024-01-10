#R

#Prepare the Gray et al dataset for joining and integration in iHBCA.

#conda env milo



#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(plyr)
library(dplyr)
library(Matrix)

#load data
gray_counts <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/original/GSE180878_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_human.csv',
                          sep=',', header = TRUE, row.names = 1)
gray_meta <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/original/GSE180878_Li_Brugge_10XscRNAseq_Metadata_human.csv', 
                        sep=',', header=TRUE)

#make seperate counts, rowdata and coldata of set formatting
dir.create('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted')


#rowdata
gray_rowData <- data.frame(row.names = rownames(gray_counts),
                           'symbol' = rownames(gray_counts))
gray_rowData$gene_id <- mapIds(org.Hs.eg.db,
                                  keys=gray_rowData$symbol,
                                  column="ENSEMBL",
                                  keytype="SYMBOL",
                                  multiVals="first")
#rownames(gray_rowdata) <- gray_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(gray_rowData))
write.csv(gray_rowData, '/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/gray_features.csv')

#counts
colnames(gray_counts) <- sapply(strsplit(colnames(gray_counts), '\\.'), function(x) paste0(x[1], '-', x[2]))
gray_counts <- as.matrix(gray_counts)
gray_counts <- as(gray_counts, "sparseMatrix")
writeMM(gray_counts, file = '/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/gray_counts.mtx')


#coldata

#Grab and format the coldata info
rownames(gray_meta) <- gray_meta$cellID
gray_meta <- gray_meta[colnames(gray_counts), ]

patients <- sapply(strsplit(gray_meta$cellID, '_'), function(x) x[1])
gray_colData <- data.frame(row.names = gray_meta$cellID,
                           'cellID' = gray_meta$cellID,
                           'patientID' = patients,
                           'batch' = patients,
                           'level0' = mapvalues(gray_meta$Major.subtype, 
                                                from = c("AV", "BA", "HS", "Fibroblast", "VasLymph", "Immune"),
                                                to = c("Epithelial", "Epithelial", "Epithelial", "Stroma", "Stroma", "Immune")),
                           'level1' = gray_meta$Major.subtype,
                           'level2' = gray_meta$Cell.subtype)

print(head(gray_colData))
write.csv(gray_colData, '/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/gray_phenodata.csv')
