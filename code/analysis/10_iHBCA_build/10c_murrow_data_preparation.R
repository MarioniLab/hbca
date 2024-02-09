#R

#Prepare the Murrow et al dataset for joining and integration in iHBCA.

#conda env milo


#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(plyr)
library(dplyr)
library(Matrix)




#load data
srt <- readRDS('/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/original/GSE198732_breast.data.rds')

#make seperate counts, rowdata and coldata of set formatting
dir.create('/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted')

#rowdata
murrow_rowData <- data.frame(row.names = rownames(srt@assays$RNA@counts),
                           'symbol' = rownames(srt@assays$RNA@counts))
murrow_rowData$gene_id <- mapIds(org.Hs.eg.db,
                               keys=murrow_rowData$symbol,
                               column="ENSEMBL",
                               keytype="SYMBOL",
                               multiVals="first")
#rownames(murrow_rowdata) <- murrow_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(murrow_rowData))
write.csv(murrow_rowData, '/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_features.csv')

#coldata

#Grab and format the coldata info
murrow_colData <- srt@meta.data

#save a raw colData
write.csv(murrow_colData, '/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_full_phenodata.csv')

#subset and save common colData
murrow_colData_short <-   data.frame(row.names = rownames(murrow_colData),
                                     'cellID' = rownames(murrow_colData),
                                     'patientID' = murrow_colData$Sample,
                                     'batch' = murrow_colData$Batch,
                                     'level0' = mapvalues(murrow_colData$Type,
                                                          from = c("HRpos_Luminal", "Basal", "Secretory_Luminal", 
                                                                   "Fibroblast",  "Vascular_Endothelial", "Vascular_Accessory",
                                                                   "Lymphocyte", "Lymphatic_Endothelial", "Macrophage"),
                                                          to = c("Epithelial", "Epithelial", "Epithelial", "Stroma", "Stroma", "Stroma", "Immune", "Stroma", "Immune")),
                                     'level1' = murrow_colData$Type,
                                     'level2' = NA)

print(head(murrow_colData_short))
write.csv(murrow_colData_short, '/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_phenodata.csv')

#counts
murrow_counts <- srt@assays$RNA@counts
writeMM(murrow_counts, file = '/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_counts.mtx')
