#R

#Prepare the Nee et al dataset for joining and integration in iHBCA.

#conda env milo




#libraries
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(plyr)
library(dplyr)
library(Matrix)




#load data
mat_list <- list()
mat_list[['ctrl1']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320152_scRNA_ctrl1_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl2']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320153_scRNA_ctrl2_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl3']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320154_scRNA_ctrl3_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl4']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320155_scRNA_ctrl4_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl5']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320156_scRNA_ctrl5_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl6']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320157_scRNA_ctrl6_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl7']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320158_scRNA_ctrl7_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl8']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320159_scRNA_ctrl8_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl9']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320160_scRNA_ctrl9_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl10']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320161_scRNA_ctrl10_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['ctrl11']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320162_scRNA_ctrl11_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca1']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320163_scRNA_brca1_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca2']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320164_scRNA_brca2_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca3']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320165_scRNA_brca3_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca4']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320166_scRNA_brca4_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca5']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320167_scRNA_brca5_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca6']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320168_scRNA_brca6_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca7']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320169_scRNA_brca7_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca8']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320170_scRNA_brca8_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca9']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320171_scRNA_brca9_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca10']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320172_scRNA_brca10_matrix.txt', sep=' ', header=TRUE, row.names=1)
mat_list[['brca11']] <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/GSM5320173_scRNA_brca11_matrix.txt', sep=' ', header=TRUE, row.names=1)



#make seperate counts, rowdata and coldata of set formatting
dir.create('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted')


#counts
nee_counts <- do.call(cbind, mat_list)

#fix barcodes
bcs <- sapply(strsplit(colnames(nee_counts), '\\.'), function(x) {
  if (x[2] == 'UCI_Patients'){
    return(paste0(x[2], ' ', x[3], '-', x[4]))
  } else {
    return(paste0(x[2], '-', x[3]))
  }
})

nee_counts <- as.matrix(nee_counts)
nee_counts <- as(nee_counts, "sparseMatrix")   

colnames(nee_counts) <- bcs
writeMM(nee_counts, file = '/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_counts.mtx')

#rowdata

nee_rowData <- data.frame(row.names = rownames(nee_counts),
                             'symbol' = rownames(nee_counts))
nee_rowData$gene_id <- mapIds(org.Hs.eg.db,
                                 keys=nee_rowData$symbol,
                                 column="ENSEMBL",
                                 keytype="SYMBOL",
                                 multiVals="first")
#rownames(nee_rowdata) <- nee_rowdata$gene_id #fails as some ensemble IDs are not found
print(head(nee_rowData))
write.csv(nee_rowData, '/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_features.csv')

#coldata

#Grab and format the coldata info
nee_colData <- read.table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/nee_etal_metadata.txt', sep=' ', header=TRUE, row.names=1)
nee_colData <- nee_colData[colnames(nee_counts), ]

#save a raw colData
write.csv(nee_colData, '/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_full_phenodata.csv')

#subset and save common colData
nee_colData_short <-   data.frame(row.names = rownames(nee_colData),
                                     'cellID' = rownames(nee_colData),
                                     'patientID' = nee_colData$patient_id,
                                     'batch' = nee_colData$patient_id,
                                     'level0' = mapvalues(nee_colData$cell_type_final,
                                                          from = c("Basal", "Endothelial", "Fibroblasts", "Immune", "Luminal1", "Luminal2", 
                                                                   "Lymphatic", "Pericytes"),
                                                          to = c("Epithelial", "Stroma", "Stroma", "Immune", "Epithelial", "Epithelial", 
                                                                 "Stroma", "Stroma")),
                                     'level1' = nee_colData$cell_type_final,
                                     'level2' = nee_colData$cell_state) 

print(head(nee_colData_short))
write.csv(nee_colData_short, '/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_phenodata.csv')
