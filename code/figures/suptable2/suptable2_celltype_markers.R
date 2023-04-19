# This script aims to perform the Differential expression testing between our different level1 labels and then the level2 labels within each level1 celltype.
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/figures/suptable2/suptable2_celltype_markers.R


### Load libraries
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(BiocParallel))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))


### Functions

# setup_summed_sce <- function(sce, celltype_subset){
#   if (celltype_subset != 'none'){
#     sce_sub <- sce[,sce$level1 %in% celltype_subset]
#     
#     #create pseudobulk object
#     summed <- aggregateAcrossCells(sce, ids=DataFrame(celltype=sce$level2,
#                                                       sampleID=sce$sampleID),
#                                    use.assay.type = "counts",
#                                    BPPARAM = MulticoreParam(workers = multicoreWorkers()))
#     
#   } else {
#     summed <- aggregateAcrossCells(sce, ids=DataFrame(celltype=sce$level1,
#                                                       sampleID=sce$sampleID),
#                                    use.assay.type = "counts",
#                                    BPPARAM = MulticoreParam(workers = multicoreWorkers()))
#   }
#   return(summed)
# }

setup_summed_sce <- function(sce, celltype_subset, celltype_level) {
  if (celltype_level %in% c('level2_nosub')) {
    celltype_level_column <- 'level2'
  } else {
    celltype_level_column <- celltype_level
  }
  
  if (celltype_subset != 'none') {
    sce <- sce[, sce$level1 %in% celltype_subset]
  } 
  
  sampleID_list <- unique(sce$sampleID)
  print(sampleID_list)
  
  i=1
  first <- TRUE
  while (5*i < 5 + length(sampleID_list)) {
    sampleID_to_use = sampleID_list[(5*i-4):min((5*i), length(sampleID_list))]
    print(i)
    print(sampleID_to_use)
    
    #sce$psuedobatch[sce$sampleID %in% sampleID_to_use] <- i
    sce_sub <- sce[, sce$sampleID %in% sampleID_to_use]
    sce_sub$psuedobatch <- i
    summed_sub <- aggregateAcrossCells(sce_sub, ids=DataFrame(celltype=colData(sce_sub)[, celltype_level_column],
                                                          sampleID=sce_sub$sampleID),
                                       use.assay.type = "counts",
                                       BPPARAM = MulticoreParam(workers = 16))
    #make combined sumed object
    if (first) {
      summed <- summed_sub
      first <- FALSE
    } else {
      summed <- cbind(summed, summed_sub)
    }
    
    i = i + 1
  }
  print('Done subsetted aggregation.')
  
  #combine the batches of summed into 1
  summed <- aggregateAcrossCells(summed, ids=DataFrame(celltype=colData(summed)[, celltype_level_column],
                                                       sampleID=summed$sampleID),
                                 use.assay.type = "counts",
                                 BPPARAM = MulticoreParam(workers = 16))
  
  return(summed)
}

test_1_vs_all <- function(current_summed, celltype_level, celltype, dge_pwd){
  #Only consider samples with >5 cells
  current_summed <- current_summed[,(current_summed$ncells > 5)]
  
  #setup test var
  current_summed$dge_test <- current_summed$celltype == celltype
  
  y <- DGEList(counts(current_summed), samples=colData(current_summed))
  
  design <- model.matrix(~ sample_type_coarse + dge_test, data = y$sample)
  print(head(design))
  
  keep <- filterByExpr(y, design=design)
  y <- y[keep,]
  
  #Normalize
  y <- calcNormFactors(y)
  
  #NB dispersion estimation
  y <- estimateDisp(y, design)
  
  #Plot biological coefficient of variation graph
  #plotBCV(y) # no plotting on server
  
  #Find QL neg. bin. fit
  fit <- glmQLFit(y, design, robust=TRUE)
  #plotQLDisp(fit) # no plotting on server
  
  #Run test
  res <- glmQLFTest(fit, coef=ncol(design))
  print(summary(decideTests(res)))
  
  #See top gene markers
  print(topTags(res))
  
  res$table$FDR <- p.adjust(res$table$PValue, method="BH")
  DEGs <- res$table[order(res$table$FDR),]
  DEGs_up <- DEGs[DEGs$logFC > 0,]
  DEGs_down <- DEGs[DEGs$logFC < 0,]
  print(paste0('Number of DEGs: ', sum(DEGs$FDR < 0.1))) 
  
  dir.create(paste0(dge_pwd, '/', celltype_level), showWarnings = FALSE, recursive = TRUE)
  write.csv(DEGs, file = paste0(dge_pwd, '/', celltype_level, '/dge-', celltype, '-all.csv'))
  write.csv(DEGs_up, file = paste0(dge_pwd, celltype_level, '/dge-', celltype, '-up.csv'))
  write.csv(DEGs_down, file = paste0(dge_pwd, celltype_level, '/dge-', celltype, '-down.csv'))
}




### Analysis 

#load in data
#sce <- readRDS('/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/bbknn/epi/output/input/input_sce.rds')
sce <- readRDS('/nfs/research/marioni/areed/projects/hbca/cellchat/2022-04-05/scvi_new/overview_all/output/input/HBCA_postlabelling.rds') #from cellchat prepare_sce.R script

#level2 labels
df <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
df_sub <- df[,c("cellID", "patientID", "sampleID", "tissue_condition", "before", "level1", "level2")]


dim(sce)
sum(sce$cellID %in% df_sub$cellID)

df_sub <- df_sub[df_sub$cellID %in% sce$cellID,]
df_sub <- df_sub[order(match(df_sub$cellID, sce$cellID)),]
level1_labels <- df_sub$level1
level1_labels[level1_labels %in% c('Lymphoid', 'Myeloid')] <- 'Immune'  ##neccesary as Macrophage is the only myeloid celltype
level2_labels <- df_sub$level2

sum(sce$cellID == df_sub$cellID)

sce$level1 <- level1_labels
sce$level2 <- level2_labels
# sce$level2[(sce$level1 %in% c("Lymphoid", "Myeloid"))] <- sce$level1[(sce$level1 %in% c("Lymphoid", "Myeloid"))]

#clean up celltypes
celltypes_to_ignore <- c('Doublet', '')
sce <- sce[, !(sce$level2 %in% celltypes_to_ignore)]

#change from ensemble_ID to gene_name
rownames(sce) <- rowData(sce)$X



#Do testing
for (celltype_level in c('level2', 'level1', 'level2_nosub')) {
  print(celltype_level)
  
  #First we need to subset to the groups of celltypes we are testing within
  if (celltype_level %in% c('level1')) {
    subset_list <- 'none'
    subsetornot <- 'nosubset'
  } else if (celltype_level %in% c('level2_nosub')) {
    subset_list <- 'none'
    subsetornot <- 'nosubset'
  } else {
    subset_list <-  unique(colData(sce)[, 'level1'])
    subsetornot <- 'subset'
  } 
  
  for (celltype_subset in subset_list) {
    print(celltype_subset)
    
    summed <- setup_summed_sce(sce = sce, celltype_subset = celltype_subset,
                               celltype_level = celltype_level)
    
    #Now we iterate over the relevant (sub)cell types to obtain markers.
    celltypes_to_test <- unique(summed$celltype)
    print(celltypes_to_test)
    
    for (celltype in celltypes_to_test) {
      print(celltype)
      
      test_1_vs_all(current_summed = summed, celltype_level = celltype_level, 
                    celltype = celltype, 
                    dge_pwd = paste0('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/' , subsetornot, '/'))
    } 
  }
}







