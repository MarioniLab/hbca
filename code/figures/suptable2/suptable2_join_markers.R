#Join celltype marker genes into sing data frames


suppressMessages(library(plyr))
suppressMessages(library(dplyr))


#Functions


#Analysis

metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')

for (celltype_level in c('level2_nosub', 'level1', 'level2')) {
  print(celltype_level)
  
  #First we need to subset to the groups of celltypes we are testing within
  if (celltype_level %in% c('level1')) {
    subsetornot <- 'nosubset'
    celltype_label = 'level1'
  } else if (celltype_level %in% c('level2_nosub')) {
    subsetornot <- 'nosubset'
    celltype_label = 'level2'
  } else {
    subsetornot <- 'subset'
    celltype_label = 'level2'
  }
  
  if (celltype_label == 'level1') {
    celltypes_loop <- c( "Luminal progenitor", "Luminal hormone sensing", "Basal", 
                         "Fibroblast", "Endothelial", "Vascular mural", 
                         "Lymphoid", "Myeloid")
  } else if (celltype_label == 'level2') {
    celltypes_loop <- c("LP1", "LP2", "LP3", "LP4", 
                        "LP proliferating", "HS1", "HS2", "HS3",
                        "HS4", "BSL1", "BSL2",
                        "FB1", "FB2", "FB3", "FB4", "FB5",
                        "VM1", "VM2", "VM3", "VM4", "VM5",
                        "EC venous", "EC capillary", "EC arterial", 
                        "EC angiogenic tip", "LEC1", "LEC2",
                        'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T naive',
                        'IFNG+ T', 'NK1', 'NK2', 'NK3',
                        'ILC3', 'B cell', 'Plasma cell', 'Macrophage')
  }
  
  #make joint dataframe
  marker_df <- data.frame(matrix(ncol = 0, nrow = 100)) #data.frame(rownames = 1:100)
  dge_pwd <- paste0('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/' , subsetornot, '/')
  for (celltype in celltypes_loop){
    print(celltype)
    dge <- read.csv(paste0(dge_pwd, celltype_level, '/dge-', celltype, '-up.csv'))
    
    top100_genes <- dge[1:min(100, sum(dge$FDR < 0.05)),1]
    marker_df[1:min(100, sum(dge$FDR < 0.05)), celltype] <- top100_genes
  }
  
  dir.create(paste0('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/joined/'), 
             recursive = T, showWarnings = F)
  write.csv(marker_df, file = paste0('/nfs/research/marioni/areed/projects/hbca/figures/suptable2/preliminary/dge/joined/markers_', celltype_level, '.csv'))
}