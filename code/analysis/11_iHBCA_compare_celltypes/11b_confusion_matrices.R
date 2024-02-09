#Plot nice confusion matricies from logistic regression mappings


#conda env milo

## libraries
library(scran)
library(scater)
library(umap)
library(viridis)

library(plyr)
library(dplyr)
library(ggplot2)
library(pheatmap)



#data
general_save_pw <- '/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/'
colData <- read.csv(paste0(general_save_pw, 'output_final/colData/integrated_unlabelled.csv'))


dir.create(paste0(general_save_pw, 'output_final/redux_heatmap/'), showWarnings = F, recursive = T)

##ORDERING LABELS FOR VISUALISATION EASE

order_of_celltypes_list <- list('reed' = c("LASP1", "LASP2", "LASP3", "LASP4", "LASP5",
                                           "LHS1", "LHS2", "LHS3",
                                           "BMYO1", "BMYO2", 
                                           "FB1", "FB2", "FB3", "FB4",
                                           "PV1", "PV2", "PV3", "PV4", "PV5",
                                           "VEV", "VEC", "VEA","VEAT", 
                                           "LE1", "LE2",
                                           'CD4_naive', "CD4_Th", 
                                           'CD8_Tem', 'CD8_Trm',
                                           'CD8_Tc1',
                                           'NKT', 'NK',
                                           'ILC', 
                                           'B_naive', 'B_mem_switched', 'B_mem_unswitched', 
                                           'Plasma_cell', 
                                           'Macro', 'Macro-lipo', 'DC'),
                                'kumar' = c('Lumsec-major','Lumsec-basal', 'Lumsec-myo', 
                                            'Lumsec-KIT', 'Lumsec-lac', 'Lumsec-HLA', 
                                            'Lumsec-prol',
                                            'LummHR-major', 'LummHR-active', 'LummHR-SCGB',
                                            'basal', 
                                            'Fibro-major', 'Fibro-matrix', 'Fibro-prematrix',
                                            'Fibro-SFRP4',
                                            'pericytes', 'vsmc',
                                            'Vas-arterial', 'Vas-capillary', 'Vas-venous',
                                            'Lymph-major', 'Lymph-immune', 'Lymph-valve1',
                                            'Lymph-valve2', 
                                            'CD4-activated', 'CD4-naive', 'CD4-Tem', 
                                            'CD4-Th', 'CD4-Th-like', 'CD4-Treg', 
                                            'CD8-activated', 'CD8-Tem', 'CD8-Trm',
                                            'T_prol', 'GD',
                                            'NK', 'NKT', 'NK-ILCs',
                                            'b_naive', 'bmem_switched', 'bmem_unswitched', 
                                            'plasma_IgA', 'plasma_IgG', 
                                            'Macro-IFN', 'Macro-lipo', 'Macro-m1', 
                                            'Macro-m1-CCL', 'Macro-m2', 'Macro-m2-CXCL', 
                                            'Mono-classical', 'Mono-non-classical',
                                            'mDC', 'pDC', 'cDC1', 'cDC2',  
                                            'Neutrophil',  
                                            'Mast', 
                                            'mye-prol'),
                                'nee' = c('Luminal1-ALDH1A3', 'Luminal1-LTF',
                                          'Luminal2-AREG', 'Luminal2-MUCL1',
                                          'Basal', 'Basal-Myoepithelial',
                                          'Fibroblasts',
                                          'Pericytes',
                                          'Endothelial',
                                          'Lymphatic',
                                          'Immune'),
                                'gray' = c('AP', 'BL', 
                                           'HSa', 'HSb', 
                                           'BAa', 'BAb',
                                           'F1', 'F2', 'F3', 
                                           'VL1_LE', 'VL2_VE', 'VL3_PE',
                                           'I3_Tcell', 'I2_NK',
                                           'I4_Bcell', 'I5_PlasmaCell', 
                                           'I1_Myeloid'),
                                'murrow' = c('Secretory_Luminal', 'HRpos_Luminal', 'Basal', 
                                             'Fibroblast',
                                             'Vascular_Accessory',
                                             'Vascular_Endothelial', 'Lymphatic_Endothelial',
                                             'Lymphocyte', 
                                             'Macrophage'),
                                'twigger' = c('LC1', 'LC2', 'LP', 'HR', 'BA', 
                                              'FB',
                                              'VA',
                                              'EN',
                                              'IM'),
                                'pal' = c('Epithelial',
                                          'Luminal Progenitor', 'Mature Luminal', 'Basal', 
                                          'Fibroblast',
                                          'Stroma')
)


#map on the logstic regression predicted labels
dataset_names <- unique(colData$dataset)
for (dataset_from in c(dataset_names)) {
  print(dataset_from)
  
  dataset_from_low <- tolower(dataset_from)
  
  if (dataset_from == 'HBCA'){
    dataset_from_low <- 'reed'
  }
  
  temp <- read.csv(paste0(general_save_pw, 'output_final/celltypist_demo_folder/predictions_meta_from_', dataset_from_low, '.csv'))
  rownames(temp) <- temp$X
  temp <- temp[colData$X[colData$dataset!=dataset_from],]
  colData[, paste0('map_celltype_map2', dataset_from_low)] <- paste0('True-', colData$celltype)
  colData[colData$dataset!=dataset_from, paste0('map_celltype_map2', dataset_from_low)] <- paste0(dataset_from, '-', temp$predicted_labels)
}


#save data
dir.create(paste0(general_save_pw, 'output_final/colData/'), showWarnings = F, recursive = T)
write.csv(colData, paste0(general_save_pw, 'output_final/colData/integrated_labelled.csv'))


#remove dataset labels for clearer plotting
dataset_names <- unique(colData$dataset)
for (dataset_from in c(dataset_names)) {
  dataset_from_low <- tolower(dataset_from)
  if (dataset_from == 'HBCA'){
    dataset_from_low <- 'reed'
  }
  
  dataset_tag <- paste0(dataset_from, '-')
  
  colData[, paste0('map_celltype_map2', dataset_from_low)] <- gsub(dataset_tag, '', colData[, paste0('map_celltype_map2', dataset_from_low)])
  colData[, paste0('map_celltype_map2', dataset_from_low)] <- gsub('True-', '', colData[, paste0('map_celltype_map2', dataset_from_low)])
}


#make plots
for (dataset_from in c(dataset_names)){
  for (dataset_to in dataset_names){
    print(dataset_from)
    print(dataset_to)
    
    dataset_from_low <- tolower(dataset_from)
    dataset_to_low <- tolower(dataset_to)
    
    if (dataset_from == 'HBCA'){
      dataset_from_low <- 'reed'
    }
    if (dataset_to == 'HBCA'){
      dataset_to_low <- 'reed'
    }
    
    df <- colData[colData$dataset==dataset_from,]
    
    tbl <- table(df[, paste0('map_celltype_map2', dataset_from_low)], df[, paste0('map_celltype_map2', dataset_to_low)])
    rw_names <- order_of_celltypes_list[[dataset_from_low]][order_of_celltypes_list[[dataset_from_low]] %in% rownames(tbl)]
    cl_names <- order_of_celltypes_list[[dataset_to_low]][order_of_celltypes_list[[dataset_to_low]] %in% colnames(tbl)]
    tbl <- tbl[rw_names, cl_names]
    
    #normalise rows to 0-1 range
    tbl <- as.matrix(tbl)
    tbl_norm <- tbl
    for (row in rownames(tbl)) {
      row_sum <- sum(tbl[row,])
      print(row_sum)
      tbl_norm[row,] = tbl[row,] / row_sum
    }
    
    pdf(paste0(general_save_pw, paste0('output_final/redux_heatmap/heatmap_', dataset_from_low, '2', dataset_to_low, '.pdf')))
    print(pheatmap(tbl_norm, cluster_rows=F, cluster_cols=F, scale='none',
                   color=hcl.colors(50, "Lajolla"), border_color='#00000004', fontsize=7))
    dev.off()
  }
}

