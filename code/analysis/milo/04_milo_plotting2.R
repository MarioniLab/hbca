#! make EXTRA plots for milo analysis. Join ST (sample_types)
#This time set up to use the same neighborhoods for all AR vs HR tests. 
#scVI new

suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(optparse))
suppressMessages(library(MASS))
suppressMessages(library(ggrastr))

#example input
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/' #epi_2/output/prop0.3k50'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo outputs.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Functions



# Analysis
prefix <- opt$input_path

for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM')){
  
  if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM', 'WT_WTPM-unk')){
    block_var <- 'parity_age'
  } else if (test_var == 'patient_age'){
    block_var <- 'parity'
  } else if (test_var %in% c('parity', 'WT_WTPM-unk2')) {
    block_var <- 'patient_age'
  } else {
    block_var <- 'none'
  }
  
  if (block_var %in% c('none', 'patient_age')){
    block_var_category <- 'none'
  } else {
    block_var_category <- 'parity'
  }
  
  print('New round of plotting:')
  print(test_var)
  print(block_var)
  
  milo_list <- list()
  da_list <- list()
  
  for (celltype_load in c('epi_2', 'str_2', 'imm_2')) {
    milo_list[[celltype_load]] <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/epi_nghd_', block_var_category, '.rds'))
    da_list[[celltype_load]] <- read.csv(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_testing/', test_var, '/da-block_', block_var, '-all.csv'))
    
    da_list[[celltype_load]] <- annotateNhoods(milo_list[[celltype_load]], 
                                               da_list[[celltype_load]], 
                                               coldata_col = "level1")
    da_list[[celltype_load]] <- annotateNhoods(milo_list[[celltype_load]], 
                                               da_list[[celltype_load]], 
                                               coldata_col = "level2")
    
  }
  
  da_joined <- do.call(rbind, da_list)
  da_joined$level1 <- factor(da_joined$level1,
                             levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
                                            "Fibroblast", "Endothelial", "Vascular Mural",
                                            "Lymphoid", "Myeloid")))
  da_joined$level2 <- factor(da_joined$level2,
                             levels = rev(c("LP1", "LP2", "LP3", "LP4",
                                            "LP5", "HS1", "HS2", "HS3",
                                            "HS4", "BSL1", "BSL2",
                                            "DDC1", "DDC2",
                                            "FB1", "FB2", "FB3", "FB4", "FB5",
                                            "VM1", "VM2", "VM3", "VM4", "VM5",
                                            "EC venous", "EC capillary", "EC arterial", #B9AF4A
                                            "EC angiogenic tip", "LEC1", "LEC2",
                                            'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                                            'IFNG+ T', 'NK1', 'NK2', 'NK3',
                                            'ILC3', 'B cell', 'Plasma cell', 'Macrophage')))
  
  dir.create(paste0(prefix, 'join_2/output/data/', test_var, '/'), showWarnings = F, recursive = T)
  write.csv(da_joined, file = paste0(prefix, 'join_2/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
  
  #make agrregated df to plot the mean values
  da_join_agg <- aggregate(x = da_joined,
                           by = list(da_joined$level2), 
                           FUN = mean)
  
  da_join_agg$level1 <- factor(da_join_agg$level1,
                               levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
                                              "Fibroblast", "Endothelial", "Vascular Mural",
                                              "Lymphoid", "Myeloid")))
  da_join_agg$level2 <- factor(da_join_agg$Group.1,
                               levels = rev(c("LP1", "LP2", "LP3", "LP4", 
                                              "LP5", "HS1", "HS2", "HS3",
                                              "HS4", "BSL1", "BSL2",
                                              "DDC1", "DDC2",
                                              "FB1", "FB2", "FB3", "FB4", "FB5",
                                              "VM1", "VM2", "VM3", "VM4", "VM5",
                                              "EC venous", "EC capillary", "EC arterial", #B9AF4A
                                              "EC angiogenic tip", "LEC1", "LEC2",
                                              'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                                              'IFNG+ T', 'NK1', 'NK2', 'NK3',
                                              'ILC3', 'B cell', 'Plasma cell', 'Macrophage')))
  
  #Remove celltypes not present across all types in sufficient numbers
  da_joined <- da_joined[da_joined$level2 %in% c("LP1", "LP2", "LP3", "LP4", 
                                                 "HS1", "HS2", "HS3", #no "LP5",
                                                 "HS4", "BSL1", "BSL2",
                                                 # no "DDC1", "DDC2",
                                                 "FB1", "FB2", "FB3", "FB4", "FB5",
                                                 "VM1", "VM2", "VM3", "VM4", "VM5",
                                                 "EC venous", "EC capillary", "EC arterial", 
                                                 "EC angiogenic tip", "LEC1", "LEC2",
                                                 'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                                                 'IFNG+ T', 'NK1', 'NK2', 'NK3',
                                                 'ILC3', 'B cell', 'Plasma cell', 'Macrophage'),]
  da_join_agg <- da_join_agg[da_join_agg$level2 %in% c("LP1", "LP2", "LP3", "LP4", 
                                                       "HS1", "HS2", "HS3", #no "LP5",
                                                       "HS4", "BSL1", "BSL2",
                                                       # no "DDC1", "DDC2",
                                                       "FB1", "FB2", "FB3", "FB4", "FB5",
                                                       "VM1", "VM2", "VM3", "VM4", "VM5",
                                                       "EC venous", "EC capillary", "EC arterial",
                                                       "EC angiogenic tip", "LEC1", "LEC2",
                                                       'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                                                       'IFNG+ T', 'NK1', 'NK2', 'NK3',
                                                       'ILC3', 'B cell', 'Plasma cell', 'Macrophage'),]
  
  dir.create(paste0(prefix, 'join_2/output/', test_var, '/'), showWarnings = F, recursive = T)
  if (sum(da_joined$SpatialFDR < 0.1) > 0) {
    pdf(paste0(prefix, 'join_2/output/', test_var, '/beeswarm_level1-block_', block_var, '-all.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level1") +
            scale_y_continuous(limits=c(-5, 5)) + 
            theme(axis.text = element_text(size = 15)))
    dev.off()
    pdf(paste0(prefix, 'join_2/output/', test_var, '/beeswarm_level2-block_', block_var, '-all.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level2")+
            scale_y_continuous(limits=c(-5, 5)) + 
            theme(axis.text = element_text(size = 15))) 
    dev.off()
  } else {
    pdf(paste0(prefix, 'join_2/output/', test_var, '/beeswarm_level1-block_', block_var, '-all_alpha=1.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level1", alpha=1)+
            scale_y_continuous(limits=c(-5, 5)) + 
            theme(axis.text = element_text(size = 15)))
    dev.off()
    pdf(paste0(prefix, 'join_2/output/', test_var, '/beeswarm_level2-block_', block_var, '-all_alpha=1.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level2", alpha=1)+
            scale_y_continuous(limits=c(-5, 5)) + 
            theme(axis.text = element_text(size = 15)))
    dev.off()
  }
  
}



