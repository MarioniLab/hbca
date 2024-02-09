#R

#Make more plots for milo analysis

#Use the same neighborhoods for all AR vs HR tests
#scVI embedding

#libraries
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

# #example options (testing)

# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo outputs.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#celltype vectors

celltypes <- c("LASP1", "LASP2", "LASP3", "LASP4", "LASP5",
               "LHS1", "LHS2", "LHS3",
               "BMYO1", "BMYO2", 
               "DDC1", "DDC2",
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
               'Macro', 'DC')
celltypes_reduced <- celltypes[!(celltypes %in% c('DDC1', 'DDC2', 'Macro-lipo', 'LASP5'))]

# Analysis
prefix <- opt$input_path
for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM', 'patient_age', 'parity')){  

  if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')){
    block_var <- 'parity_age'
  } else if (test_var == 'patient_age'){
    block_var <- 'parity'
  } else if (test_var %in% c('parity', 'WT_WTPM-unk2')) {
    block_var <- 'patient_age'
  } else {
    block_var <- 'none'
  }
  
  if (block_var %in% c('none')){
    block_var_category <- 'none'
  } else {
    block_var_category <- 'parity'
  }
  
  print('New round of plotting:')
  print(test_var)
  print(block_var)
  
  milo_list <- list()
  da_list <- list()
  
  for (celltype_load in c('epi', 'str', 'imm')) {
    if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')) {
      milo_list[[celltype_load]] <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/HR_nghd_', block_var_category, '.rds'))
    } else {
      milo_list[[celltype_load]] <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/AR_nghd_', block_var_category, '.rds'))
    } 
    da_list[[celltype_load]] <- read.csv(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_testing/', test_var, '/da-block_', block_var, '-all.csv'))
    
    da_list[[celltype_load]] <- annotateNhoods(milo_list[[celltype_load]], 
                                               da_list[[celltype_load]], 
                                               coldata_col = "level2")
    
  }
  
  #join the dataframes
  da_joined <- do.call(rbind, da_list)
  
  da_joined$level2 <- factor(da_joined$level2,
                             levels = rev(celltypes))
  
  dir.create(paste0(prefix, 'join/output/data/', test_var, '/'), showWarnings = F, recursive = T)
  write.csv(da_joined, file = paste0(prefix, 'join/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
  # da_joined <- read.csv(paste0(prefix, 'join/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
  
  #make agrregated df to plot the mean values
  da_join_agg <- aggregate(x = da_joined,
                           by = list(da_joined$level2), 
                           FUN = mean)
  
  da_join_agg$level2 <- factor(da_join_agg$Group.1,
                               levels = rev(celltypes))
  
  #Remove celltypes not present across all types in sufficient numbers
  da_joined <- da_joined[da_joined$level2 %in% celltypes_reduced,]
  da_join_agg <- da_join_agg[da_join_agg$level2 %in% celltypes_reduced,]
  
  if (test_var == 'patient_age'){
    max_x = 0.3
    min_x = -0.3
  } else if (test_var == 'parity') {
    max_x = 8
    min_x = -8
  } else {
    max_x = 5
    min_x = -5
  }
  
  dir.create(paste0(prefix, 'join/output/', test_var, '/'), showWarnings = F, recursive = T)
  if (sum(da_joined$SpatialFDR < 0.1) > 0) {
    pdf(paste0(prefix, 'join/output/', test_var, '/beeswarm_level2-block_', block_var, '-all.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level2")+
            scale_y_continuous(limits=c(min_x, max_x)) + 
            theme(axis.text = element_text(size = 15))) 
    dev.off()
  } else {
    pdf(paste0(prefix, 'join/output/', test_var, '/beeswarm_level2-block_', block_var, '-all_alpha=1.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_joined, group.by = "level2", alpha=1)+
            scale_y_continuous(limits=c(min_x, max_x)) + 
            theme(axis.text = element_text(size = 15)))
    dev.off()
  }
  
}

