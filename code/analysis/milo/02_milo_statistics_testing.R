#! complete statistical testing for milo analysis.
#This time set up to use the same neighborhoods for all AR vs HR tests. 
#scVI new


suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# #example input

# #Epi
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi_2/output/prop0.3k50'


# #Str
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/str_2/output/prop0.3k50'


# #Imm
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/imm_2/output/prop0.3k50'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo object and to use for saving da results.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



#Functions
milo_stat_test <- function(milo, da_pwd, block_var, sample_tp, reduced_dimension){
  design <- as.data.frame(colData(milo))[,c('sampleID', 'before', 'milo_block1', 'milo_block2', 'milo_test')]
  design <- distinct(design)
  rownames(design) <- design[, 'sampleID']
  da_out <- testNhoods(milo, 
                       design = ~ before + milo_block1 + milo_block2 + milo_test, 
                       design.df = design, 
                       fdr.weighting="graph-overlap",
                       reduced.dim = reduced_dimension)
  
  dir.create(da_pwd, showWarnings = FALSE, recursive = TRUE)
  write.csv(da_out, paste0(da_pwd, 'da-block_', block_var, '-', sample_tp, '.csv'))
  
  pdf(paste0(da_pwd, 'pvalue_histogram-block_', block_var, '-', sample_tp, '.pdf'), 
      width=12, height = 8)
  print(ggplot(da_out, aes(PValue)) + geom_histogram(bins=50))
  dev.off()
  
  pdf(paste0(da_pwd, 'nghd_volcano-block_', block_var, '-', sample_tp, '.pdf'), width=12,height = 8)
  print(ggplot(da_out, aes(logFC, -log10(SpatialFDR))) +
          geom_point() +
          geom_hline(yintercept = 1)) ## Mark significance threshold (10% FDR)
  dev.off()
  }



## Analysis

##choose the proportion/knn input with optparse

#load Data
prefix <- opt$input_path 

for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM')){
  for (block_var in c('parity_age')){
    print('New round of testing:')
    print(test_var)
    print(block_var)
    
    if (block_var %in% c('none', 'patient_age')){
      block_var_category <- 'none'
    } else {
      block_var_category <- 'parity'
    }
    
    
    milo <- readRDS(paste0(prefix, '/milo_objects/epi_nghd_', block_var_category, '.rds'))

    #test order is okay (previously had mixed up row labels)
    dir.create(paste0(prefix, '/test_umap/', test_var, '/'), showWarnings = FALSE, recursive = TRUE)
    pdf(paste0(prefix, '/test_umap/', test_var, '/test_umap_', block_var_category, '.pdf'))
    print(plotReducedDim(milo, 'UMAP', colour_by = 'level1'))
    dev.off()
    
    #Testing
    
    #test_var meta
    if (test_var == 'WT_BRCA1PM') {
      milo <- milo[, milo$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1')]
      milo$milo_test <- milo$tissue_condition
    } else if (test_var == 'WT_BRCA2PM') {
      milo <-milo[,milo$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA2')]
      milo$milo_test <-milo$tissue_condition
    } else if (test_var == 'WT_WTPM') {
      milo <- milo[,milo$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy WT')]
      milo$milo_test <- milo$tissue_condition
    } 
    
    #block_var meta
    if (block_var == 'none') {
      milo$milo_block <- 'none'
    } else if (block_var == 'parity_age') {
      milo <- milo[, milo$parity %in% c('0','1','2','3','4')]
      milo$milo_block1[milo$parity %in% c('1','2','3','4')] <- 'Parous'
      milo$milo_block1[milo$parity %in% c('0')] <- 'Nulliparous'
      milo$milo_block2 <- milo$patient_age
    } else {
      print(paste0('Incorrect input: block_var = ', block_var))
      next
    }
    
    milo_stat_test(milo = milo,
                   da_pwd = paste0(prefix, '/milo_testing/', test_var, '/'),
                   block_var = block_var,
                   sample_tp = 'all',
                   reduced_dimension = 'scVI')
  }
}

















