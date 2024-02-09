#R

#Complete statistical testing for milo analysis.

#Use the same neighborhoods for all AR vs HR tests
#scVI embedding



suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

# #example options (testing)

# #Epi
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/epi/output/prop0.3k50'


# #Str
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/str/output/prop0.3k50'


# #Imm
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/imm/output/prop0.3k50'


## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo object and to use for saving da results.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



#Functions
milo_stat_test <- function(milo, da_pwd, block_var, sample_tp, reduced_dimension){
  if (block_var == 'none') {
    design <- as.data.frame(colData(milo))[,c('sampleID', 'before', 'milo_test')]
    design <- distinct(design)
    rownames(design) <- design[, 'sampleID']
    da_out <- testNhoods(milo, 
                         design = ~ before + milo_test, 
                         design.df = design, 
                         fdr.weighting="graph-overlap",
                         reduced.dim = reduced_dimension)
  } else if (block_var == 'parity_age') {
    design <- as.data.frame(colData(milo))[,c('sampleID', 'before', 'milo_block1', 'milo_block2', 'milo_test')]
    design <- distinct(design)
    rownames(design) <- design[, 'sampleID']
    da_out <- testNhoods(milo, 
                         design = ~ before + milo_block1 + milo_block2 + milo_test, 
                         design.df = design, 
                         fdr.weighting="graph-overlap",
                         reduced.dim = reduced_dimension)
  } else {
    design <- as.data.frame(colData(milo))[,c('sampleID', 'before', 'milo_block', 'milo_test')]
    design <- distinct(design)
    rownames(design) <- design[, 'sampleID']
    da_out <- testNhoods(milo, 
                         design = ~ before + milo_block + milo_test, 
                         design.df = design, 
                         fdr.weighting="graph-overlap",
                         reduced.dim = reduced_dimension)
  }
  
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

for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM', 'patient_age', 'parity')){
  for (block_var in c('parity', 'patient_age', 'parity_age')){ 
    print('New round of testing:')
    print(test_var)
    print(block_var)
    
    if (block_var %in% c('none')){
      block_var_category <- 'none'
    } else {
      block_var_category <- 'parity'
    }
    
    #skip if test and block are not testable
    if (test_var == block_var) {
      next
    } else if ((test_var %in% c('patient_age', 'parity')) & (block_var == 'parity_age')) {
      next
    } else if ((test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')) & !(block_var %in% c('parity_age'))){
      next
    }
    
    if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')) {
      milo <- readRDS(paste0(prefix, '/milo_objects/HR_nghd_', block_var_category, '.rds'))
    } else {
      milo <- readRDS(paste0(prefix, '/milo_objects/AR_nghd_', block_var_category, '.rds'))
    } 
    
    
    #Testing
    
    #test_var meta
    if (test_var == 'WT_BRCA1PM') {
      milo <- milo[, milo$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1')]
      milo$milo_test <- milo$tissue_condition
    } else if (test_var == 'WT_BRCA2PM') {
      milo <- milo[, milo$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA2')]
      milo$milo_test <- milo$tissue_condition
    } else if (test_var == 'patient_age') {
      milo <- milo[, milo$tissue_condition %in% c('Mammoplasty WT')]
      milo$milo_test <- as.double(milo$patient_age)
    } else if (test_var == 'parity') {
      milo <- milo[, milo$tissue_condition %in% c('Mammoplasty WT')]
      milo <- milo[, milo$parity != 'unknown']
      milo$milo_test <- milo$parity != '0'
    } else {
      print(paste0('Incorrect input: test_var = ', test_var))
      next
    }
    #block_var meta
    if (block_var == 'none') {
      milo$milo_block <- 'none'
    } else if (block_var == 'patient_age') {
      milo$milo_block <- milo$patient_age
    } else if (block_var == 'parity') {
      milo <- milo[, milo$parity %in% c('0','1','2','3','4')]
      milo$milo_block[milo$parity %in% c('1','2','3','4')] <- 'Parous'
      milo$milo_block[milo$parity %in% c('0')] <- 'Nulliparous'
    } else if (block_var == 'parity_age') {
      milo <- milo[, milo$parity %in% c('0','1','2','3','4')]
      milo$milo_block1[milo$parity %in% c('1','2','3','4')] <- 'Parous'
      milo$milo_block1[milo$parity %in% c('0')] <- 'Nulliparous'
      milo$milo_block2 <- milo$patient_age
    }else {
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
