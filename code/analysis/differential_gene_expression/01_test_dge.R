#! complete differential gene expression testing.
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/dge/scvi_new/test_dge.R

suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))


# Example inputs

# #LP
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi/output/input/input_sce.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/lp/output/'
# opt$celltype_subset = 'LP1,LP2,LP3,LP4,LP_proliferating' 

#BSL
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi/output/input/input_sce.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/bsl/output/'
# opt$celltype_subset = 'BSL1,BSL2'


## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve sce object.'), 
                   make_option(c('--output_pwd'),
                               type='character',
                               help='Pathway to directory used for saving dge results.'),
                   make_option(c('--celltype_subset'),
                               type='character',
                               help='Selection of cells to use for testing (comma seperated labels). Use "None" for no subsetting'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# Functions
doDGE <- function(summed, dge_pwd, block_var){
  #this function will test for differentially expressed genes across our predefined dge_test variable
  
  samples_to_test <- unique(summed$sampleID[(summed$dge_test != 'NA') & !is.na(summed$dge_test)])
  current <- summed[, summed$sampleID %in% samples_to_test]
  
  #print(current$ncells)
  #print(table(current$level2,current$sampleID))
  
  # PREVIOUSLY USED IN MOUSE ATLAS TESTING - not require here?
  # Remove any samples that are more than two standard deviations below the mean counts 
  # OR less than 20 counts
  # for (celltype in unique(current$level2)){
  #   #want upper 97.5% (of assumed normal distribution) giving the 1.96 sd's.
  #   threshold <- mean(current$ncells[current$level2 == celltype]) - 1.96 * sd(current$ncells[current$level2 == celltype])
  #   #need to account for case when sd returns NA due to only one valid cell.
  #   if (is.na(threshold)){
  #     threshold <- 0
  #   }
  #   current <- current[,((current$ncells > max(threshold, 20)) & (current$level2 == celltype)) | (current$level2 != celltype)] 
  # }
  
  #Only consider samples with >5 cells
  current <- current[,(current$ncells > 5)]
  
  #print(table(current$level2,current$sampleID))
  
  #Make DGEList object
  y <- DGEList(counts(current), samples=colData(current))
  
  # FROM MICE TESTING NOT NEEDED NOW?
  # #Remove samples where no comparison can be made across Age
  # for (celltype in unique(current$level2)){
  #   temp1 <- sum(y$samples$level2 == celltype & y$samples$dge_test == 0)
  #   temp2 <- sum(y$samples$level2 == celltype & y$samples$dge_test == 1)
  #   if ((temp1 == 0) | (temp2 == 0)){
  #     y <- y[,y$samples$level2 != celltype]
  #     current <- current[, current$level2 != celltype]
  #   }
  # }
  
  #current$level2 <- unique(current$level2)
  #y$samples$level2 <- unique(y$samples$level2)
  #print(table(unique(current$level2), unique(current$sampleID)))  
  
  if (dim(current)[2] == 0){
    print('Too few cells of this type')
    return()
  }
  
  ##Stat modeling
  if (length(unique(summed$level2)) == 1) {
    if (block_var == 'none') {
      design <- model.matrix(~ before + dge_test, data = y$sample)
    } else if (block_var == 'parity_age'){
      design <- model.matrix(~ before + dge_block1 + dge_block2 + dge_test, data = y$sample) 
    } else {
      design <- model.matrix(~ before + dge_block + dge_test, data = y$sample)
    }
  } else{
    if (block_var == 'none') {
      design <- model.matrix(~ before + level2 + dge_test, data = y$sample)
    } else if (block_var == 'parity_age'){
      design <- model.matrix(~ before + dge_block1 + dge_block2 + dge_test, data = y$sample)
    } else {
      design <- model.matrix(~ before + level2 + dge_block + dge_test, data = y$sample)
    }
  }
  print(head(design))
  
  #Removing lowly expressed genes
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
  
  dir.create(dge_pwd, showWarnings = FALSE, recursive = TRUE)
  write.csv(DEGs, file = paste0(dge_pwd, '/dge-block_', block_var, '-all.csv'))
  write.csv(DEGs_up, file = paste0(dge_pwd, '/dge-block_', block_var, '-up.csv'))
  write.csv(DEGs_down, file = paste0(dge_pwd, '/dge-block_', block_var, '-down.csv'))
  
  return()
}

setup_and_test <- function(summed, test_var, block_var, prefix){
  #set up test_var coldata 
  #(used similarly to that in milo_statistics_testing_explore.R)
  
  if (test_var == 'parity') {
    summed <- summed[, summed$parity != 'unknown']
    summed$dge_test[summed$parity %in% c('1','2','3','4')] <- 'Parous'
    summed$dge_test[summed$parity %in% c('0')] <- 'Nulliparous'
  } else if (test_var == 'patient_age') {
    summed$dge_test <- summed$patient_age
  } else if (test_var == 'WT_BRCA1PM') {
    summed <- summed[, summed$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1')]
    summed$dge_test <- summed$tissue_condition
  } else if (test_var == 'WT_BRCA2PM') {
    summed <-summed[,summed$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA2')]
    summed$dge_test <-summed$tissue_condition
  } else if (test_var == 'BRCA1PM_BRCA1C') {
    summed <- summed[, summed$tissue_condition %in% c('Contralateral BRCA1', 'Mastectomy BRCA1')]
    summed$dge_test <- summed$tissue_condition
  } else if (test_var == 'WT_BRCA1all') {
    summed <- summed[, summed$tissue_condition %in% c('Mammoplasty WT', 'Contralateral BRCA1', 'Mastectomy BRCA1')]
    summed$dge_test[summed$tissue_condition %in% c('Mammoplasty WT')] <- 'WT'
    summed$dge_test[summed$tissue_condition %in% c('Contralateral BRCA1', 'Mastectomy BRCA1')] <- 'BRCA1'
    summed$dge_test <- factor(summed$dge_test, levels = c('WT', 'BRCA1'))
  } else if (test_var == 'Menopause_status') {
    summed <- summed[, !is.na(summed$Menopause_status)]
    summed$dge_test <- summed$Menopause_status
    summed$dge_test[summed$dge_test == 'Post_(surgically_induced)'] <- 'Post'
  } else if (test_var == 'Body_mass_index') {
    summed <- summed[, !is.na(summed$Body_mass_index)]
    summed$dge_test <- as.double(summed$Body_mass_index)
  } else if (test_var == 'Smoking_status') {
    summed <- summed[, !is.na(summed$Smoking_status)]
    summed$dge_test <- summed$Smoking_status
  } else if (test_var == 'HRT_use') {
    summed <- summed[, !is.na(summed$HRT_use)]
    summed$dge_test <- mapvalues(summed$HRT_use, from =c('Never', 'Past', 'Present'), to = c(FALSE, TRUE, TRUE))
  } else if (test_var == 'OCP_use') {
    summed <- summed[, !is.na(summed$OCP_use)]
    summed$dge_test <- mapvalues(summed$OCP_use, from =c('Never', 'No', 'Past'), to = c(FALSE, FALSE, TRUE))
  } else if (test_var == 'Image_scans') {
    summed <- summed[, !is.na(summed$Image_scans)]
    summed$dge_test <- as.double(summed$Image_scans)
  } else {
    return(paste0('Incorrect input: test_var = ', test_var))
  }
  
  #setup block_var coldata
  if (block_var == 'parity') {
    summed$dge_block <- 'Unk' #I dont think I can just leave the unknown as NA/NULL so I will make them their own category
    summed$dge_block[summed$parity %in% c('1','2','3','4')] <- 'Parous'
    summed$dge_block[summed$parity %in% c('0')] <- 'Nulliparous'
  } else if (block_var == 'patient_age') {
    summed$dge_block <- summed$patient_age
  } else if (block_var == 'tissue_condition') {
    summed$dge_block <- summed$tissue_condition
  } else if (block_var == 'parity_age') {
    summed$dge_block1 <- 'Unk' #I dont think I can just leave the unknown as NA/NULL so I will make them their own category
    summed$dge_block1[summed$parity %in% c('1','2','3','4')] <- 'Parous'
    summed$dge_block1[summed$parity %in% c('0')] <- 'Nulliparous'
    summed$dge_block2 <- summed$patient_age
  } else if (block_var == 'none') {
    summed$dge_block <- 'none'
  } else {
    print(paste0('Incorrect input: block_var = ', block_var))
    return()
  }
  
  print('testing:')
  print(test_var)
  print(block_var)
  
  doDGE(summed = summed,
        dge_pwd = paste0(prefix, '/dge_testing/', test_var),
        block_var = block_var) 
}



# Analysis 

#load in data
sce <- readRDS(opt$input_path)
rownames(sce) <- rowData(sce)$X

#clean up celltypes
celltypes_to_ignore <- c('Doublet')
sce <- sce[, !(sce$level2 %in% celltypes_to_ignore)]

if (opt$celltype_subset != 'None') {
  celltypes_to_test <- strsplit(opt$celltype_subset, split= ',')[[1]]
  celltypes_to_test_nounderscore <- sub('_', ' ', celltypes_to_test)
  sce <- sce[, sce$level2 %in% c(celltypes_to_test, celltypes_to_test_nounderscore)]
}

#Extra coldata (BMI etc.) already in colData(sce)

#create pseudobulk object
summed <- aggregateAcrossCells(sce, ids=DataFrame(level2=sce$level2,
                                                  sampleID=sce$sampleID))

#make test/block variable lists to test over
test_var_list <- list('parity', 'patient_age', 'WT_BRCA1PM', 'BRCA1PM_BRCA1C', 'WT_BRCA2PM', "Menopause_status", "Body_mass_index", 'Smoking_status', "HRT_use", "OCP_use", "Image_scans")
block_var_list <- list('none', 'parity_age', 'parity', 'patient_age', 'tissue_condition')

#complete tests
for (test_var in test_var_list){
  for (block_var in block_var_list){
    print('New round of test:')
    print(test_var)
    print(block_var)
    if (test_var == block_var) {
      next
    } else if ((test_var %in% c('parity', 'patient_age')) & (block_var == 'parity_age')) {
      next
    } else if ((block_var == 'tissue_condition') & (test_var %in% c('WT_BRCA1PM', 'WT_BRCA1all', 'BRCA1PM_BRCA1C', 'WT_BRCA2PM'))) {
      next
    } else if ((block_var != 'none') & (test_var %in% c("Menopause_status", "Body_mass_index", 'Smoking_status', "HRT_use", "OCP_use", "Image_scans"))) { #these have too few patients to consider blocking
      next
    } else {
      print('test = TRUE')
      setup_and_test(summed = summed, test_var = test_var, block_var = block_var,
                     prefix = opt$output_pwd)
    }
  }
}


#test inputs
# test_var = 'parity'
# block_var = 'none'
# prefix = opt$output_pwd
# dge_pwd = paste0(prefix, '/dge_testing/', test_var)

