#! complete differential gene expression testing.
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/dge/scvi_revision/revision_test_dge.R

suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))


# Example inputs

# #LP
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/lp/output/'
# opt$celltype_subset = 'LASP1,LASP2,LASP3,LASP4,LASP5' 

#LHS
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/hs/output/'
# opt$celltype_subset = 'LHS1,LHS2,LHS3'

#BSL
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/bsl/output/'
# opt$celltype_subset = 'BMYO1,BMYO2'

#Macro
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/imm/sce_imm_annotated.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/macro/output/'
# opt$celltype_subset = 'Macro'

#DC
# opt = list()
# opt$input_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/imm/sce_imm_annotated.rds'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/dc/output/'
# opt$celltype_subset = 'DC'

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
  
  #OLD - Only consider samples with >5 cells (except for DC which are too rare). Ignore DC as they are too rare in dataset to test. only a few per sample max.
  # if (opt$celltype_subset != 'DC'){
  #     current <- current[,(current$ncells > 5)]
  # }
  current <- current[,(current$ncells > 5)]
  
  #print(table(current$level2,current$sampleID))
  
  #Make DGEList object
  y <- DGEList(counts(current), samples=colData(current))
  
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
    summed$dge_block <- 'Unk' #Make them their own category
    summed$dge_block[summed$parity %in% c('1','2','3','4')] <- 'Parous'
    summed$dge_block[summed$parity %in% c('0')] <- 'Nulliparous'
  } else if (block_var == 'patient_age') {
    summed$dge_block <- summed$patient_age
  } else if (block_var == 'tissue_condition') {
    summed$dge_block <- summed$tissue_condition
  } else if (block_var == 'parity_age') {
    summed$dge_block1 <- 'Unk' #Make them their own category
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
test_var_list <- list('WT_BRCA1PM', 'WT_BRCA2PM')
block_var_list <- list('parity_age')

#complete tests
for (test_var in test_var_list){
  for (block_var in block_var_list){
    print('New round of test:')
    print(test_var)
    print(block_var)
    setup_and_test(summed = summed, test_var = test_var, block_var = block_var,
                   prefix = opt$output_pwd)
  }
}

