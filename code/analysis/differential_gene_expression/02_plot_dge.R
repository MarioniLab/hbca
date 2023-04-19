#! plot differential gene expression results from 'test_dge.R'.
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/dge/scvi_new/plot_dge.R

suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(EnhancedVolcano))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(ggrastr))


# Example inputs

# #LP
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi/output/input/input_sce.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/lp/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/lp/output/'
# opt$celltype_subset = 'LP1,LP2,LP3,LP4,LP_proliferating'

#BSL
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi/output/input/input_sce.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/bsl/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/bsl/output/'
# opt$celltype_subset = 'BSL1,BSL2'


## options
option_list = list(make_option(c('--input_sce_path'),
                               type='character',
                               help='Pathway to retrieve sce object.'), 
                   make_option(c('--input_dge_pwd'),
                               type='character',
                               help='Pathway to directory used for saving dge results.'),
                   make_option(c('--output_pwd'),
                               type='character',
                               help='Pathway to directory used for saving output plots.'),
                   make_option(c('--celltype_subset'),
                               type='character',
                               help='Selection of cells to use for testing (comma seperated labels). Use "None" for no subsetting'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Functions
doVolcano <- function(DEGs, test_var, block_var, FCcutoff = 1, maxoverlapsConnectors = 20, no_label=FALSE){
  ttl <- paste0('DEGs testing ', test_var, ' - block: ', block_var, '.')
  
  if (test_var == 'patient_age') {
    FCcutoff <- 0.01 #approximately equivelent to 1 logFC over lifetime
  }
  if (no_label) {
    Volcano <- EnhancedVolcano(DEGs,
                               lab = NA,
                               x = 'logFC',
                               y = 'FDR',
                               title = ttl,
                               subtitle = '',
                               subtitleLabSize = 2,
                               legendPosition = "bottom",
                               pointSize = 3.0,
                               FCcutoff = FCcutoff,
                               pCutoff = 10e-2,
                               col = c("grey", "forestgreen", "steelblue", "red"),
                               #legendVisible = FALSE,
                               drawConnectors = FALSE,
                               typeConnectors = 'open',
                               #raster=TRUE
                               )
  } else {
    Volcano <- EnhancedVolcano(DEGs,
                               lab = DEGs$X,
                               x = 'logFC',
                               y = 'FDR',
                               title = ttl,
                               subtitle = '',
                               subtitleLabSize = 2,
                               legendPosition = "bottom",
                               pointSize = 3.0,
                               labSize = 2.0,
                               FCcutoff = FCcutoff,
                               pCutoff = 10e-2,
                               col = c("grey", "forestgreen", "steelblue", "red"),
                               #legendVisible = FALSE,
                               drawConnectors = TRUE,
                               typeConnectors = 'open',
                               maxoverlapsConnectors = maxoverlapsConnectors)
  }
  
  return(Volcano)
}

setup_and_plot <- function(sce, DEGs, test_var, block_var, prefix){
  volcano_plot1 <- doVolcano(DEGs = DEGs,
                             test_var = test_var,
                             block_var = block_var)
  volcano_plot1 <- rasterize(volcano_plot1, layers='GeomPoint', dpi=500)
  volcano_plot2 <- doVolcano(DEGs = DEGs,
                             test_var = test_var,
                             block_var = block_var,
                             no_label = TRUE)
  volcano_plot2 <- rasterize(volcano_plot2, layers='GeomPoint', dpi=500)
  
  volcano_pwd = paste0(prefix, '/dge_plotting/', test_var)
  dir.create(volcano_pwd, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(volcano_pwd, '/Volcano-', test_var, '-block_', block_var, '.pdf'), width=12,height = 8)
  print(volcano_plot1)
  dev.off()
  pdf(paste0(volcano_pwd, '/Volcano-', test_var, '-block_', block_var, '_nolabel.pdf'), width=12,height = 8)
  print(volcano_plot2)
  dev.off()

  #plot UMAPs for some of the top dge genes
  num_genes=10
  genes_up <- DEGs[DEGs$logFC>0,]$X[1:num_genes]
  genes_down <- DEGs[DEGs$logFC<0,]$X[1:num_genes]

  # dir.create(paste0(prefix, '/umap/', test_var, '/'), showWarnings = FALSE, recursive = TRUE)
  # for (i in 1:num_genes) {
  #   png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-UpReg_', i, '_', genes_up[i], '.png'))
  #   print(plotReducedDim(sce, 'UMAP', colour_by = genes_up[i], text_by='level2'))
  #   dev.off()
  #   png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-DownReg_', i, '_', genes_down[i], '.png'))
  #   print(plotReducedDim(sce, 'UMAP', colour_by = genes_down[i], text_by='level2'))
  #   dev.off()
  # }
}
 


# Analysis 

# #testing
# test_var='WT_BRCA1all'
# block_var='parity'


#load in data
sce <- readRDS(opt$input_sce_path)
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


#make test/block variable lists to test over
test_var_list <- list('parity', 'patient_age', 'WT_BRCA1PM', 'BRCA1PM_BRCA1C', 'WT_BRCA2PM', "Menopause_status", "Body_mass_index", 'Smoking_status', "HRT_use", "OCP_use", "Image_scans")
block_var_list <- list('none', 'parity_age', 'parity', 'patient_age', 'tissue_condition')

#complete tests
for (test_var in test_var_list){
  for (block_var in block_var_list){
    print('New round of plot:')
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
      print('plotting = TRUE')
      DEGs <- read.csv(paste0(opt$input_dge_pwd, '/dge_testing/', test_var, '/dge-block_', block_var, '-all.csv'))
      setup_and_plot(sce = sce, DEGs = DEGs, test_var = test_var, block_var = block_var,
                     prefix = opt$output_pwd)
    }
  }
}


#test inputs
# test_var = 'parity'
# block_var = 'none'
# prefix = opt$output_pwd
# dge_pwd = paste0(prefix, '/dge_plotting/', test_var)



###copying

# setup_and_plot <- function(sce, DEGs, test_var, block_var, prefix){
#   volcano_plot1 <- doVolcano(DEGs = DEGs,
#                              test_var = test_var,
#                              block_var = block_var) 
#   volcano_pwd = paste0(prefix, '/dge_plotting/', test_var)
#   dir.create(volcano_pwd, showWarnings = FALSE, recursive = TRUE)
#   pdf(paste0(volcano_pwd, '/Volcano-', test_var, '-block_', block_var, '.pdf'), width=12,height = 8)
#   print(volcano_plot1)
#   dev.off()
#   num_genes=10
#   genes_up <- DEGs[DEGs$logFC>0,]$X[1:num_genes]
#   genes_down <- DEGs[DEGs$logFC<0,]$X[1:num_genes]
#   dir.create(paste0(prefix, '/umap/', test_var, '/'), showWarnings = FALSE, recursive = TRUE)
#   for (i in 1:num_genes) {
#     png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-UpReg_', i, '_', genes_up[i], '.png'))
#     print(plotReducedDim(sce, 'UMAP', colour_by = genes_up[i], text_by='level2'))
#     dev.off()
#     png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-DownReg_', i, '_', genes_down[i], '.png'))
#     print(plotReducedDim(sce, 'UMAP', colour_by = genes_down[i], text_by='level2'))
#     dev.off()
#   }
# }


##OLD

# setup_and_plot <- function(sce, DEGs, test_var, block_var, prefix){
#   if (test_var == 'parity') {
#     sce <- sce[, sce$parity != 'unknown']
#     sce$dge_test[sce$parity %in% c('1','2','3','4')] <- 'Parous'
#     sce$dge_test[sce$parity %in% c('0')] <- 'Nulliparous'
#   } else if (test_var == 'patient_age') {
#     sce$dge_test <- sce$patient_age
#   } else if (test_var == 'WT_BRCA1PM') {
#     sce <- sce[, sce$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1')]
#     sce$dge_test <- sce$tissue_condition
#   } else if (test_var == 'WT_BRCA2PM') {
#     sce <-sce[,sce$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA2')]
#     sce$dge_test <-sce$tissue_condition
#   } else if (test_var == 'BRCA1PM_BRCA1C') {
#     sce <- sce[, sce$tissue_condition %in% c('Contralateral BRCA1', 'Mastectomy BRCA1')]
#     sce$dge_test <- sce$tissue_condition
#   } else if (test_var == 'WT_BRCA1all') {
#     sce <- sce[, sce$tissue_condition %in% c('Mammoplasty WT', 'Contralateral BRCA1', 'Mastectomy BRCA1')]
#     sce$dge_test[sce$tissue_condition %in% c('Mammoplasty WT')] <- 'WT'
#     sce$dge_test[sce$tissue_condition %in% c('Contralateral BRCA1', 'Mastectomy BRCA1')] <- 'BRCA1'
#     sce$dge_test <- factor(sce$dge_test, levels = c('WT', 'BRCA1'))
#   } else if (test_var == 'Menopause_status') {
#     sce <- sce[, !is.na(sce$Menopause_status)]
#     sce$dge_test <- sce$Menopause_status
#     sce$dge_test[sce$dge_test == 'Post_(surgically_induced)'] <- 'Post'
#   } else if (test_var == 'Body_mass_index') {
#     sce <- sce[, !is.na(sce$Body_mass_index)]
#     sce$dge_test <- as.double(sce$Body_mass_index)
#   } else if (test_var == 'Smoking_status') {
#     sce <- sce[, !is.na(sce$Smoking_status)]
#     sce$dge_test <- sce$Smoking_status
#   } else if (test_var == 'HRT_use') {
#     sce <- sce[, !is.na(sce$HRT_use)]
#     sce$dge_test <- mapvalues(sce$HRT_use, from =c('Never', 'Past', 'Present'), to = c(FALSE, TRUE, TRUE))
#   } else if (test_var == 'OCP_use') {
#     sce <- sce[, !is.na(sce$OCP_use)]
#     sce$dge_test <- mapvalues(sce$OCP_use, from =c('Never', 'No', 'Past'), to = c(FALSE, FALSE, TRUE))
#   } else if (test_var == 'Image_scans') {
#     sce <- sce[, !is.na(sce$Image_scans)]
#     sce$dge_test <- as.double(sce$Image_scans)
#   } else {
#     return(paste0('Incorrect input: test_var = ', test_var))
#   }
#   
#   #setup block_var coldata
#   if (block_var == 'parity') {
#     sce$dge_block <- 'Unk' #I dont think I can just leave the unknown as NA/NULL so I will make them their own category
#     sce$dge_block[sce$parity %in% c('1','2','3','4')] <- 'Parous'
#     sce$dge_block[sce$parity %in% c('0')] <- 'Nulliparous'
#   } else if (block_var == 'patient_age') {
#     sce$dge_block <- sce$patient_age
#   } else if (block_var == 'tissue_condition') {
#     sce$dge_block <- sce$tissue_condition
#   } else if (block_var == 'parity_age') {
#     summed$dge_block1 <- 'Unk' #I dont think I can just leave the unknown as NA/NULL so I will make them their own category
#     summed$dge_block1[summed$parity %in% c('1','2','3','4')] <- 'Parous'
#     summed$dge_block1[summed$parity %in% c('0')] <- 'Nulliparous'
#     summed$dge_block2 <- summed$patient_age
#   } else if (block_var == 'none') {
#     sce$dge_block <- 'none'
#   } else {
#     print(paste0('Incorrect input: block_var = ', block_var))
#     return()
#   }
#   
#   print('plotting:')
#   print(test_var)
#   print(block_var)
#   
#   volcano_plot1 <- doVolcano(DEGs = DEGs,
#                              test_var = test_var,
#                              block_var = block_var) 
#   
#   volcano_pwd = paste0(prefix, '/dge_plotting/', test_var)
#   dir.create(volcano_pwd, showWarnings = FALSE, recursive = TRUE)
#   pdf(paste0(volcano_pwd, '/Volcano-', test_var, '-block_', block_var, '.pdf'), width=12,height = 8)
#   print(volcano_plot1)
#   dev.off()
#   
#   #plot UMAPs for some of the top dge genes
#   num_genes=10
#   genes_up <- DEGs[DEGs$logFC>0,]$X[1:num_genes]
#   genes_down <- DEGs[DEGs$logFC<0,]$X[1:num_genes]
#   
#   dir.create(paste0(prefix, '/umap/', test_var, '/'), showWarnings = FALSE, recursive = TRUE)
#   for (i in 1:num_genes) {
#     png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-UpReg_', i, '_', genes_up[i], '.png'))
#     print(plotReducedDim(sce, 'UMAP', colour_by = genes_up[i], text_by='level2'))
#     dev.off()
#     png(paste0(prefix, '/umap/', test_var, '/', 'Block_', block_var, '-DownReg_', i, '_', genes_down[i], '.png'))
#     print(plotReducedDim(sce, 'UMAP', colour_by = genes_down[i], text_by='level2'))
#     dev.off()
#   }
# }






