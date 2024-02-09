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
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/lp/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/lp/output/'
# opt$celltype_subset = 'LASP1,LASP2,LASP3,LASP4,LASP5'

#LHS
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/hs/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/hs/output/'
# opt$celltype_subset = 'LHS1,LHS2,LHS3'

#BSL
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/bsl/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/bsl/output/'
# opt$celltype_subset = 'BMYO1,BMYO2'

#Macro
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/imm/sce_imm_annotated.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/macro/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/macro/output/'
# opt$celltype_subset = 'Macro'

#Macro
# opt = list()
# opt$input_sce_path = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/imm/sce_imm_annotated.rds'
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/dc/output/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2023-06-21/scvi/dc/output/'
# opt$celltype_subset = 'DC'


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
    FCcutoff <- 0.01 #approximately equivalent to 1 logFC over lifetime
    x_buffer <- 0.1
  } else if (test_var == 'Body_mass_index') {
    FCcutoff <- 0.1 #equivalent to 1 logFC over 10 BMI shift
    x_buffer <- 0.15
  } else {
    x_buffer <- 1.5
  }
  if (no_label) {
    Volcano <- EnhancedVolcano(DEGs,
                               lab = NA,
                               x = 'logFC',
                               y = 'FDR',
                               xlim = c(min(DEGs[['logFC']], na.rm = TRUE) - x_buffer, 
                                        max(DEGs[['logFC']], na.rm = TRUE) + x_buffer),
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
                               xlim = c(min(DEGs[['logFC']], na.rm = TRUE) - x_buffer, 
                                        max(DEGs[['logFC']], na.rm = TRUE) + x_buffer),
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

}
 


# Analysis 

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
test_var_list <- list('WT_BRCA1PM','WT_BRCA2PM')
block_var_list <- list('parity_age')

#complete tests
for (test_var in test_var_list){
  for (block_var in block_var_list){
    print('New round of plot:')
    print(test_var)
    print(block_var)
    DEGs <- read.csv(paste0(opt$input_dge_pwd, '/dge_testing/', test_var, '/dge-block_', block_var, '-all.csv'))
    setup_and_plot(sce = sce, DEGs = DEGs, test_var = test_var, block_var = block_var,
                   prefix = opt$output_pwd)
  }
}



