#! run inferCNV 


# load packages
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(infercnv))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(tidyr))

# #Example options
opt <- list()
opt$output_subdir <- 'all_vs_individual_ddc_and_control'
opt$analysis_mode <- 'subclusters'
opt$BayesMaxPNormal <- 0.5
opt$output_dir <- paste0('/nfs/research/marioni/areed/projects/hbca/infercnv/2023-06-21/scvi/output/', opt$output_subdir, '/')
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

#Load data
sce <- readRDS('/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/all/sce_all_annotated.rds')

sce <- sce[, !(sce$level2 %in% c('Doublet', 'stripped_nuclei', ''))]
sce$test_donors <- sce$patientID %in% c('2975PM', 'LS15-02672N', 'LS13-3045N', 'LS15-2510N')

set.seed(1) #reproducibility when sampling.
cellIDs_mammoplasty_all <- sce$cellID[sce$tissue_condition == 'Mammoplasty WT' & !(sce$test_donors) & !(sce$level2 %in% c('DDC1', 'DDC2'))]
cellIDs_mammoplasty <- cellIDs_mammoplasty_all[sample(length(cellIDs_mammoplasty_all), 20000)]


#check sampled cell type spread:
print(table(sce$level1[sce$cellID %in% cellIDs_mammoplasty]))

#Have the option to take only ddc epithelium, but I think it makes sense to include the ddc stroma etc as a reference of this donors general expression patterns over the chromosomes.
#sce_ddc <- sce[, (sce$test_donors & sce$level0 %in% 'Epithelial') | sce$cellID %in% cellIDs_mammoplasty] 
sce_ddc <- sce[, sce$test_donors | sce$cellID %in% cellIDs_mammoplasty] 

gene_order <- read.table('/nfs/research/marioni/areed/projects/hbca/misc/genes_pos_refdata-gex-GRCh38-2020-A.txt', sep='\t', header=FALSE, row.names=1)

level1_celltypes <- c("Luminal adaptive secretory precurser", "Luminal hormone sensing", "Basal-myoepithelial", 
                      "Fibroblast", "Vascular endothelial", "Lymphatic endothelial", "Perivascular", 
                      "Lymphoid", "Myeloid")

#run infercnv
patient_list <- c('2975PM', 'LS15-02672N', 'LS13-3045N', 'LS15-2510N')
for (patient in patient_list) {
  print(patient)
  
  sce_sub <- sce_ddc[, sce_ddc$patientID %in% c(patient) | sce_ddc$cellID %in% cellIDs_mammoplasty]
  
  if (patient == 'LS15-02672N'){
    sce_sub <- sce_sub[, sce_sub$level2 != 'DDC1'] #there is one classified DDC1 cell which causes errors.
    ref_group_names <- c(level1_celltypes) 
    #ref_group_names <- c(level1_celltypes, paste0(patient, '.', level1_celltypes[c(4:7,9)])) 
  }
  if (patient == '2975PM'){
    sce_sub <- sce_sub[, sce_sub$level2 != 'DDC2'] #there are two classified DDC2 cells which causes errors.
    ref_group_names <- c(level1_celltypes) 
    #ref_group_names <- c(level1_celltypes, paste0(patient, '.', level1_celltypes[4:9]))
  }
  if (patient == 'LS13-3045N'){
    sce_sub <- sce_sub[, !(sce_sub$level2 %in% c('DDC1', 'DDC2'))]
    ref_group_names <- c(level1_celltypes) #, paste0(patient, '.', level1_celltypes[c(4:9)]))
  }
  if (patient == 'LS15-2510N'){
    sce_sub <- sce_sub[, !(sce_sub$level2 %in% c('DDC1', 'DDC2'))]
    ref_group_names <- c(level1_celltypes) #, paste0(patient, '.', level1_celltypes[c(4:9)]))
  }
  
  ##Prepare for inferCNV
  #I want to see the differences between LP, BSL, HS (as controls) and then DDC1/2
  sce_sub$test <- sce_sub$level1
  sce_sub$test[sce_sub$level2 %in% c('DDC1', 'DDC2')] <- sce_sub$level2[sce_sub$level2 %in% c('DDC1', 'DDC2')]
  sce_sub$tested_cells_ddc <- (sce_sub$level0 == 'Epithelial') & (sce_sub$patientID %in% c('2975PM', 'LS15-02672N'))
  sce_sub$tested_cells_control <- (sce_sub$level0 == 'Epithelial') & (sce_sub$patientID %in% c('LS13-3045N', 'LS15-2510N'))
  sce_sub$tested_cells <- sce_sub$tested_cells_ddc | sce_sub$tested_cells_control
  sce_sub$test[sce_sub$tested_cells_ddc] <- paste0('Test-', sce_sub$test[sce_sub$tested_cells_ddc])
  sce_sub$test[sce_sub$tested_cells_control] <- paste0('Control-', sce_sub$test[sce_sub$tested_cells_control])
  
  #for ease of visualisation subsample the maximum number of cells for each tested cell type to maximum 700 cells ~ number of DDC1.
  sce_sub$keep <- TRUE
  set.seed(1)
  for (test_celltype in unique(sce_sub$test[sce_sub$tested_cells])) {
    number_of_test_cells <- sum(sce_sub$test == test_celltype)
    if (number_of_test_cells > 700){
      sample_700 <- sample(1:number_of_test_cells, 700)
      sce_sub$keep[sce_sub$test == test_celltype][!(1:number_of_test_cells %in% sample_700)] <- FALSE
    }
  }
  
  #remove any mammoplasty sampled cells classified as DDC1/2. 
  #These should be very rare but can cause errors (there was one DDC1 cell selected in 2975PM donor).
  sce_sub <- sce_sub[, !(!sce_sub$tested_cells_ddc & (sce_sub$test %in% c('DDC1', 'DDC2')))]
  
  #check this worked correctly
  print('Table 1:')
  print(table(sce_sub$tested_cells, sce_sub$test))
  print('Table 2:')
  print(table(sce_sub$keep, sce_sub$test))
  
  sce_sub <- sce_sub[, sce_sub$keep]
  print('Table 3:')
  print(table(sce_sub$tested_cells, sce_sub$test))
  
  annotations_df <- data.frame(row.names = colnames(sce_sub), 'V2' = sce_sub$test)
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts(sce_sub)),
                                       annotations_file=annotations_df,
                                       delim="\t",
                                       gene_order_file=gene_order,
                                       ref_group_names=ref_group_names,
                                       chr_exclude = c('GL000009.2', 'GL000194.1', 'GL000195.1', 
                                                       'GL000213.1', 'GL000218.1', 'GL000219.1', 
                                                       'KI270711.1', 'KI270713.1', 'KI270721.1', 
                                                       'KI270726.1', 'KI270727.1', 'KI270728.1', 
                                                       'KI270731.1', 'KI270734.1', 'chrM'))
  
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir=paste0(opt$output_dir, patient),  # dir is auto-created for storing outputs
                                cluster_by_groups=T,   # cluster
                                denoise=T,
                                HMM=T,
                                BayesMaxPNormal = opt$BayesMaxPNormal,
                                analysis_mode = opt$analysis_mode)
  
  
  if (opt$analysis_mode == 'samples') {
    genes_file <- paste0(opt$output_dir, patient,
                         '/HMM_CNV_predictions.HMMi6.hmm_mode-',
                         opt$analysis_mode,
                         '.Pnorm_',
                         opt$BayesMaxPNormal, 
                         '.pred_cnv_genes.dat')
    genes_out_file <- paste0(opt$output_dir, patient,
                             '/HMM_CNV_predictions.HMMi6.hmm_mode-',
                             opt$analysis_mode,
                             '.Pnorm_',
                             opt$BayesMaxPNormal, 
                             '.pred_cnv_gene_symbols.dat')
  } else {
    genes_file <- paste0(opt$output_dir, patient,
                         '/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-',
                         opt$analysis_mode,
                         '.Pnorm_',
                         opt$BayesMaxPNormal, 
                         '.pred_cnv_genes.dat')
    genes_out_file <- paste0(opt$output_dir, patient,
                             '/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-',
                             opt$analysis_mode,
                             '.Pnorm_',
                             opt$BayesMaxPNormal, 
                             '.pred_cnv_gene_symbols.dat')
  }
  
  if (file.exists(genes_file)) {
    CNV_gene_by_chr <- read.table(genes_file, 
                                  header=TRUE, sep='\t')
    
    gene_symbols <- rowData(sce_sub)[CNV_gene_by_chr$gene,]$X
    
    CNV_gene_by_chr_symbols <- CNV_gene_by_chr
    CNV_gene_by_chr_symbols$gene_symbol <- gene_symbols
    
    write.table(CNV_gene_by_chr_symbols,
                file = genes_out_file)
  }
  
  
}

