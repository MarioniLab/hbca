#!/usr/bin/env Rscript
# Complete emptyDrops and general quality control on single cell RNA-seq data

# input:
#  - run_id
#  - sample_id
#  - out_dir
#  - barcode_file
#  - feature_file
#  - count_matrix_file
#  - mito_file
#  - tf_file
#  - lower (bound for emptydrops)
#  - niters (for emptydrops pvals)
#  - sig (for FDR of emptydrops)
#  - umi_nmads
#  - detected_nmads
#  - mt_nmads
#  - min_count_depth
#  - min_num_genes
#  - random_seed




# #Example inputs
# opt <- list()
# opt$sample_id <- 'SLX-19902-20449_SIGAA4'
# 
# #out_dir <- '/nfs/research/marioni/areed/projects/hbca/tmp_analysis/qualitycontrol'
# out_dir <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/tmp_analysis/qualitycontrol/2022-03-11'
# run_dir <- '/nfs/research/marioni/asteif/projects/hbca/downsample/2021-07-10/20K/output/downsample'
# #run_dir <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/downsample/2021-07-10/20K/output/downsample'
# 
# opt$barcode_file <- paste0(run_dir, '/', opt$sample_id, '.downsample_barcodes_raw.tsv.gz')
# opt$feature_file <- paste0(run_dir, '/', opt$sample_id, '.downsample_features.tsv.gz')
# opt$count_matrix <- paste0(run_dir, '/', opt$sample_id, '.downsample_matrix_raw.mtx.gz')
# 
# #opt$mito_file <- '/nfs/research/marioni/areed/projects/hbca/code/single_cell_gene_expression/qualitycontrol/config/MitoGenes.txt'
# #opt$tf_file <- '/nfs/research/marioni/areed/projects/hbca/code/single_cell_gene_expression/qualitycontrol/config/TFcheckpoint_WithENSID.txt'
# opt$mito_file <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/qualitycontrol/config/MitoGenes.txt'
# opt$tf_file <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/qualitycontrol/config/TFcheckpoint_WithENSID.txt'
# opt$metadata_file <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/metadata/2021-08-10_master_table.csv'
# #opt$metadata_file <- "/nfs/research/marioni/areed/projects/hbca/metadata/2021-08-10_master_table.csv"
# 
# opt$lower <- 100
# opt$niters <- 10000
# opt$sig <- 0.001
# opt$umi_nmads <- 3
# opt$detected_nmads <-3
# opt$mt_nmads <- 3
# opt$min_count_depth <- 600
# opt$min_num_genes <- 0
# opt$random_seed <- 42
# 
# opt$out_sce <- paste0(out_dir, '/', 'sce', '/', opt$sample_id, '.', 'sce.rds')
# opt$out_pD <- paste0(out_dir, '/', 'pD', '/', opt$sample_id, '.', 'pD.rds')
# 
# opt$out_sce.alldrops <- paste0(out_dir, '/', 'sce.alldrops', '/', opt$sample_id, '.', 'sce.alldrops.rds')
# opt$out_pD.alldrops <- paste0(out_dir, '/', 'pD.alldrops', '/', opt$sample_id, '.', 'pD.alldrops.rds')
# 
# opt$features_out_pw <- paste0(out_dir, '/', 'final_qc_data', '/', opt$sample_id, '.downsample_features.tsv.gz')
# opt$barcodes_out_pw <- paste0(out_dir, '/', 'final_qc_data', '/', opt$sample_id, '.downsample_barcodes_raw.tsv.gz')
# opt$counts_out_pw <- paste0(out_dir, '/', 'final_qc_data', '/', opt$sample_id, '.downsample_matrix_raw.mtx.gz')



##TODO: 
##remove genes that are not expressed - from joined sce?. Think about latter in the normalisation.
##Ask whether I should use min_count_depth & min_num_genes as well as the isOutlier MADs.
##Optimise the sig value




#Load the required packages
suppressMessages(library(optparse))
suppressMessages(library(scran))
suppressMessages(library(DropletUtils))
suppressMessages(library(Matrix))
suppressMessages(library(scater))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(tidyverse))
suppressMessages(library(readr))
# 

#options
option_list = list(make_option(c('--sample_id'),
                               type='character',
                               help='Unique identifier for the sample.'),
                   make_option(c('--out_dir'),
                               type='character',
                               help='Path to output directory'),
                   make_option(c('--barcode_file'),
                               type='character',
                               help='Path to barcode file.'),
                   make_option(c('--feature_file'),
                               type='character',
                               help='Path to feature file.'),
                   make_option(c('--count_matrix'),
                               type='character',
                               help='Path to count matrix.'),
                   make_option(c('--mito_file'),
                               type='character',
                               help='Path to mitochondrial genes file.'),
                   make_option(c('--tf_file'),
                               type='character',
                               help='Path to transcription factor list file.'),
                   make_option(c('--metadata_file'),
                                type='character',
                                help='Path to metadata .csv file.'),
                   make_option(c('--lower'),
                               type='integer',
                               default=100,
                               help='Lower bound on the total UMI count for emptyDrops.'),
                   make_option(c('--niters'),
                               type='integer',
                               default=10000,
                               help='The number of iterations for the Monte Carlo p-value calculations in emptyDrops.'),
                   make_option(c('--sig'),
                               type='double',
                               default = 0.001,
                               help='Threshold for FDR to qualify as a non-empty droplet in emptyDrops.'),
                   make_option(c('--umi_nmads'),
                               type='double',
                               default=3,
                               help='Threshold for umi MADs.'),
                   make_option(c('--detected_nmads'),
                               type='double',
                               default=3,
                               help='Threshold for detected genes MADs.'),
                   make_option(c('--mt_nmads'),
                               type='double',
                               default=3,
                               help='Threshold for mitochondrial percentage MADs.'),
                   make_option(c('--min_count_depth'),
                               type='integer',
                               default=500,
                               help='Minimum count depth per barcode.'),
                   make_option(c('--min_num_genes'),
                               type='integer',
                               default=300,
                               help='Minimum number of genes per barcode.'),
                   make_option(c('--max_mito'),
                               type='double',
                               default=0.15,
                               help='Maximum mitochondrial percentage allowed through QC.'),
                   make_option(c('--random_seed'),
                               type='integer',
                               default=42,
                               help='Random seed for to set for analysis.'),
                   make_option(c('--out_sce'),
                               type='character',
                               help='Path to sce output.'),
                   make_option(c('--out_pD'),
                               type='character',
                               help='Path to pD output.'),
                   make_option(c('--out_sce.alldrops'),
                               type='character',
                               help='Path to sce.alldrops output.'),
                   make_option(c('--out_pD.alldrops'),
                               type='character',
                               help='Path to pD.alldrops output.'),
                   make_option(c('--features_out_pw'),
                               type='character',
                               help='Path to features output (.tsv.gz).'),
                   make_option(c('--barcodes_out_pw'),
                               type='character',
                               help='Path to barcodes output (.tsv.gz).'),
                   make_option(c('--counts_out_pw'),
                               type='character',
                               help='Path to counts output (.mtx.gz).'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#Functions
read_count_data <- function(barcode_file, feature_file, count_matrix,
                            mito_file, tf_file, metadata_file, sample_id){
  # Read the count matrix and return a single cell experiment object
  counts <- readMM(count_matrix)

  features <- read.delim(feature_file,
                         header = FALSE,
                         stringsAsFactors = FALSE)
  colnames(features) <- c('gene_id', 'gene_name', 'gene_type')

  barcodes <- read.delim(barcode_file,
                         header = FALSE,
                         stringsAsFactors = FALSE)
  colnames(barcodes) <- c('barcode')

  #Also add in some rowdata about ensemblID, TranscriptionFactors and Mitochondrial genes.
  rwData <- features

  #read mito and TF genes and add to rowData
  mitoGenes <- read.table(mito_file)
  rwData$mito <- rwData$gene_id %in% mitoGenes$V1

  TFGenes <- read.table(tf_file, sep="\t")
  TFEnsemblID <- TFGenes$V2 #select Human Data EnsemblID's
  rwData$tf <- rwData$gene_id %in% TFEnsemblID

  #Fix colnames
  rownames(rwData) <- uniquifyFeatureNames(rwData$gene_id, rwData$gene_name)

  colnames(counts) <- paste0(sample_id, '_', barcodes$barcode)
  rownames(counts) <- rownames(rwData)

  #Create colData df and add meta data
  mast_cols <- c('sample_id', 'patient_id', 'sample_type', 'replicate',
                 'processing_date', 'patient_age', 'age_range', 'parity', 'tissue_condition')
  metadata <- as.data.frame(read.csv(metadata_file))[, mast_cols]
  smpl_data <- metadata[metadata$sample_id == sample_id,]

  clData <- data.frame(barcode = barcodes$barcode)
  rownames(clData) <- paste0(sample_id, '_', barcodes$barcode)
  clData$sample_id <- sample_id
  #now add metadata from mastertable
  clData$patient_id <- smpl_data$patient_id
  clData$sample_type <- smpl_data$sample_type
  clData$replicate <- as.numeric(smpl_data$replicate)
  clData$processing_date <- smpl_data$processing_date
  clData$patient_age <- as.numeric(smpl_data$patient_age)
  clData$age_range <- smpl_data$age_range
  clData$parity <- smpl_data$parity
  clData$tissue_condition <- smpl_data$tissue_condition

  #Create sce
  sce <- SingleCellExperiment(assays=list(counts=counts), rowData = rwData, colData = clData)

  return(sce)
}

complete_emptyDrops <- function(sce, lower=100, niters=10000, sig=0.001, random_seed){
  # Use emptyDrops to determine cells from ambient RNA
  
  set.seed(random_seed)
  e.out <- emptyDrops(counts(sce), lower = lower, niters = niters)
  e.out.all <- emptyDrops(counts(sce), lower = lower, niters = niters, test.ambient = TRUE)
  is.cell <- e.out$FDR <= sig
  is.cell[is.na(is.cell)] <- FALSE #droplets with total umi < lower have NA for all e.out data
  
  #keeping only the non-empty droplets in sce.cleaned but mark the 'cells' in sce colData.
  sce.cleaned <- sce[,is.cell]
  sce$total <- e.out.all$Total
  sce$is.cell <- is.cell
  sce$limited <- e.out.all$Limited 
  sce$drops_PValue <- e.out.all$PValue
  sce$drops_FDR <- e.out.all$FDR
  
  #Limited will tell us if a lower pval could have been attained with higher niters parameter. 
  #If niters is sufficiently large this matrix should be upper triangular.
  tbl <- table(Limited=e.out$Limited, Significant=is.cell)
  
  out <- list('SCE_cleaned' = sce.cleaned,
              'SCE_alldrops' = sce,
              'Table' = tbl)
  
  return(out)
}

do_perCellQC <- function(sce, umi_nmads, detected_nmads, mt_nmads, 
                         min_count_depth, min_num_genes, max_mito){
  # Add PerCellQC metrics to sce and mark cells of insufficient quality using the relative amount of nmads for each metric
  sce <- addPerCellQC(sce, subsets=list(Mito=rowData(sce)$mito))
  
  #Removal by hard cut offs.
  # counts <- counts(sce)
  # count_depth <- colSums(counts)
  # num_genes <- colSums(counts > 0)
  # df = data.frame(count_depth=count_depth, num_genes=num_genes, subsets_Mito)
  pass_cutoff <- (sce$detected >= min_num_genes) & (sce$sum >= min_count_depth) & (sce$subsets_Mito_percent <= max_mito) #| sce$sum > 1000 )
  
  #Removal by mean absolute deviance from those that pass the cutoff.
  umi.outlier <- isOutlier(sce$sum[pass_cutoff], nmads=umi_nmads, type="lower", log=TRUE)
  d.outlier <- isOutlier(sce$detected[pass_cutoff], nmads=detected_nmads, type="lower", log=TRUE)
  mt.outlier <- isOutlier(sce$subsets_Mito_percent[pass_cutoff], nmads=mt_nmads, type="higher")
  
  sce$pass_cutoff <- pass_cutoff
  sce$pass_umi[pass_cutoff] <- !umi.outlier
  sce$pass_umi[!pass_cutoff] <- FALSE #already removed
  sce$pass_detected[pass_cutoff] <- !d.outlier
  sce$pass_detected[!pass_cutoff] <- FALSE #already removed
  sce$pass_mt[pass_cutoff] <- !mt.outlier #| sce$sum[pass_cutoff] > 1000
  sce$pass_mt[!pass_cutoff] <- FALSE #already removed
  
  sce$pass_all <- sce$pass_umi & sce$pass_detected & sce$pass_mt
  
  return(sce)
}

save_data <- function(sce.cells, sce.alldrops, out_sce, out_sce.alldrops, out_pD, 
                      out_pD.alldrops, sample_id, features_out_pw, barcodes_out_pw,
                      counts_out_pw){
  pD <- as.data.frame(colData(sce.cells)[,c('barcode', 'sample_id','sum','detected','subsets_Mito_percent','pass_all','patient_id', 
                                            'sample_type', 'replicate', 'processing_date', 'patient_age', 'age_range', 'parity', 
                                            'tissue_condition')])
  pD.alldrops <- as.data.frame(colData(sce.alldrops)[,c('barcode', 'sample_id','is.cell', 'patient_id', 'sample_type', 'replicate',
                                                        'processing_date', 'patient_age', 'age_range', 'parity', 'tissue_condition')])
  
  #Create directories in-case they do not exist.
  dir.create(dirname(out_sce), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(out_sce.alldrops), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(out_pD), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(out_pD.alldrops), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(features_out_pw), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(barcodes_out_pw), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(counts_out_pw), showWarnings = FALSE, recursive = TRUE)
  
  saveRDS(sce.cells, file = out_sce)
  saveRDS(sce.alldrops, file = out_sce.alldrops)
  write.csv(pD, file = out_pD)
  write.csv(pD.alldrops, file = out_pD.alldrops)
  
  #Save the quality controlled data (remove cells of bad quality) in the same format as CellRanger output.
  sce.qc <- sce.cells[,sce.cells$pass_all]
  
  features_out <- as.data.frame(rowData(sce.qc)[,c('gene_id', 'gene_name', 'gene_type')])
  rownames(features_out) <- NULL
  barcodes_out <- as.data.frame(colData(sce.qc)[,c('barcode')])
  counts_out <- counts(sce.qc) #as.matrix(counts(sce.qc))
  rownames(counts_out) <- NULL
  colnames(counts_out) <- NULL
  
  write_tsv(features_out, file = features_out_pw, col_names = FALSE) #tsv.gz
  write_tsv(barcodes_out, file = barcodes_out_pw, col_names = FALSE) #tsv.gz
  writeMM(counts_out, file = counts_out_pw) #mtx.gz
}






#Complete the Analysis
sce <- read_count_data(barcode_file = opt$barcode_file,
                       feature_file = opt$feature_file,
                       count_matrix = opt$count_matrix,
                       mito_file = opt$mito_file,
                       tf_file = opt$tf_file,
                       metadata_file = opt$metadata_file,
                       sample_id = opt$sample_id)
sce

out <- complete_emptyDrops(sce = sce, 
                           lower = opt$lower,
                           niters = opt$niters,
                           sig = opt$sig,
                           random_seed = opt$random_seed)
sce.cells <- out[[1]]
sce.alldrops <- out[[2]]
tbl <- out[[3]]

print('done emptyDrops')

sce.QC <- do_perCellQC(sce = sce.cells,
                       umi_nmads = opt$umi_nmads,
                       detected_nmads = opt$detected_nmads,
                       mt_nmads = opt$mt_nmads,
                       min_count_depth = opt$min_count_depth,
                       min_num_genes = opt$min_num_genes,
                       max_mito = opt$max_mito)

print('done QC')

save_data(sce.cells = sce.QC, 
          sce.alldrops = sce.alldrops, 
          out_sce = opt$out_sce, 
          out_sce.alldrops = opt$out_sce.alldrops,
          out_pD = opt$out_pD, 
          out_pD.alldrops = opt$out_pD.alldrops,
          sample_id = opt$sample_id,
          features_out_pw = opt$features_out_pw, 
          barcodes_out_pw = opt$barcodes_out_pw,
          counts_out_pw = opt$counts_out_pw)

print('done saving')







