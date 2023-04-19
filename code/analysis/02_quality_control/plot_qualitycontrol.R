#!/usr/bin/env Rscript
# Plot QC graphs for an individual sample.

# input:
#  - run_id  (WHY?? is this just the date?)
#  - sample_id
#  - pD_file (colData saved by qualitycontrol.R)
#  - random_seed
#  - out_QC_plots

# #Example inputs
# opt <- list()
# opt$sample_id <- 'SLX-19902-20449_SIGAA4'
# 
# #out_dir <- '/nfs/research/marioni/areed/projects/hbca/tmp_analysis/qualitycontrol'
# #out_dir <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-11/20K/output'
# out_dir <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/tmp_analysis/qualitycontrol'
# 
# opt$pD_file <- paste0(out_dir, '/', 'pD', '/', opt$sample_id, '.', 'pD.csv')
# opt$sce.alldrops_file <- paste0(out_dir, '/', 'sce.alldrops/', opt$sample_id, '.sce.alldrops.rds')
# opt$random_seed <- 42
# opt$bins <- 20
# opt$lower <- 100

# opt$out_QC_plots <- paste(out_dir, 'plots', paste0(sample_id, '.qualitycontrol.pdf'), sep='/')
# opt$out_drops_plots <- paste(out_dir, 'plots', paste0(sample_id, '.alldrops.pdf'), sep='/')



# TODO: try understand the purpose of the histogram plots and specifically the need to cut <= lower. 

#Load the required packages
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DropletUtils))

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(tidyverse))

#Options
option_list = list(make_option(c('--sample_id'),
                               type='character',
                               help='Unique identifier for the sample.'),
                  make_option(c('--pD_file'),
                              type='character',
                              help='Path to colData (pD) directory.'),
                  make_option(c('--sce.alldrops_file'),
                              type='character',
                              help='Path to alldrops colData (pD.alldrops) directory.'),
                   make_option(c('--random_seed'),
                               type='integer',
                               default=42,
                               help='Random seed to set for analysis.'),
                  make_option(c('--bins'),
                              type='integer',
                              default=20,
                              help='Number of bins used for emptyDrops histogram plots.'),
                  make_option(c('--lower'),
                              type='integer',
                              default=100,
                              help='Lower bound on the total UMI count used in emptyDrops.'),
                   make_option(c('--out_QC_plots'),
                               type = 'character',
                               help = 'Output path for quality control plots.'),
                  make_option(c('--out_drops_plots'),
                              type = 'character',
                              help = 'Output path for emptyDrops quality control plots.'))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Functions
plot_emptyDrops <- function(sce.alldrops, sample_id, out_drops_plots, random_seed, bins, lower){
  #Plot some quality control plots for the emptydroplets analysis
  #Some of the plots are based off of those used by Karsten for Tumourigenesis2021 paper.
  
  pD.alldrops <- as.data.frame(colData(sce.alldrops))
  #Shuffle
  set.seed(random_seed)
  pD.alldrops <- pD.alldrops[sample(nrow(pD.alldrops)),]
  
  bcrank <- barcodeRanks(counts(sce.alldrops))
  #only showing unique points for plotting speed.
  uniq <- !duplicated(bcrank$rank)
  bcdf <- data.frame('rank' = bcrank$rank[uniq],
                     'total' = bcrank$total[uniq])
  
  bcplot <- ggplot(bcdf, aes(x=rank, y=total)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(yintercept=metadata(bcrank)$inflection,
               lty='dashed', colour='blue') +
    geom_hline(yintercept=metadata(bcrank)$knee,
               lty='dashed') + 
    ggtitle(paste0(sample_id, ' Barcode Ranks')) + 
    xlab('Barcode Rank') +
    ylab('Total UMIs')
  
  #Make Pval histogram for droplets with barcodes < lower and check for uniform 
  #distribution to check model is suitable - no overdispersion (emptyDrops Docs).
  pD.alldrops <- pD.alldrops[pD.alldrops$total <= lower & pD.alldrops$total > 0,]

  #PValue histogram
  PVal_hist <- ggplot(pD.alldrops, aes(x=drops_PValue)) +
    geom_histogram(bins=bins) +
    xlab('PValue') +
    ggtitle(paste0(sample_id, ' PValue Histogram'))
  
  #Save plots
  dir.create(dirname(out_drops_plots), showWarnings = FALSE, recursive = TRUE)
  
  pdf(file = out_drops_plots)
  print(bcplot)
  print(PVal_hist)
  dev.off()
  
}

plot_QC <- function(pD, sample_id, out_QC_plots, random_seed) {
  #Plot some of the quality control metrics for a single sample
  
  # I will also shuffle the rows of pD and columns of m for plotting purpose and to prevent any funny effects from the ordering
  set.seed(random_seed)
  pD <- pD[sample(nrow(pD)),]
  
  #Violin plot of UMI counts/genes/mito% per cell for the sample
  umi_violin <- ggplot(pD, aes(x=sample_id,y=sum, colour = pass_all)) +
    geom_quasirandom() + 
    scale_y_log10() +
    ylab("umi") +
    ggtitle("Total UMI Counts") +
    theme_bw()
  detected_violin <- ggplot(pD, aes(x=sample_id,y=detected, colour = pass_all)) +
    geom_quasirandom() + 
    scale_y_log10() +
    ylab("detected") +
    ggtitle("Detected Genes") +
    theme_bw()
  mt_violin <- ggplot(pD, aes(x=sample_id,y=subsets_Mito_percent, colour = pass_all)) +
    geom_quasirandom() + 
    scale_y_log10() +
    ylab("subsets_Mito_percent") +
    ggtitle("Mito Percent") +
    theme_bw()
  
  #set up clean pD and cell counts before and after.
  pD.pass <- pD[pD$pass_all,]
  num_pass <- dim(pD.pass)[1]
  num_cells <- dim(pD)[1]
  
  
  #Scatter plots for correlation between metrics
  #UMI count and gene count per cell
  umi_detected_scatter <- ggplot(pD, aes(x=sum, y=detected, color=subsets_Mito_percent)) +
    geom_point(shape=1) +
    scale_color_gradient(low="gray47", high="red") +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle(paste0("Pre-QC %Mito vs. UMIs (ncells = ", as.character(num_cells), ')'))
  
  #UMI and gene count per cell with clean pD
  umi_detected_scatter_clean <- ggplot(pD.pass, aes(x=sum, y=detected, color=subsets_Mito_percent)) +
    geom_point(shape=1) +
    scale_color_gradient(low="gray47", high="red") +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle(paste0("Post-QC %Mito vs. UMIs (ncells = ", as.character(num_pass), ')'))
  
  
  #Lastly as a sanity check lets do DGE between lost cells and kept. 
  #Shouldn't see marker genes upregulated in lost.
  # NOT DONE CURRENTLY BUT MAYBE LATER - NOT SURE IF ITS HELPFUL?
  
  #Save
  dir.create(dirname(out_QC_plots), showWarnings = FALSE, recursive = TRUE)
  
  pdf(file = paste0(out_QC_plots))
  print(umi_violin)
  print(detected_violin)
  print(mt_violin)
  print(umi_detected_scatter)
  print(umi_detected_scatter_clean)
  dev.off()
  
  #Removed plots for now
  #mito_detected_scatter
  #mito_umi_scatter
  #mito_detected_scatter %+% pD.pass
  #mito_umi_scatter %+% pD.pass
}


#Complete plots
sce.alldrops <- readRDS(file = opt$sce.alldrops_file)
plot_emptyDrops(sce.alldrops <- sce.alldrops, 
                sample_id = opt$sample_id, 
                out_drops_plots = opt$out_drops_plots, 
                random_seed = opt$random_seed, 
                bins = opt$bins, 
                lower = opt$lower)

pD <- read.csv(opt$pD_file)
plot_QC(pD = pD, 
        sample_id = opt$sample_id, 
        out_QC_plots = opt$out_QC_plots, 
        random_seed = opt$random_seed)




















