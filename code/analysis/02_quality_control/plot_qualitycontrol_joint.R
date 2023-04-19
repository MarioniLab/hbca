#!/usr/bin/env Rscript
# Plot QC graphs over all samples.

# input:
#  - run_id  (WHY?? is this just the date?)
#  - samples_list.
#  - pD_pwd (Path to directory storing colData saved by qualitycontrol.R)
#  - random_seed
#  - out_QC_plots_joint

# #Example inputs
# opt <- list()
# opt$mastertable_pw <- '/nfs/research/marioni/areed/projects/hbca/metadata/2021-08-10_master_table.csv'

# 
# out_dir <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-11/20K/output'
# #out_dir <- '/home/austin/OneDrive/WTKLAB/Projects/hbca/codon/tmp_analysis/qualitycontrol'
# 
# opt$pD_pwd <- paste0(out_dir, '/', 'pD')
# opt$random_seed <- 42
# opt$out_QC_plots_joint <- paste(out_dir, 'QC_plots_joint', sep='/')

##Real example inputs
# opt <- list()
# opt$mastertable_pw <- '/nfs/research/marioni/areed/projects/hbca/metadata/2021-08-10_master_table.csv'
# opt$pD_pwd <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/pD'
# opt$random_seed <- 42
# opt$out_QC_plots_joint_tissuecondition <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/plots/allsamples.qualitycontrol.tissuecondition.pdf'
# opt$out_QC_plots_joint_tissuecondition_clean <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/plots/allsamples.qualitycontrol.tissuecondition.clean.pdf'
# opt$out_QC_plots_joint_sampletype <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/plots/allsamples.qualitycontrol.sampletype.pdf'
# opt$out_QC_plots_joint_sampletype_clean <- '/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/output/plots/allsamples.qualitycontrol.sampletype.clean.pdf'
# 


#Load the required packages
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DropletUtils))

suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(tidyverse))

#Options
option_list = list(make_option(c('--mastertable_pw'),
                               type='character',
                               help='Path to master table .csv file.'),
                   make_option(c('--pD_pwd'),
                               type='character',
                               help='Path to colData (pD) directory.'),
                   make_option(c('--random_seed'),
                               type='integer',
                               default=42,
                               help='Random seed for to set for analysis.'),
                   make_option(c('--out_QC_plots_joint_tissuecondition'),
                               type = 'character',
                               help = 'Output path for joint quality control plots in .pdf format.'),
                   make_option(c('--out_QC_plots_joint_tissuecondition_clean'),
                               type = 'character',
                               help = 'Output path for joint quality control plots in .pdf format.'),
                   make_option(c('--out_QC_plots_joint_sampletype'),
                               type = 'character',
                               help = 'Output path for joint quality control plots in .pdf format.'),
                   make_option(c('--out_QC_plots_joint_sampletype_clean'),
                               type = 'character',
                               help = 'Output path for joint quality control plots in .pdf format.'))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Functions
plot_QC_joint <- function(mastertable_pw, pD_pwd, random_seed, 
                          out_QC_plots_joint_tissuecondition, 
                          out_QC_plots_joint_tissuecondition_clean, 
                          out_QC_plots_joint_sampletype, 
                          out_QC_plots_joint_sampletype_clean){
  #Plot some quality control plots over all samples.
  
  #Can only input a character string so break up string to get individual items.
  mst_tbl <- as.data.frame(read.csv(file = mastertable_pw))
  samples_list <- mst_tbl$sample_id
  #samples_list <- lapply(strsplit(as.character(samples_list), "[][']"), function(x) x[nzchar(x)])[[1]]
  
  print(head(mst_tbl))
  print(samples_list)
  
  #read in and join the data frames
  pD_list <- list()
  pD_joint <- data.frame()
  for (smpl in samples_list){
    if (file.exists(paste0(pD_pwd, '/', smpl, '.pD.csv'))){
      pD_temp <- as.data.frame(read.csv(file = paste0(pD_pwd, '/', smpl, '.pD.csv')))
      pD_joint <- rbind(pD_joint, pD_temp)
    }
  }
  
  print(head(pD_joint)) 
  print(dim(pD_joint)) #1651763      15
  
  
  #This line is only required currently as the tissue_condition colData is missing.
  #pD_joint <- merge(pD_joint, mst_tbl[,c('sample_id','tissue_condition', 'sample_type')], by='sample_id', sort=FALSE)
  
  print(head(pD_joint))
  
  #Set up colours for the tissue_conditions consistent with the colours used in 'plot_joint_normalization.py'.
  colours_tc <- c('Mammoplasty WT' = '#f1ce63', 
                  'Mastectomy BRCA1' = '#4e79a7', 
                  'Mastectomy BRCA2' = '#a0cbe8', 
                  'Mastectomy WT' = '#f28e2b',  
                  'Mastectomy unknown' = '#bab0ac', 
                  'Contralateral BRCA1' =  '#ff9d9a')
  
  #decide on colours for sample_type
  #colours_st <- c()
  
  #make a clean pD
  pD_joint_clean <- pD_joint[pD_joint$pass_all, ]
  
  #Make violins coloured by tissue_condition
  umi_violin <- ggplot(pD_joint, aes(x=sample_id,y=sum, fill=tissue_condition)) +
    geom_violin() +
    scale_y_log10() +
    ylab("umi") +
    ggtitle("Total UMI Counts by Sample") +
    scale_color_manual(values = colours_tc) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  detected_violin <- ggplot(pD_joint, aes(x=sample_id,y=detected, fill=tissue_condition)) +
    geom_violin() +
    scale_y_log10() +
    ylab("detected") +
    ggtitle("Total Genes by Sample") +
    scale_color_manual(values = colours_tc) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  mt_violin <- ggplot(pD_joint, aes(x=sample_id,y=subsets_Mito_percent, fill=tissue_condition)) +
    geom_violin() +
    scale_y_log10() +
    ylab("subsets_Mito_percent") +
    ggtitle("Mito Percent by Sample") +
    scale_color_manual(values = colours_tc) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  #Make violins coloured by sample_type
  umi_violin2 <- ggplot(pD_joint, aes(x=sample_id,y=sum, fill=sample_type)) +
    geom_violin() +
    scale_y_log10() +
    ylab("umi") +
    ggtitle("Total UMI Counts by Sample") +
    #scale_color_manual(values = colours_st) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  detected_violin2 <- ggplot(pD_joint, aes(x=sample_id,y=detected, fill=sample_type)) +
    geom_violin() +
    scale_y_log10() +
    ylab("detected") +
    ggtitle("Total Genes by Sample") +
    #scale_color_manual(values = colours_st) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  mt_violin2 <- ggplot(pD_joint, aes(x=sample_id,y=subsets_Mito_percent, fill=sample_type)) +
    geom_violin() +
    scale_y_log10() +
    ylab("subsets_Mito_percent") +
    ggtitle("Mito Percent by Sample") +
    #scale_color_manual(values = colours_st) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  #save all these plots
  pdf(file = out_QC_plots_joint_tissuecondition, width = 35, height= 7)
  print(umi_violin)
  print(detected_violin)
  print(mt_violin)
  dev.off()
  
  pdf(file = out_QC_plots_joint_tissuecondition_clean, width = 35, height= 7)
  print(umi_violin %+% pD_joint_clean)
  print(detected_violin %+% pD_joint_clean)
  print(mt_violin %+% pD_joint_clean)
  dev.off()
  
  pdf(file = out_QC_plots_joint_sampletype, width = 35, height= 7)
  print(umi_violin2)
  print(detected_violin2)
  print(mt_violin2)
  dev.off()
  
  pdf(file = out_QC_plots_joint_sampletype_clean, width = 35, height= 7)
  print(umi_violin2 %+% pD_joint_clean)
  print(detected_violin2 %+% pD_joint_clean)
  print(mt_violin2 %+% pD_joint_clean)
  dev.off()
}


#Complete plots
plot_QC_joint(mastertable_pw = opt$mastertable_pw, 
              pD_pwd = opt$pD_pwd, 
              random_seed = opt$random_seed, 
              out_QC_plots_joint_tissuecondition = opt$out_QC_plots_joint_tissuecondition,
              out_QC_plots_joint_tissuecondition_clean = opt$out_QC_plots_joint_tissuecondition_clean,
              out_QC_plots_joint_sampletype = opt$out_QC_plots_joint_sampletype,
              out_QC_plots_joint_sampletype_clean = opt$out_QC_plots_joint_sampletype_clean)
  
  
  
  
  