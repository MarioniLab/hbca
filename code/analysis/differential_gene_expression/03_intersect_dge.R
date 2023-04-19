#! Look at overall DGE gene intersections.


suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(edgeR))
suppressMessages(library(EnhancedVolcano))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))


# Example inputs

# #LP/HS/BSL
# opt = list()
# opt$input_dge_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/'
# opt$output_pwd = '/nfs/research/marioni/areed/projects/hbca/dge/2022-04-05/scvi_new/intersection/output/'




## options
option_list = list(make_option(c('--input_dge_pwd'),
                               type='character',
                               help='Pathway to directory used for saving dge results.'),
                   make_option(c('--output_pwd'),
                               type='character',
                               help='Pathway to directory used for saving output plots.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



#Analysis
block_var = 'parity_age'

lp_dge_br1 = read.csv(paste0(opt$input_dge_pwd, 'lp/output/dge_testing/WT_BRCA1PM/dge-block_', block_var, '-all.csv'))
lp_dge_br2 = read.csv(paste0(opt$input_dge_pwd, 'lp/output/dge_testing/WT_BRCA2PM/dge-block_', block_var, '-all.csv'))

hs_dge_br1 = read.csv(paste0(opt$input_dge_pwd, 'hs/output/dge_testing/WT_BRCA1PM/dge-block_', block_var, '-all.csv'))
hs_dge_br2 = read.csv(paste0(opt$input_dge_pwd, 'hs/output/dge_testing/WT_BRCA2PM/dge-block_', block_var, '-all.csv'))

bsl_dge_br1 = read.csv(paste0(opt$input_dge_pwd, 'bsl/output/dge_testing/WT_BRCA1PM/dge-block_', block_var, '-all.csv'))
bsl_dge_br2 = read.csv(paste0(opt$input_dge_pwd, 'bsl/output/dge_testing/WT_BRCA2PM/dge-block_', block_var, '-all.csv'))

dge_list = list('lp_br1' = lp_dge_br1, 'lp_br2' = lp_dge_br2, 
                'hs_br1' = hs_dge_br1, 'hs_br2' = hs_dge_br2, 
                'bsl_br1' = bsl_dge_br1, 'bsl_br2' = bsl_dge_br2)
upreg_list = list()
for (dge_name in names(dge_list)){
  dge <- dge_list[[dge_name]]
  upreg_list[[dge_name]] <- dge$X[(dge$logFC > 0) & (dge$FDR < 0.05)]
  print(length(dge$X[(dge$logFC > 1) & (dge$FDR < 0.05)]))
}  

lp_up <- intersect(upreg_list[['lp_br1']], upreg_list[['lp_br2']])
hs_up <- intersect(upreg_list[['hs_br1']], upreg_list[['hs_br2']])
bsl_up <- intersect(upreg_list[['bsl_br1']], upreg_list[['bsl_br2']])

br1_up <- intersect(intersect(upreg_list[['lp_br1']], upreg_list[['hs_br1']]), upreg_list[['bsl_br1']])
br2_up <- intersect(upreg_list[['lp_br2']], upreg_list[['hs_br2']])

all_up <- intersect(lp_up, hs_up)
all_up <- intersect(all_up, br1_up)
all_up <- intersect(all_up, br2_up)
















