##This script joins together the three sce subsets (EPI, STR, IMM)
# C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/code/single_cell_gene_expression/cellchat/scvi_new/prepare_sce.R

##Libraries
library(scran)
library(scater)

library(plyr)
library(dplyr)




##Analysis

#load data
sce_epi <- readRDS('/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi/output/input/input_sce.rds')
sce_str <- readRDS('/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/str/output/input/input_sce.rds')
sce_imm <- readRDS('/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/imm/output/input/input_sce.rds')

#fix metadata
sce_imm$level0 <- 'Immune'
colData(sce_imm) <- colData(sce_imm)[, colnames(colData(sce_epi))]

#same genes
universe = intersect(rownames(sce_epi), rownames(sce_str))
universe = intersect(universe, rownames(sce_imm))

sce_list = list('epi' = sce_epi, 'str' = sce_str, 'imm' = sce_imm)

#More chages to make these compatible
for (sce_name in names(sce_list)){
  print(sce_name)
  sce <- sce_list[[sce_name]]
  rowData(sce) <- rowData(sce)[, c('X', 'gene_ids', 'feature_types', 'genome')]
  sce <- sce[universe,]
  sce_list[[sce_name]] <- sce
  
  #check
  print(sce)
}


sce_all <- do.call(cbind, sce_list)

dir.create('/nfs/research/marioni/areed/projects/hbca/cellchat/2022-04-05/scvi_new/overview_all/output/input/', recursive = T, showWarnings = F)
saveRDS(sce_all, '/nfs/research/marioni/areed/projects/hbca/cellchat/2022-04-05/scvi_new/overview_all/output/input/HBCA_postlabelling.rds')

