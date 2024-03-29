#!/usr/bin/env Rscript
# Complete overview of signalling pathways for the detailed cell types.

# #Example inputs

# #Level2
# opt <- list()
# opt$sce_pw <- '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/all/sce_all_annotated.rds'
# opt$celltype_level <- 'level2'
# opt$out_pwd <- '/nfs/research/marioni/areed/projects/hbca/cellchat/2023-06-21/scvi/brca_epistrimm/output'
# opt$workers <- 32

#Load required packages
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(CellChat))

suppressMessages(library(optparse))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(RColorBrewer))
suppressMessages(library(NMF))
suppressMessages(library(ggalluvial))
suppressMessages(library(ggplot2))

options(future.globals.maxSize= 9000*1024^2) #9000mb as 300000/32 = 9375

#options
option_list = list(make_option(c('--sce_pw'),
                               type='character',
                               help='Pathway to read in sce.'),
                   make_option(c('--celltype_level'),
                               type='character',
                               help='Cell type level to use for grouping.'),
                   make_option(c('--workers'),
                               type='integer',
                               default=32,
                               help='Number of workers to parallelize over.'),
                   make_option(c('--out_pwd'),
                               type='character',
                               help='Pathway to save cellchat results.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Functions


# Analysis


#Load data
sce <- readRDS(opt$sce_pw)

sce$level2[(sce$level1 %in% c("Lymphoid", "Myeloid"))] <- sce$level1[(sce$level1 %in% c("Lymphoid", "Myeloid"))]

print(table(sce$level2))

#fix rownames (cellchat fails without this!)
rownames(sce) <- rowData(sce)$X

#Subset data
celltypes_to_analyse <- c("LASP1", "LASP2", "LASP3", "LASP4",
                          "LHS1", "LHS2", "LHS3",
                          "BMYO1", "BMYO2", 
                          "FB1", "FB2", "FB3", "FB4",
                          "PV1", "PV2", "PV3", "PV4", "PV5",
                          "VEV", "VEC", "VEA","VEAT", 
                          "LE1", "LE2",
                          'Lymphoid', 'Myeloid')
sce <- sce[, colData(sce)[,opt$celltype_level] %in% celltypes_to_analyse]
sce$CellChatGroup <- factor(colData(sce)[,opt$celltype_level], 
                            levels = celltypes_to_analyse)

print(table(sce$CellChatGroup))
print(sum(is.na(sce$CellChatGroup)))



sce.sub.wt <- sce[, sce$tissue_condition == 'Mammoplasty WT']
sce.sub.brca1 <- sce[, sce$tissue_condition == 'Mastectomy BRCA1']
sce.sub.brca2 <- sce[, sce$tissue_condition == 'Mastectomy BRCA2']

#colours
level2_colour_dictionary = c('LASP1'= "#DA80DA", 'LASP2'= "#815481", 'LASP3'= "#C040C0", 'LASP4'= "#E1AFE1",
                             'LHS1'= "#EDABB9", 'LHS2'= "#EB5C79", 'LHS3'= "#A06A75",
                             'BMYO1'= "#EB675E", 'BMYO2'= "#A23E36",
                             'FB1'= "#DFA38A", 'FB2'= "#8C3612", 'FB3'= "#623623", 'FB4'= "#916350",
                             'PV1'= "#F8770B", 'PV2'= "#E09E3A", 'PV3'= "#CD7225", 'PV4'= "#FFC990", 'PV5'= "#AC5812",
                             'VEV'= "#FEE083", 'VEC'= "#897538", 'VEA'= "#E7B419", 'VEAT'= "#BCA048",
                             'LE1'= "#6F8BE2", 'LE2'= "#3053BC",
                             "Lymphoid" = '#AAB256', 
                             "Myeloid" = '#9FC5E8') 


#Create CellChat Object of each dataset before merging
set.seed(42)

#We saw subsetting was unneccessary and had little overall effect (so maybe remove?)
cellchat.wt <- createCellChat(sce.sub.wt, group.by = 'CellChatGroup')
cells_to_sample <- sample(1:dim(sce.sub.wt)[2], dim(sce.sub.brca1)[2])
cellchat.wt_sub1 <- createCellChat(sce.sub.wt[, cells_to_sample[order(cells_to_sample)]], group.by = 'CellChatGroup')
cells_to_sample <- sample(1:dim(sce.sub.wt)[2], dim(sce.sub.brca2)[2])
cellchat.wt_sub2 <- createCellChat(sce.sub.wt[, cells_to_sample[order(cells_to_sample)]], group.by = 'CellChatGroup')
cellchat.brca1 <- createCellChat(sce.sub.brca1, group.by = 'CellChatGroup')
cellchat.brca2 <- createCellChat(sce.sub.brca2, group.by = 'CellChatGroup')
cellchat.list <- list('wt' = cellchat.wt, 'wt_sub1' = cellchat.wt_sub1, 'wt_sub2' = cellchat.wt_sub2, 'brca1' = cellchat.brca1, 'brca2' = cellchat.brca2)

for (genotype in c('wt', 'wt_sub1', 'wt_sub2', 'brca1', 'brca2')){
  print(genotype)
  cellchat <- cellchat.list[[genotype]]
  
  #Set Idents and see group sizes
  # lvls <- unique(cellchat@meta$CellChatGroup) #this messes up order which we don't want
  # cellchat@meta$CellChatGroup <- as.character(cellchat@meta$CellChatGroup)
  # cellchat@meta$CellChatGroup <- factor(cellchat@meta$CellChatGroup, levels = lvls)
  cellchat <- setIdent(cellchat, ident.use = "CellChatGroup") # set "CellChatGroup" as default cell identity
  groupSize <- as.numeric(table(cellchat@idents))
  
  #Set up interaction database
  CellChatDB <- CellChatDB.human  # use CellChatDB.mouse if running on mouse data
  cellchat@DB <- CellChatDB
  
  #Preprocessing
  cellchat <- subsetData(cellchat) 
  future::plan("multiprocess", workers = 32)
  cellchat <- identifyOverExpressedGenes(cellchat) #takes time
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.human) # project gene expression data onto PPI network (optional)
  
  #Inference of cell-cell communication network
  future::plan('sequential')
  cellchat <- computeCommunProb(cellchat) #takes a very long time
  cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter out if there are only few number of cells in certain cell groups
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  #Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  #Compute centrality a metric used in some of the plots (plots summarising all signalling over each celltype).
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  #Visualise
  dir.create(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/setup_overview_all_cells_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
  dir.create(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/percelltype/'), recursive = TRUE, showWarnings = FALSE)
  mat <- cellchat@net$weight
  for (i in 1:nrow(mat)) {
    celltype_plotted <- colnames(mat)[i]
    print(celltype_plotted)
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    png(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/percelltype/setup_overview_', celltype_plotted, '_', opt$celltype_level, '_sub', opt$cell_number, '.png'))
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     edge.weight.max = max(mat), title.name = rownames(mat)[i], 
                     color.use = level2_colour_dictionary)
    dev.off()
  }
  
  dir.create(paste0(opt$out_pwd, '/', genotype, '/overview_plots/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', genotype, '/overview_plots/sources_targets_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  print(netAnalysis_signalingRole_scatter(cellchat))
  dev.off()
  
  # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                           width = 10, height = 16,
                                           color.use = level2_colour_dictionary)
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", 
                                           width = 10, height = 16,
                                           color.use = level2_colour_dictionary)
  pdf(paste0(opt$out_pwd, '/', genotype, '/overview_plots/main_signal_pathways_outgoing_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'), width = 10, height = 16)
  print(ht1)
  dev.off()
  pdf(paste0(opt$out_pwd, '/', genotype, '/overview_plots/main_signal_pathways_incoming_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'), width = 10, height = 16)
  print(ht2)
  dev.off()
  
  #See major active pathways
  pathways <- cellchat@netP$pathways
  print(pathways)
  
  dir.create(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/perpathway/'), recursive = TRUE, showWarnings = FALSE)
  for (pathway.select in pathways){
    png(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/perpathway/pathway_overview_', pathway.select, '_', opt$celltype_level, '_sub', opt$cell_number, '.png'))
    netVisual_aggregate(cellchat, signaling = pathway.select, layout = "circle", color.use = level2_colour_dictionary)
    dev.off()
    pdf(paste0(opt$out_pwd, '/', genotype, '/interaction_plot/perpathway/pathway_overview_', pathway.select, '_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
    netVisual_aggregate(cellchat, signaling = pathway.select, layout = "circle", color.use = level2_colour_dictionary)
    dev.off()
  } 
  
  dir.create(paste0(opt$out_pwd, '/', genotype, '/signalling_role/perpathway/'), recursive = TRUE, showWarnings = FALSE)
  for (pathway.select in pathways){
    png(paste0(opt$out_pwd, '/', genotype, '/signalling_role/perpathway/signalling_role_', pathway.select, '_', opt$celltype_level, '_sub', opt$cell_number, '.png'))
    netAnalysis_signalingRole_network(cellchat, signaling = pathway.select, width = 8, height = 2.5, font.size = 10)
    dev.off()
  }
  
  #save
  cellchat.list[[genotype]] <- cellchat
  
  dir.create(paste0(opt$out_pwd, '/rds/'), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cellchat, file = paste0(opt$out_pwd, '/rds/cellchat_', genotype, '_setup.rds'))
}

#Make grouped objects for comparisons

#WT mammoplasty vs BRCA1 mastectomy
wt_brca1.list <- list('wt' = cellchat.list[['wt_sub1']], 'brca1' = cellchat.list[['brca1']])
cellchat.wt_brca1 <- mergeCellChat(wt_brca1.list, add.names = names(wt_brca1.list), merge.data = TRUE)

#WT mammoplasty vs BRCA2 mastectomy
wt_brca2.list <- list('wt' = cellchat.list[['wt_sub2']], 'brca2' = cellchat.list[['brca2']])
cellchat.wt_brca2 <- mergeCellChat(wt_brca2.list, add.names = names(wt_brca2.list), merge.data = TRUE)


#Save merged objects
saveRDS(cellchat.wt_brca1, file = paste0(opt$out_pwd, '/rds/cellchat_merged_wtbrca1.rds'))
saveRDS(cellchat.wt_brca2, file = paste0(opt$out_pwd, '/rds/cellchat_merged_wtbrca2.rds'))


































