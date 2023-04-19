#!/usr/bin/env Rscript
# Complete overview of signalling pathways for the detailed cell types.

# #Example inputs

# #Level2
# opt <- list()
# opt$sce_pw <- '/nfs/research/marioni/areed/projects/hbca/cellchat/2022-04-05/scvi_new/brca_epistrimm/output/input/HBCA_postlabelling.rds'
# opt$celltype_level <- 'level2'
# opt$out_pwd <- '/nfs/research/marioni/areed/projects/hbca/cellchat/2022-04-05/scvi_new/brca_epistrimm/output'
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
option_list = list(make_option(c('--celltype_level'),
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
cellchat.wt_sub1 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_wt_sub1_setup.rds'))
cellchat.wt_sub2 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_wt_sub2_setup.rds'))
cellchat.brca1 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_brca1_setup.rds'))
cellchat.brca2 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_brca2_setup.rds'))

#need this for later plots
# cellchat.wt_sub1 <- netAnalysis_computeCentrality(cellchat.wt_sub1, slot.name = "net")
# cellchat.wt_sub2 <- netAnalysis_computeCentrality(cellchat.wt_sub1, slot.name = "net")
# cellchat.brca1 <- netAnalysis_computeCentrality(cellchat.brca1, slot.name = "net")
# cellchat.brca2 <- netAnalysis_computeCentrality(cellchat.brca2, slot.name = "net")
# cellchat.wt_sub1 <- netAnalysis_computeCentrality(cellchat.wt_sub1, slot.name = "netP")
# cellchat.wt_sub2 <- netAnalysis_computeCentrality(cellchat.wt_sub1, slot.name = "netP")
# cellchat.brca1 <- netAnalysis_computeCentrality(cellchat.brca1, slot.name = "netP")
# cellchat.brca2 <- netAnalysis_computeCentrality(cellchat.brca2, slot.name = "netP")

wt_brca1.list <- list('wt' = cellchat.wt_sub1, 'brca1' = cellchat.brca1)
wt_brca2.list <- list('wt' = cellchat.wt_sub2, 'brca2' = cellchat.brca2)

cellchat.wt_brca1 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_merged_wtbrca1.rds'))
cellchat.wt_brca2 <- readRDS(paste0(opt$out_pwd, '/rds/cellchat_merged_wtbrca2.rds'))
cellchat.list <- list('wt_brca1' = cellchat.wt_brca1, 'wt_brca2' = cellchat.wt_brca2)

#colours
level2_colour_dictionary = c('LP1' = "#DA80DA", 'LP2' = "#815481", 'LP3' = "#C040C0", 'LP4' = "#E1AFE1", 
                             'HS1' = "#EDABB9", 'HS2' = "#EB5C79", 'HS3' = "#A06A75", 'HS4' = "#C00028",
                             'BSL1' = "#EB675E", 'BSL2' = "#A23E36",
                             'FB1' = "#DFA38A", 'FB2' = "#8C3612", 'FB3' = "#623623", 'FB4' = "#916350", 'FB5' = "#DAC3C3",
                             'VM1' = "#F8770B", 'VM2' = "#E09E3A", 'VM3' = "#CD7225", 'VM4' = "#FFC990", 'VM5' = "#AC5812",
                             'EC venous' = "#FEE083", 'EC capillary' = "#897538", 'EC arterial' = "#E7B419", 'EC angiogenic tip' = "#BCA048",
                             'LEC1' = "#6F8BE2", 'LEC2' = "#3053BC",
                             "Lymphoid" = '#9FC5E8', 
                             "Myeloid" = '#AAB256') 

for (test_type in c('wt_brca1','wt_brca2')){
  print(test_type)
  cellchat <- cellchat.list[[test_type]]
  
  #Compare number of interactions and their strength
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  
  dir.create(paste0(opt$out_pwd, '/', test_type, '/overview_plots'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', test_type, '/overview_plots/compare_interactions_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  print(gg1 + gg2)
  dev.off()
  
  #Differential number of interactions or interaction strength among different cell populations. 
  #Red means up in WT, blue means up in the test data set BRCA1/2.
  dir.create(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/circleplot_overall_diffinteraction_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  dev.off()
  
  #Same summary as a heatmap
  gg1 <- netVisual_heatmap(cellchat)
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  pdf(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/heatmap_overall_diffinteraction_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  print(gg1 + gg2)
  dev.off()
  
  #Compare the overall information flow of each signaling pathway
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  dir.create(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_overall_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
  print(gg1 + gg2)
  dev.off()
  pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_overall_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_rel.pdf')) #relaltive boxplot only
  print(gg1)
  dev.off()
  
  #Use specified tagets only
  if (test_type == 'wt_brca1') {
    targets = c(1:4)
    target_name = 'LP'
  } else {
    targets = 24
    target_name = 'ECangio'
  }
  
  
  gg1 <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = T, do.flip = T, do.stat = TRUE)
  gg1_flip <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = T, do.flip = F, do.stat = TRUE)
  gg2 <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = F, do.flip = T, do.stat = TRUE)
  dir.create(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_target', target_name, '_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
  print(gg1 + gg2)
  dev.off()
  pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_target', target_name, '_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_rel.pdf')) #relaltive boxplot only (flipped)
  print(gg1_flip)
  dev.off()
  
  
  #Identify dysfunctional signaling by comparing the communication probabities
  gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Decreased signaling in WT", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in WT", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/dotplot_diffcommunprob_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
  print(gg1 + gg2)
  dev.off()
}

print('Done loop 1.')

cellchat.wt_brca.list <- list('wt_brca1' = wt_brca1.list, 'wt_brca2' = wt_brca2.list)
for (list_name in c('wt_brca1','wt_brca2')) {
  print(list_name)
  test.list <- cellchat.wt_brca.list[[list_name]]
  
  #plot the changes in signalling on a celltype specific basis
  # lvls <- unique(test.list[[1]]@meta$CellChatGroup)
  # dir.create(paste0(opt$out_pwd, '/', list_name, '/scatter_signalling_changes/percelltype/'), recursive = TRUE, showWarnings = FALSE)
  # for (celltype_plot in lvls) {
  #   print(celltype_plot)
  #   # png(paste0(opt$out_pwd, '/', list_name, '/scatter_signalling_changes/percelltype/scatter_net_', celltype_plot, '_', opt$celltype_level, '_sub', opt$cell_number, '.png'))
  #   # print(netAnalysis_signalingChanges_scatter(test.list, idents.use = celltype_plot, slot.name = 'net')) #number interactions
  #   # dev.off()
  #   png(paste0(opt$out_pwd, '/', list_name, '/scatter_signalling_changes/percelltype/scatter_', celltype_plot, '_', opt$celltype_level, '_sub', opt$cell_number, '.png'))
  #   print(netAnalysis_signalingChanges_scatter(test.list, 
  #                                              idents.use = celltype_plot, 
  #                                              slot.name = 'netP',
  #                                              #color.use = level2_colour_dictionary
  #                                              )) #strength of interactions
  #   dev.off()
  # }
  
  #plot of overall incoming outgoing signalling scatter per celltype
  num.link <- sapply(test.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(test.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(test.list[[i]], 
                                                 title = names(test.list)[i],
                                                 color.use = level2_colour_dictionary,
                                                 weight.MinMax = weight.MinMax)
  }
  
  dir.create(paste0(opt$out_pwd, '/', list_name, '/overview_plots/'), recursive = TRUE, showWarnings = FALSE)
  pdf(paste0(opt$out_pwd, '/', list_name, '/overview_plots/scatter_signalling_role_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
  print(patchwork::wrap_plots(plots = gg))
  dev.off()
}


#Pathways most incoming to 'EC angiogenesis'
ec_angio_pathways <- c('SEMA3', 'SEMA6', 'VEGF', 'PECAM1', 'CALCR', 'CD46', 'ESAM', 'ANGPT', 'VISFATIN', 'ncWNT', 'ANGPTL') #removed TNF, CD96, CDH5 for ncWNT and ANGPTL


ht1 <- netAnalysis_signalingRole_heatmap(cellchat,
                                         signaling = ec_angio_pathways,
                                         color.use = level2_colour_dictionary,
                                         pattern = "outgoing", 
                                         width = 10,
                                         height = 16)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,
                                         signaling = ec_angio_pathways,
                                         color.use = level2_colour_dictionary,
                                         pattern = "incoming", 
                                         width = 10,
                                         height = 16)
pdf(paste0(opt$out_pwd, '/brca2/overview_plots/ec_angiogenisis_pathways_outgoing_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'), width = 10, height = 16)
print(ht1)
dev.off()
pdf(paste0(opt$out_pwd, '/brca2/overview_plots/ec_angiogenesis_pathways_incoming_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'), width = 10, height = 16)
print(ht2)
dev.off()



##compare all three groups

wt_brca1_brca2.list <- list('wt' = cellchat.wt_sub1, 'brca1' = cellchat.brca1, 'brca2' = cellchat.brca2)

num.link <- sapply(wt_brca1_brca2.list , function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(wt_brca1_brca2.list )) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(wt_brca1_brca2.list [[i]], 
                                               title = names(wt_brca1_brca2.list )[i],
                                               color.use = level2_colour_dictionary,
                                               weight.MinMax = weight.MinMax)
}

dir.create(paste0(opt$out_pwd, '/wt_brca1_brca2/overview_plots/'), recursive = TRUE, showWarnings = FALSE)
pdf(paste0(opt$out_pwd, '/wt_brca1_brca2/overview_plots/scatter_signalling_role_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
print(patchwork::wrap_plots(plots = gg))
dev.off()



# ###Copying
# 
# for (test_type in c('wt_brca1','wt_brca2')){
#   print(test_type)
#   cellchat <- cellchat.list[[test_type]]
#   gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
#   gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
#   dir.create(paste0(opt$out_pwd, '/', test_type, '/overview_plots'), recursive = TRUE, showWarnings = FALSE)
#   pdf(paste0(opt$out_pwd, '/', test_type, '/overview_plots/compare_interactions_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
#   print(gg1 + gg2)
#   dev.off()
#   dir.create(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/'), recursive = TRUE, showWarnings = FALSE)
#   pdf(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/circleplot_overall_diffinteraction_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
#   netVisual_diffInteraction(cellchat, weight.scale = T)
#   netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
#   dev.off()
#   gg1 <- netVisual_heatmap(cellchat)
#   gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#   pdf(paste0(opt$out_pwd, '/', test_type, '/interaction_plot/heatmap_overall_diffinteraction_', opt$celltype_level, '_sub', opt$cell_number, '.pdf'))
#   print(gg1 + gg2)
#   dev.off()
#   gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
#   gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
#   dir.create(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/'), recursive = TRUE, showWarnings = FALSE)
#   pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_overall_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
#   print(gg1 + gg2)
#   dev.off()
#   pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_overall_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_rel.pdf')) #relaltive boxplot only
#   print(gg1)
#   dev.off
#   if (test_type == 'wt_brca1') {
#     targets = c(1:4)
#     target_name = 'LP'
#   } else {
#     targets = 24
#     target_name = 'ECangio'
#   }
#   gg1 <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = T, do.flip = T, do.stat = TRUE)
#   gg1_flip <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = T, do.flip = F, do.stat = TRUE)
#   gg2 <- rankNet(cellchat, targets.use=targets, mode = "comparison", stacked = F, do.flip = T, do.stat = TRUE)
#   dir.create(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/'), recursive = TRUE, showWarnings = FALSE)
#   pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_target', target_name, '_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
#   print(gg1 + gg2)
#   dev.off()
#   pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/barplot_target', target_name, '_diffinfoflow_', opt$celltype_level, '_sub', opt$cell_number, '_rel.pdf'), width=10, height=4) #relaltive boxplot only (flipped)
#   print(gg1_flip)
#   dev.off()
#   gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Decreased signaling in WT", angle.x = 45, remove.isolate = T)
#   gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in WT", angle.x = 45, remove.isolate = T)
#   pdf(paste0(opt$out_pwd, '/', test_type, '/differential_signalling/dotplot_diffcommunprob_', opt$celltype_level, '_sub', opt$cell_number, '_both.pdf')) #both relative and true size boxplots
#   print(gg1 + gg2)
#   dev.off()
# }





