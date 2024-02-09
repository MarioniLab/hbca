#R

#Make more plots (x2) for milo analysis

#Use the same neighborhoods for all AR vs HR tests
#scVI embedding

#libraries
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(optparse))
suppressMessages(library(MASS))
suppressMessages(library(ggrastr))
suppressMessages(library(ggbeeswarm))

# #example options (testing)

# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo outputs.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Functions
celltype_list <- c("LASP1", "LASP2", "LASP3", "LASP4", "LASP5",
                   "LHS1", "LHS2", "LHS3",
                   "BMYO1", "BMYO2", 
                   "DDC1", "DDC2",
                   "FB1", "FB2", "FB3", "FB4",
                   "PV1", "PV2", "PV3", "PV4", "PV5",
                   "VEV", "VEC", "VEA","VEAT", 
                   "LE1", "LE2",
                   'CD4_naive', "CD4_Th", 
                   'CD8_Tem', 'CD8_Trm',
                   'CD8_Tc1',
                   'NKT', 'NK',
                   'ILC', 
                   'B_naive', 'B_mem_switched', 'B_mem_unswitched', 
                   'Plasma_cell', 
                   'Macro', 'Macro-lipo', 'DC')
celltype_list_reduced <- celltype_list[!(celltype_list %in% c('DDC1', 'DDC2', 'LASP5'))]


plotDAbeeswarm2 <- function(da.res, da.agg, celltype_list, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  #make consistent levels
  celltype_list_short <- celltype_list[(celltype_list %in% unique(da.res$level2)) & (celltype_list %in% unique(da.agg$level2))]
  da.res$level2 <- factor(da.res$level2, levels = rev(celltype_list_short))
  da.agg$level2 <- factor(da.agg$level2, levels = rev(celltype_list_short))
  da.agg$group_by <- factor(da.agg$level2, levels = rev(celltype_list_short))
  
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
      # stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }
  
  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, levels=unique(group_by)))
    # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }
  
  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }
  
  # Get position with ggbeeswarm
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )
  beeswarm_pos2 <- ggplot_build(
      ggplot(data = da.agg, aes(level2, logFC)) +
      geom_point(group = NA, size=1)   ##New line y=level2, x=logFC
  )
  
  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y
  pos_x2 <- beeswarm_pos2$data[[1]]$y #da.agg$level2
  pos_y2 <- rep(0, length(beeswarm_pos2$data[[1]]$x))  #da.agg
  
  n_groups <- unique(da.res$group_by) %>% length()

  da.res2 <- da.res %>%      
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y=pos_y)
  da.agg2 <- da.agg %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    mutate(pos_x = pos_x2, pos_y=pos_y2)
  
  ggplot(data = da.res2, aes(pos_x, pos_y, color=logFC_color)) +
    scale_color_gradient2() +
    guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1,n_groups),
      labels = setNames(levels(da.res$group_by), seq(1,n_groups))
    ) +
    scale_y_continuous(c(-5,5)) + ##New line
    geom_point() +
    #geom_path(data = da.agg2, aes(size=1, colour = '#000000')) +  ##New line y=level2, x=logFC
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))
}


# Analysis

prefix <- opt$input_path

da_joined_list <- list()
da_agg_list <- list()

cell_order_all_rev <- rev(celltype_list)
cell_order_rev <- rev(celltype_list_reduced)

for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM', 'patient_age', 'parity')){
  if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')){
    block_var <- 'parity_age'
  } else if (test_var == 'patient_age'){
    block_var <- 'parity'
  } else if (test_var == 'parity') {
    block_var <- 'patient_age'
  } else {
    block_var <- 'none'
  }
  
  if (block_var %in% c('none', 'patient_age')){
    block_var_category <- 'none'
  } else {
    block_var_category <- 'parity'
  }
  
  da_joined_list[[test_var]] <- da_joined <- read.csv(paste0(prefix, 'join/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
  da_joined <- da_joined[da_joined$level2 %in% cell_order_rev,]
  

  da_join_agg <- da_joined_list[[test_var]] %>% 
    group_by(level2) %>% 
    summarise(mean = mean(logFC), sd = sd(logFC))
  
  da_join_agg$level2 <- factor(da_join_agg$level2,
                               levels = cell_order_rev)
  da_join_agg <- arrange(da_join_agg, da_join_agg$level2)
  
  da_agg_list[[test_var]] <- da_join_agg
}

#Ageing Parity
da_agg_age <- da_agg_list[['patient_age']]
# da_agg_age$mean_og <- da_agg_age$mean
da_agg_age$mean <- da_agg_age$mean * 20 #Not sure this makes any sense really
da_agg_age$level2 <- factor(da_agg_age$level2,
                            levels = cell_order_rev)
da_agg_age$type <- 'age'
da_agg_parity <- da_agg_list[['parity']]
da_agg_parity$level2 <- factor(da_agg_parity$level2,
                            levels = cell_order_rev)
da_agg_parity$type <- 'parity'
da_agg_joint_AP <- rbind(da_agg_age, da_agg_parity)
da_agg_joint_AP <- da_agg_joint_AP[!is.na(da_agg_joint_AP$level2),]

#genotype
da_agg_br1 <- da_agg_list[['WT_BRCA1PM']]
da_agg_br1$level2 <- factor(da_agg_br1$level2,
                            levels = cell_order_rev)
da_agg_br1$type <- 'HR-BR1'
da_agg_br2 <- da_agg_list[['WT_BRCA2PM']]
da_agg_br2$level2 <- factor(da_agg_br2$level2,
                            levels = cell_order_rev)
da_agg_br2$type <- 'HR-BR2'
da_agg_WTPM <- da_agg_list[['WT_WTPM']]
da_agg_WTPM$level2 <- factor(da_agg_WTPM$level2,
                            levels = cell_order_rev)
da_agg_WTPM$type <- 'HR-unk'
da_agg_joint <- rbind(da_agg_br1, da_agg_br2, da_agg_WTPM)
da_agg_joint <- da_agg_joint[!is.na(da_agg_joint$level2),]

#Make a 'milo signature plot for BRCA1 and BRCA2

dir.create(paste0(prefix, 'join/output/aggregated/all/'), showWarnings = F, recursive = T)

#This fails
# pdf(paste0(prefix, 'join/output/aggregated/all/beeswarm_level2-block_', block_var, '-all4.pdf'), width=12, height = 8)
# print(ggplot() +
#         geom_point(da_agg_br1,
#                    mapping = aes(x = level2, y = mean, size=3,
#                                  color = '#4e79a7')) +
#         geom_path(da_agg_br1,
#                   mapping = aes(x = level2, y = mean, size=1, 
#                                 color = '#4e79a7', group = NA)
#                   ) +
#         geom_point(da_agg_br2,
#                    mapping = aes(x = level2, y = mean, size=3,
#                                  color = '#a0cbe8')) +
#         geom_path(da_agg_br2,
#                   mapping = aes(x = level2, y = mean, size=1, 
#                                 color = '#a0cbe8', group = NA)
#                   ) +
#         geom_point(da_agg_WTPM,
#                    mapping = aes(x = level2, y = mean, size=3,
#                                  color = '#bab0ac')) +
#         geom_path(da_agg_WTPM,
#                   mapping = aes(x = level2, y = mean, size=1, 
#                                 color = '#bab0ac',group = NA)
#                   ) +
#         coord_flip() +
#         ylim(-5, 5))
# dev.off()

#Age Parity
pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_AgeParity.pdf'), width=12, height = 8)
print(ggplot(da_agg_joint_AP, mapping = aes(y = level2, x = mean, 
                                         color = type, group = type)) +
        geom_point(size=3) +
        geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), width=.1) +
        geom_path(size=1) +
        scale_x_continuous(-5, 5) + 
        scale_colour_manual(values = c('age' = 'yellow', 'parity' = 'pink')) + 
        theme_minimal())
dev.off()


#genotype
# pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_all.pdf'), width=12, height = 8)
# print(ggplot(rbind(da_agg_joint, da_agg_joint_AP), mapping = aes(y = level2, x = mean, 
#                                          color = type, group = type)) +
#         geom_point(size=3) +
#         geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), width=.1) +
#         geom_path(size=1) +
#         scale_x_continuous(-5, 5) + 
#         scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8', 'HR-unk' = '#BAB0AC', 'age' = 'yellow', 'parity' = 'pink')))
# dev.off()
# 
# 
# pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_BR1BR2WT.pdf'), width=12, height = 8)
# print(ggplot(da_agg_joint, mapping = aes(y = level2, x = mean, 
#                                          color = type, group = type)) +
#         geom_point(size=3) +
#         geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), width=.1) +
#         geom_path(size=1) +
#         scale_x_continuous(-5, 5) + 
#         scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8', 'HR-unk' = '#BAB0AC')))
# dev.off()

pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_BR1BR2_retry.pdf'), width=12, height = 8)
print(ggplot(da_agg_joint[da_agg_joint$type != 'HR-unk',], mapping = aes(y = level2, x = mean, 
                                                                         color = type, group = type)) +
        geom_point(size=3) +
        geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), size=0.3, alpha = 0.5) +
        geom_path(size=1) +
        xlim(-5,5) +
        #scale_x_continuous(-5, 5) + 
        scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8')) + 
        theme_minimal()
)
dev.off()

pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_BR1BR2_zoomedin_retry.pdf'), width=12, height = 8)
print(ggplot(da_agg_joint[da_agg_joint$type != 'HR-unk',], mapping = aes(y = level2, x = mean, 
                                                                         color = type, group = type)) +
        geom_point(size=3) +
        geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), size=0.3, alpha = 0.75) +
        geom_path(size=1) +
        xlim(-3,3) +
        #scale_x_continuous(-5, 5) + 
        scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8')) + 
        theme_minimal()
)
dev.off()

pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_BRCA1.pdf'), width=12, height = 8)
print(ggplot() +
        geom_point(da_agg_list[['WT_BRCA1PM']],
                   mapping = aes(y = level2, x = mean, size=3,
                                 color = '#4e79a7')) +
        geom_path(da_agg_list[['WT_BRCA1PM']],
                  mapping = aes(y = level2, x = mean, size=1, 
                                color = '#4e79a7'),
                  group = NA) +
        scale_x_continuous(-5, 5))
dev.off()
pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_BRCA2.pdf'), width=12, height = 8)
print(ggplot() +
        geom_point(da_agg_list[['WT_BRCA2PM']],
                   mapping = aes(y = level2, x = mean, size=3,
                                 color = '#a0cbe8')) +
        geom_path(da_agg_list[['WT_BRCA2PM']],
                  mapping = aes(y = level2, x = mean, size=1, 
                                color = '#a0cbe8'),
                  group = NA) +
        scale_x_continuous(-5, 5))
dev.off()
# pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_WTPM.pdf'), width=12, height = 8)
# print(ggplot() +
#         geom_point(da_agg_list[['WT_WTPM']],
#                    mapping = aes(y = level2, x = mean, size=3,
#                                  color = '#bab0ac')) +
#         geom_path(da_agg_list[['WT_WTPM']],
#                   mapping = aes(y = level2, x = mean, size=1, 
#                                 color = '#bab0ac'),
#                   group = NA) +
#         scale_x_continuous(-5, 5))
# dev.off()







# #Look at WTPM patients singularly
# 
# test_var <- 'WT_WTPM'
# 
# milo_list <- list()
# da_list <- list()
# da_1715PM_list <- list()
# da_2101PM_list <- list()
# 
# for (celltype_load in c('epi', 'str', 'imm')) {
#   milo_list[[celltype_load]] <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/epi_nghd_', test_var, '_', block_var_category, '.rds'))
#   da_list[[celltype_load]] <- read.csv(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_testing/', test_var, '/da-block_', block_var, '-all.csv'))
# 
#   da_list[[celltype_load]] <- annotateNhoods(milo_list[[celltype_load]],
#                                              da_list[[celltype_load]],
#                                              coldata_col = "level1")
#   da_list[[celltype_load]] <- annotateNhoods(milo_list[[celltype_load]],
#                                              da_list[[celltype_load]],
#                                              coldata_col = "level2")
#   
#   #For later exploration also get the individual neighbourhoods per patient (HR-unk)
#   ct_nhoodCounts <- milo_list[[celltype_load]]@nhoodCounts
#   
#   #1715PM
#   patientID_select = '1715PM'
#   sampleIDs <- unique(milo_list[[celltype_load]]$sampleID[milo_list[[celltype_load]]$patientID == patientID_select])
#   imm_nhoodCounts_sub <- imm_nhoodCounts[, colnames(imm_nhoodCounts) %in% sampleIDs]
#   imm_nhoodCounts_incl <- da_list[[celltype_load]]$Nhood[rowSums(data.frame(as.matrix(imm_nhoodCounts_sub))) > 0]
#   
#   da_1715PM_list[[celltype_load]] <- da_list[[celltype_load]][da_list[[celltype_load]]$Nhood %in% imm_nhoodCounts_incl, ]
#   
#   #2101PM
#   patientID_select = '2101PM'
#   sampleIDs <- unique(milo_list[[celltype_load]]$sampleID[milo_list[[celltype_load]]$patientID == patientID_select])
#   imm_nhoodCounts_sub <- imm_nhoodCounts[, colnames(imm_nhoodCounts) %in% sampleIDs]
#   imm_nhoodCounts_incl <- da_list[[celltype_load]]$Nhood[rowSums(data.frame(as.matrix(imm_nhoodCounts_sub))) > 0]
#   
#   da_2101PM_list[[celltype_load]] <- da_list[[celltype_load]][da_list[[celltype_load]]$Nhood %in% imm_nhoodCounts_incl, ]
# }
# 
# 
# da_joined_1715PM <- do.call(rbind, da_1715PM_list)
# da_joined_2101PM <- do.call(rbind, da_2101PM_list)
# 
# # da_joined_1715PM$level1 <- factor(da_joined_1715PM$level1,
# #                            levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
# #                                           "Fibroblast", "Endothelial", "Vascular mural",
# #                                           "Lymphoid", "Myeloid")))
# da_joined_1715PM$level2 <- factor(da_joined_1715PM$level2,
#                            levels = cell_order_all_rev)
# # da_joined_2101PM$level1 <- factor(da_joined_2101PM$level1,
# #                                   levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
# #                                                  "Fibroblast", "Endothelial", "Vascular mural",
# #                                                  "Lymphoid", "Myeloid")))
# da_joined_2101PM$level2 <- factor(da_joined_2101PM$level2,
#                                   levels = cell_order_all_rev)
# 
# #save
# write.csv(da_joined_1715PM, file = paste0(prefix, 'join/output/data/', test_var, '/da_joined-block_', block_var, '-all_1715PM.csv'))
# write.csv(da_joined_2101PM, file = paste0(prefix, 'join/output/data/', test_var, '/da_joined-block_', block_var, '-all_2101PM.csv'))
# 
# 
# #make aggregated df to plot the mean values
# da_join_1715PM_agg <- aggregate(x = da_joined_1715PM,
#                          by = list(da_joined_1715PM$level2), 
#                          FUN = mean)
# # da_join_1715PM_agg$level1 <- factor(da_join_1715PM_agg$level1,
# #                              levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
# #                                             "Fibroblast", "Endothelial", "Vascular mural",
# #                                             "Lymphoid", "Myeloid")))
# da_join_1715PM_agg$level2 <- factor(da_join_1715PM_agg$level2,
#                              levels = cell_order_all_rev)
# 
# da_join_2101PM_agg <- aggregate(x = da_joined_2101PM,
#                                 by = list(da_joined_2101PM$level2), 
#                                 FUN = mean)
# # da_join_2101PM_agg$level1 <- factor(da_join_2101PM_agg$level1,
# #                                     levels = rev(c("Luminal progenitor", "Luminal hormone sensing", "Basal",
# #                                                    "Fibroblast", "Endothelial", "Vascular mural",
# #                                                    "Lymphoid", "Myeloid")))
# da_join_2101PM_agg$level2 <- factor(da_join_2101PM_agg$level2,
#                                     levels = cell_order_all_rev)
# 
# #Remove celltypes not present across all types in sufficient numbers
# da_joined_1715PM <- da_joined_1715PM[da_joined_1715PM$level2 %in% cell_order_rev,]
# da_join_1715PM_agg <- da_join_1715PM_agg[da_join_1715PM_agg$level2 %in% cell_order_rev,]
# da_joined_2101PM <- da_joined_2101PM[da_joined_2101PM$level2 %in% cell_order_rev,]
# da_join_2101PM_agg <- da_join_2101PM_agg[da_join_2101PM_agg$level2 %in% cell_order_rev,]
# 
# 
# pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_1715PM.pdf'), width=12, height = 8)
# print(ggplot() +
#         geom_point(da_join_1715PM_agg,
#                    mapping = aes(y = level2, x = mean, size=3)) +
#         geom_path(da_join_1715PM_agg,
#                   mapping = aes(y = level2, x = mean, size=1),
#                   group = NA) +
#         scale_x_continuous(-5, 5))
# dev.off()
# 
# pdf(paste0(prefix, 'join/output/aggregated/all/milo_sig_2101PM.pdf'), width=12, height = 8)
# print(ggplot() +
#         geom_point(da_join_2101PM_agg,
#                    mapping = aes(y = level2, x = mean, size=3)) +
#         geom_path(da_join_2101PM_agg,
#                   mapping = aes(y = level2, x = mean, size=1),
#                   group = NA) +
#         scale_x_continuous(-5, 5))
# dev.off()










