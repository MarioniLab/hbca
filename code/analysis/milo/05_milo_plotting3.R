#! make EXTRA plots (2) for milo analysis. Join ST (sample_types)
#scVI new

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
suppressMessages(library(purrr))

#example input
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/' #epi/output/prop0.3k50'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo outputs.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Functions
celltype_list <- c("LP1", "LP2", "LP3", "LP4", #"LP5": '#C967C9', #LP5=Doublet now labelled so
                   "LP5", "HS1", "HS2", "HS3",
                   "HS4", "BSL1", "BSL2",
                   "DDC1", "DDC2",
                   "FB1", "FB2", "FB3", "FB4", "FB5",
                   "VM1", "VM2", "VM3", "VM4", "VM5",
                   "EC venous", "EC capillary", "EC arterial", 
                   "EC angiogenic tip", "LEC1", "LEC2",
                   'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T',
                   'IFNG+ T', 'NK1', 'NK2', 'NK3',
                   'ILC3', 'B cell', 'Plasma cell', 'Macrophage')
cell_order_rev <- rev(c("LP1", "LP2", "LP3", "LP4", 
                        #"LP proliferating", 
                        "HS1", "HS2", "HS3",
                        "HS4", "BSL1", "BSL2",
                        #"Other epithelial (1)", "Other epithelial (2)",
                        "FB1", "FB2", "FB3", "FB4", "FB5",
                        "VM1", "VM2", "VM3", "VM4", "VM5",
                        "EC venous", "EC capillary", "EC arterial", 
                        "EC angiogenic tip", "LEC1", "LEC2",
                        'CD8T 1', 'CD8T 2', 'CD8T 3', 'CD4T', #CD4T naive
                        'IFNG+ T', 'NK1', 'NK2', 'NK3',
                        'ILC3', 'B cell', 'Plasma cell', 'Macrophage'))

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
      #geom_path(data = da.agg, aes(size=1, colour = '#000000')) +  ##New line y=level2, x=logFC
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
  
  # da.res %>%
  #   mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  #   mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
  #   arrange(group_by) %>%
  #   mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
  #   mutate(pos_x = pos_x, pos_y=pos_y) %>%
  # ggplot(aes(pos_x, pos_y, color=logFC_color)) +
  # scale_color_gradient2() +
  # guides(color="none") +
  # xlab(group.by) + ylab("Log Fold Change") +
  # scale_x_continuous(
  #   breaks = seq(1,n_groups),
  #   labels = setNames(levels(da.res$group_by), seq(1,n_groups))
  # ) +
  # scale_y_continuous(c(-5,5)) + ##New line
  # geom_point() +
  # geom_path(data = da.agg, aes(size=1, colour = '#000000')) +  ##New line y=level2, x=logFC
  # coord_flip() +
  # theme_bw(base_size=22) +
  # theme(strip.text.y =  element_text(angle=0))
  
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

agg_joined_list <- list()
for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM')){
  block_var <- 'parity_age'
  block_var_category <- 'parity'
  if (test_var %in% c('WT_BRCA1PM')){
    tc_subset <- c('Mammoplasty WT', 'Mastectomy BRCA1')
  } else {
    tc_subset <- c('Mammoplasty WT', 'Mastectomy BRCA2')
  }
  
  #Create agg dataset for each celltype
  first = TRUE
  for (celltype_load in c('imm_2', 'epi_2', 'str_2')) {
    print(test_var)
    print(block_var)
    print(celltype_load)
    
    #Load data
    milo_subset <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/epi_nghd_', block_var_category, '.rds'))
    da_joined <- read.csv(paste0(prefix, 'join_2/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
    da_joined$milo_subset <- unlist(map(strsplit(as.character(da_joined$X.1), "[.]"), 1)) 
    
    metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
    rownames(metadata_all) <- metadata_all$cellID
    metadata_sub <- metadata_all[, c('level1','level2', 'patientID', 'tissue_condition')]
    metadata_sub$mean_logFC <- NaN
    
    #get the nhoods that each cell is in so we can work out its average logFC (across neighbourhoods)
    nhoods <- milo_subset@nhoods
    for (cellID_select in rownames(nhoods)) {
      milo_subset_nhoods <- (1:dim(nhoods)[2])[nhoods[cellID_select,] > 0] #colnames(nhoods)[nhoods[cellID_select,] > 0]
      level2_select <- milo_subset$level2[milo_subset$cellID == cellID_select]
      cell_logFCs <- da_joined$logFC[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load) & (da_joined$level2 == level2_select)] #make sure we are getting the right neighborhoods to average over.
      #cell_logFCs <- da_joined$logFC[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load)] #seems like most cells don't live in a 
      #cell_level2s <- da_joined$level2[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load) & (da_joined$level2 == level2_select)] 
      
      cell_mean_logFC <- mean(cell_logFCs)
      metadata_sub[cellID_select, 'mean_logFC'] <- cell_mean_logFC
      
      #print(cell_mean_logFC)
      #print(metadata_sub[cellID_select, 'mean_logFC'])
    }
    
    #Aggregate to find mean and sd per level2 subcluster:
    metadata_sub2 <- metadata_sub[!is.na(metadata_sub$mean_logFC) & (metadata_sub$tissue_condition %in% tc_subset), ]
    
    metadata_agg <- metadata_sub2 %>% 
      group_by(level2, patientID) %>% 
      summarise(mean_lvl2_patientID = mean(mean_logFC)) #, sd = sd(logFC)
    
    metadata_agg2 <- metadata_agg %>% 
      group_by(level2) %>% 
      summarise(mean_lvl2 = mean(mean_lvl2_patientID), sd = sd(mean_lvl2_patientID)) 
    
    metadata_agg2 <- metadata_agg2[!is.na(metadata_agg2$level2),]
    
    if (first == TRUE) {
      joined_agg_celltype <- metadata_agg2 
      first <- FALSE
    } else {
      joined_agg_celltype <- rbind(joined_agg_celltype, metadata_agg2)
    }
  }
  
  joined_agg_celltype$level2 <- factor(joined_agg_celltype$level2,
                                       levels = cell_order_rev)
  joined_agg_celltype <- arrange(joined_agg_celltype, joined_agg_celltype$level2)
  
  #Save
  write.csv(joined_agg_celltype, file = paste0(prefix, 'join_2/output/data/', test_var, '/agg_joined-block_', block_var, '-all.csv'))
  
  agg_joined_list[[test_var]] <- joined_agg_celltype
}


#make plots
da_agg_br1 <- agg_joined_list[['WT_BRCA1PM']]
da_agg_br1$level2 <- factor(da_agg_br1$level2,
                            levels = cell_order_rev)
da_agg_br1$type <- 'HR-BR1'
da_agg_br2 <- agg_joined_list[['WT_BRCA2PM']]
da_agg_br2$level2 <- factor(da_agg_br2$level2,
                            levels = cell_order_rev)
da_agg_br2$type <- 'HR-BR2'
da_agg_joint <- rbind(da_agg_br1, da_agg_br2)
da_agg_joint <- da_agg_joint[!is.na(da_agg_joint$level2),]


dir.create(paste0(prefix, 'join_2/output/aggregated/all/'), showWarnings = F, recursive = T)

pdf(paste0(prefix, 'join_2/output/aggregated/all/milo_sig_BR1BR2_retry.pdf'), width=12, height = 8)
print(ggplot(da_agg_joint[da_agg_joint$type != 'HR-unk',], mapping = aes(y = level2, x = mean_lvl2, 
                                                                         color = type, group = type)) +
        geom_point(size=3) +
        geom_errorbarh(aes(xmin=mean_lvl2-sd, xmax=mean_lvl2+sd), size=0.3, alpha = 0.5) +
        geom_path(size=1) +
        xlim(-5,5) +
        #scale_x_continuous(-5, 5) + 
        scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8')) + 
        theme_minimal()
        )
dev.off()

pdf(paste0(prefix, 'join_2/output/aggregated/all/milo_sig_BR1BR2_zoomedin_retry.pdf'), width=12, height = 8)
print(ggplot(da_agg_joint[da_agg_joint$type != 'HR-unk',], mapping = aes(y = level2, x = mean_lvl2, 
                                                                         color = type, group = type)) +
        geom_point(size=3) +
        geom_errorbarh(aes(xmin=mean_lvl2-sd, xmax=mean_lvl2+sd), size=0.3, alpha = 0.75) +
        geom_path(size=1) +
        xlim(-3,3) +
        #scale_x_continuous(-5, 5) + 
        scale_colour_manual(values = c('HR-BR1' = '#4E79A7', 'HR-BR2' = '#A0CBE8')) + 
        theme_minimal()
)
dev.off()






#rereading data
# for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM')){
#   block_var <- 'parity_age'
#   block_var_category <- 'parity'
# 
#   joined_agg_celltype <- read.csv(file = paste0(prefix, 'join_2/output/data/', test_var, '/agg_joined-block_', block_var, '-all.csv'))
# 
#   agg_joined_list[[test_var]] <- joined_agg_celltype
#   }













# #testing
# celltype_load <- 'imm_2'
# block_var <- 'parity_age'
# block_var_category <- 'parity'
# milo_subset <- readRDS(paste0(prefix, celltype_load, '/output/prop0.3k50/milo_objects/epi_nghd_', block_var_category, '.rds'))
# 
# 
# test_var <- 'WT_BRCA1PM'
# 
# 
# #Load data
# da_joined <- read.csv(paste0(prefix, 'join_2/output/data/', test_var, '/da_joined-block_', block_var, '-all.csv'))
# da_joined$milo_subset <- unlist(map(strsplit(as.character(da_joined$X.1), "[.]"), 1)) 
# 
# metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
# rownames(metadata_all) <- metadata_all$cellID
# metadata_sub <- metadata_all[, c('level1','level2', 'patientID', 'tissue_condition')]
# metadata_sub$mean_logFC <- NaN
# 
# #get the nhoods that each cell is in so we can work out its average logFC
# 
# 
# #AM TESTETIG ON CELLS THAT ARE LATER REMOVED - fix this for a speed up
# nhoods <- milo_subset@nhoods
# for (cellID_select in rownames(nhoods)) {
#   milo_subset_nhoods <- (1:dim(nhoods)[2])[nhoods[cellID_select,] > 0] #colnames(nhoods)[nhoods[cellID_select,] > 0]
#   level2_select <- milo_subset$level2[milo_subset$cellID == cellID_select]
#   cell_logFCs <- da_joined$logFC[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load) & (da_joined$level2 == level2_select)] #make sure we are getting the right neighborhoods to average over.
#   #cell_logFCs <- da_joined$logFC[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load)] #seems like most cells don't live in a 
#   #cell_level2s <- da_joined$level2[(da_joined$Nhood %in% milo_subset_nhoods) & (da_joined$milo_subset == celltype_load) & (da_joined$level2 == level2_select)] 
#   
#   cell_mean_logFC <- mean(cell_logFCs)
#   metadata_sub[cellID_select, 'mean_logFC'] <- cell_mean_logFC
#   
#   print(cell_mean_logFC)
#   #print(metadata_sub[cellID_select, 'mean_logFC'])
# }
# 
# 
# #Make plots of 'milo_signature' averaging over level2 group per patient
# metadata_sub2 <- metadata_sub[!is.na(metadata_sub$mean_logFC) & (metadata_sub$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1')), ]
# 
# metadata_agg <- metadata_sub2 %>% 
#   group_by(level2, patientID) %>% 
#   summarise(mean_lvl2_patientID = mean(mean_logFC)) #, sd = sd(logFC)
# 
# metadata_agg2 <- metadata_agg %>% 
#   group_by(level2) %>% 
#   summarise(mean_lvl2 = mean(mean_lvl2_patientID), sd = sd(mean_lvl2_patientID)) 
# 
# metadata_agg2 <- metadata_agg2[!is.na(metadata_agg2$level2),]
# metadata_agg2$level2 <- factor(metadata_agg2$level2,
#                              levels = cell_order_rev)
# metadata_agg2 <- arrange(metadata_agg2, metadata_agg2$level2)
# 
# 

















