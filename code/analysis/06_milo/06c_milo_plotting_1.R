#R

#Make basic plots for milo analysis

#Use the same neighborhoods for all AR vs HR tests
#scVI embedding

#libraries
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))
suppressMessages(library(igraph))

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(optparse))
suppressMessages(library(MASS))
suppressMessages(library(ggrastr))

# #example options (testing)

# #Epi
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/epi/output/prop0.3k50'

# #Str
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/str/output/prop0.3k50'

# #Imm
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/imm/output/prop0.3k50'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve milo object and to use for saving da results.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#Functions
plot_milo <- function(milo, da_out, plot_pwd, block_var, sample_tp){
  milo <- edited.buildNhoodGraph(milo, overlap=20) #20 overlap reduced nhood connection overplotting
  
  #make plots
  cells_to_plot <- (milo$milo_test != 'NA') & !is.na(milo$milo_test)
  umap_pl <- plotReducedDim(milo[,cells_to_plot], dimred = "UMAP", colour_by="milo_test", point_size=0.5) +  guides(fill="none")
  umap_pl_text <- plotReducedDim(milo[,cells_to_plot], dimred = "UMAP", colour_by="milo_test", text_by="level2", point_size=0.5) +  guides(fill="none")
  nh_graph_pl <- plotNhoodGraphDA(milo, da_out, res_column = 'logFC', layout="UMAP", alpha=0.1, node_stroke=0) #'umap'
  
  #rastorize
  umap_pl_rast <- rasterize(umap_pl, layers='GeomPoint', dpi=500)
  umap_pl_text_rast <- rasterize(umap_pl_text, layers='GeomPoint', dpi=500)
  nh_graph_pl_rast <- rasterize(nh_graph_pl, layers=c('GeomPoint', 'GeomEdgeSegment'), dpi=500)
  
  #save no text UMAPs
  dir.create(plot_pwd, showWarnings = FALSE, recursive = TRUE)
  pdf(paste0(plot_pwd, 'umap_nghd-block_', block_var, '-', sample_tp, '.pdf'), width=20,height = 8)
  print(umap_pl + nh_graph_pl +  plot_layout(guides="collect"))
  dev.off()
  
  pdf(paste0(plot_pwd, 'umap_nghd-block_', block_var, '-', sample_tp, '_rasterise.pdf'), width=20, height = 8)
  print(umap_pl_rast + nh_graph_pl_rast)# +  plot_layout(guides="collect"))
  dev.off()
  
  png(paste0(plot_pwd, 'umap_nghd-block_', block_var, '-', sample_tp, '.png'))
  print(umap_pl + nh_graph_pl +  plot_layout(guides="collect"))
  dev.off()
  
  #save with text (level2) UMAPs
  pdf(paste0(plot_pwd, 'umap_nghd_level2-block_', block_var, '-', sample_tp, '_text.pdf'), width=20,height = 8)
  print(umap_pl_text + nh_graph_pl +  plot_layout(guides="collect"))
  dev.off()
  
  pdf(paste0(plot_pwd, 'umap_nghd_level2-block_', block_var, '-', sample_tp, '_rasterise_text.pdf'), width=20, height = 8)
  print(umap_pl_text_rast + nh_graph_pl_rast)# +  plot_layout(guides="collect"))
  dev.off()
  
  png(paste0(plot_pwd, 'umap_nghd_level2-block_', block_var, '-', sample_tp, '_text.png'))
  print(umap_pl_text + nh_graph_pl +  plot_layout(guides="collect"))
  dev.off()
  
  #beeswarm plot
  da_out <- annotateNhoods(milo, da_out, coldata_col = "level2")
  if (sum(da_out$SpatialFDR < 0.1) > 0) {
    pdf(paste0(plot_pwd, 'beeswarm_level2-block_', block_var, '-', sample_tp, '.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_out, group.by = "level2"))
    dev.off()
  } else {
    pdf(paste0(plot_pwd, 'beeswarm_level2-block_', block_var, '-', sample_tp, '_alpha=1.pdf'), width=12, height = 8)
    print(plotDAbeeswarm(da_out, group.by = "level2", alpha=1))
    dev.off()
  }
}

#Milo Functions
#Had to make edits as there was an error that occured in the stromal compartment due to large number of cells/nhoods

edited.build_nhood_adjacency <- function(nhoods, overlap=1){
  nh_intersect_mat <- Matrix::crossprod(nhoods)
  
  #This line causes error "In int2i(as.integer(i), n) : NAs introduced by coercion to integer range"
  #Due to (i think) exceeding R max integer value
  #nh_intersect_mat[nh_intersect_mat < overlap] <- 0
  
  #Instead do per column to reduced the sizes of vector at each step (this takes longer so add if overlap>1 to run - mostly use overlap<1 anyway).
  if (overlap > 1) {
    for (column in 1:dim(nh_intersect_mat)[2]) {
      nh_intersect_mat[nh_intersect_mat[, column] < overlap, column] <- 0
    }
  }
  
  rownames(nh_intersect_mat) <- colnames(nhoods)
  colnames(nh_intersect_mat) <- colnames(nhoods)
  return(nh_intersect_mat)
}

edited.buildNhoodGraph <- function(x, overlap=1){
  
  if(!is(x, "Milo")){
    stop("Not a valid Milo object")
  }
  
  # are neighbourhoods calculated?
  if(ncol(nhoods(x)) == 1 & nrow(nhoods(x)) == 1){
    stop("No neighbourhoods found - run makeNhoods first")
  }
  
  ## Build adjacency matrix for nhoods
  nh_intersect_mat <- edited.build_nhood_adjacency(nhoods(x), overlap=overlap)
  
  # add to slot if empty
  nhoodAdjacency(x) <- nh_intersect_mat
  
  ## Make igraph object
  ig <- graph.adjacency(nh_intersect_mat, mode="undirected", weighted=TRUE)
  nhood_sizes <- colSums(nhoods(x))
  ig <- set_vertex_attr(ig, name = 'size', value = nhood_sizes[vertex.attributes(ig)$name])
  ## Add to nhoodGraph slot in milo object
  nhoodGraph(x) <- ig
  return(x)
}

## Analysis


##choose k and prop with optparse

#load Data
prefix <- opt$input_path



for (test_var in c('WT_BRCA1PM', 'WT_BRCA2PM', 'patient_age', 'parity')){
  for (block_var in c('parity', 'patient_age', 'parity_age')){ 
    print('New round of testing:')
    print(test_var)
    print(block_var)
    
    if (block_var %in% c('none')){
      block_var_category <- 'none'
    } else {
      block_var_category <- 'parity'
    }
    
    #skip if test and block are not testable
    if (test_var == block_var) {
      next
    } else if ((test_var %in% c('patient_age', 'parity')) & (block_var == 'parity_age')) {
      next
    } else if ((test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')) & !(block_var %in% c('parity_age'))){
      next
    }
    
    if (test_var %in% c('WT_BRCA1PM', 'WT_BRCA2PM')) {
      milo <- readRDS(paste0(prefix, '/milo_objects/HR_nghd_', block_var_category, '.rds'))
    } else {
      milo <- readRDS(paste0(prefix, '/milo_objects/AR_nghd_', block_var_category, '.rds'))
    } 
    
    
    # Plotting
    
    #test_var meta
    if (test_var == 'WT_BRCA1PM') {
      milo$milo_test <- 'NA'
      milo$milo_test[milo$tissue_condition %in% c('Mammoplasty WT')] <- 'Mammoplasty WT'
      milo$milo_test[milo$tissue_condition %in% c( 'Mastectomy BRCA1')] <- 'Mastectomy BRCA1'
    } else if (test_var == 'WT_BRCA2PM') {
      milo$milo_test <- 'NA'
      milo$milo_test[milo$tissue_condition %in% c('Mammoplasty WT')] <- 'Mammoplasty WT'
      milo$milo_test[milo$tissue_condition %in% c( 'Mastectomy BRCA2')] <- 'Mastectomy BRCA2'
    } else if (test_var == 'patient_age') {
      milo <- milo[,milo$tissue_condition %in% c('Mammoplasty WT')]
      milo$milo_test <- as.double(milo$patient_age)
    } else if (test_var == 'parity') {
      milo <- milo[,milo$tissue_condition %in% c('Mammoplasty WT')]
      milo <- milo[, milo$parity != 'unknown']
      milo$milo_test <- milo$parity != '0'
    } else {
      print(paste0('Incorrect input: test_var = ', test_var))
      next
    }
    #block_var meta
    if (block_var == 'none') {
      milo$milo_block <- 'none'
    } else if (block_var == 'patient_age') {
      milo$milo_block <- milo$patient_age
    } else if (block_var == 'parity') {
      milo <- milo[, milo$parity %in% c('0','1','2','3','4')]
      milo$milo_block[milo$parity %in% c('1','2','3','4')] <- 'Parous'
      milo$milo_block[milo$parity %in% c('0')] <- 'Nulliparous'
    } else if (block_var == 'parity_age') {
      milo <- milo[, milo$parity %in% c('0','1','2','3','4')]
      milo$milo_block1[milo$parity %in% c('1','2','3','4')] <- 'Parous'
      milo$milo_block1[milo$parity %in% c('0')] <- 'Nulliparous'
      milo$milo_block2 <- milo$patient_age
    }else {
      print(paste0('Incorrect input: block_var = ', block_var))
      next
    }
    
    da_out <- read.csv(paste0(prefix, '/milo_testing/', test_var, '/da-block_', block_var, '-all.csv'))
    
    plot_milo(milo = milo,
              da_out = da_out,
              plot_pwd = paste0(prefix, '/milo_plots/', test_var, '/'),
              block_var = block_var,
              sample_tp = 'all')
  }
}
