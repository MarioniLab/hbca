#! Do milo analysis on HBCA subsets run without splitting ST's (sample_types - organoid/supernatant).
#This time set up to use the same neighborhoods for all AR vs HR tests. 
#scVI new


suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(optparse))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

# #Epi
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/initial/output/outs/epi'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/epi_2/output'

# #Str
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/initial/output/outs/str'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/str_2/output'

# #Imm
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/initial/output/outs/imm'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2022-04-05/scvi_new/imm_2/output'

## options
option_list = list(make_option(c('--input_path'),
                               type='character',
                               help='Pathway to retrieve outs data from.'),
                   make_option(c('--output_path'),
                               type='character',
                               help='Pathway to output folder.'))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


## functions
read_outs <- function(outs_pwd){
  
  #read in counts matrices
  counts_m <- t(readMM(paste0(outs_pwd, '/counts.mtx')))
  logcounts_m <-  t(readMM(paste0(outs_pwd, '/logcounts.mtx')))
  
  #metadata
  clData <- read.csv(paste0(outs_pwd, '/colData.csv')) 
  rownames(clData) <- clData$cellID
  rwData <- read.csv(paste0(outs_pwd, '/rowData.csv')) 
  rownames(rwData) <- rwData$gene_ids
  
  #reducedDims
  # pca <- read.csv(paste0(outs_pwd, '/pca.csv'))
  # pca <- subset(as.data.frame(pca), select=-c(X))
  # rownames(pca) <- clData$cellID
  scvi <- read.csv(paste0(outs_pwd, '/scVI.csv'))
  scvi <- subset(as.data.frame(scvi), select=-c(X))
  rownames(scvi) <- clData$cellID
  
  umap <- read.csv(paste0(outs_pwd, '/umap.csv'))
  umap <- subset(as.data.frame(umap), select=-c(X))
  rownames(umap) <- clData$cellID
  
  #create sce
  sce <- SingleCellExperiment(assays=SimpleList('counts' = counts_m, 
                                                'logcounts' = logcounts_m), 
                              rowData = rwData, 
                              colData = clData, 
                              reducedDims = list('scVI'=scvi, 'UMAP'=umap))
  
  return(sce)
}

setup_milo <- function(sce, subset_tf=FALSE, subset_category=NULL, 
                       subset_selection=NULL, k_param, d_param, prop, 
                       reduced_dimension, samples_var, out_nghd_size_pw, 
                       random_seed){
  set.seed(random_seed)
  
  if (!(subset_tf %in% c(FALSE, 'FALSE', 'False'))){
    sce <- sce[, colData(sce)[,subset_category] %in% subset_selection]
  }
  
  #Create Milo object
  milo <- Milo(sce)
  milo
  
  #Fix colData
  colData(milo)[, samples_var] <- factor(colData(milo)[, samples_var])
  
  #KNN graph
  milo <- buildGraph(milo, k = k_param, d = d_param, transposed = TRUE, 
                     reduced.dim = reduced_dimension)
  
  #Define representative neighbourhoods of KNN
  milo <- makeNhoods(milo, prop = prop, k = k_param, d=d_param, refined = TRUE, 
                     reduced_dims = reduced_dimension, refinement_scheme = 'graph') #use refinement_scheme = 'graph' for speed (devel branch specific)
  
  #plot distribution of neighbourhood sizes (want mean to be about 5x number of samples)
  dir.create(dirname(out_nghd_size_pw), showWarnings = FALSE, recursive = TRUE)
  pdf(out_nghd_size_pw, width=12, height = 8)
  print(plotNhoodSizeHist(milo))
  dev.off()
  
  #Counting cells in Nghds
  milo <- countCells(milo, meta.data = colData(milo), samples = samples_var)
  
  #Compute nghd connectivity #not required with refinement_scheme = 'graph'
  #milo <- calcNhoodDistance(milo, d=d_param, reduced.dim = reduced_dimension) #not needed as we use refinement_scheme="graph"
  
  return(milo)
}


## Analysis

#setup sce
print(opt$input_path)
print(opt$output_path)
sce <- read_outs(opt$input_path) 

#make sure to have level1 as celltype label
if (is.null(sce$level1)) {
  sce$level1 <- sce$celltype
}

df <- read.csv('/nfs/research/marioni/areed/projects/hbca/metadata/2022-05-17_all_metadata.csv')
df$X = NULL
df <- df[,c('sample_id', "Menopause_status.Pre.peri.post.", "Body_mass_index", 'Smoking_status_.past.present.no.never.NA.', "HRT_use_.past.present.never.NA.", "OCP_use.past.present.no.never.NA.")]
names(df) <- c('sampleID', "Menopause_status", "Body_mass_index", 'Smoking_status', "HRT_use", "OCP_use")
joint_clData_all <- left_join(as.data.frame(colData(sce)), df, by='sampleID')
new_cols <- c("Menopause_status", "Body_mass_index", 'Smoking_status', "HRT_use", "OCP_use")
for (col in new_cols) {
  colData(sce)[, col] <- joint_clData_all[, col]
}

metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
rownames(metadata_all) <- metadata_all$cellID
metadata_sub <- metadata_all[sce$cellID, c('level1','level2')]
sce$level1 <- metadata_sub$level1
sce$level2 <- metadata_sub$level2

dir.create(paste0(opt$output_path, '/input/'), showWarnings = FALSE, recursive = TRUE)
saveRDS(sce, paste0(opt$output_path, '/input/input_sce.rds')) 
# sce <- readRDS(paste0(opt$output_path, '/input/input_sce.rds'))

#remove bad quality cells/doublets
bad_celltypes <- c('Doublet')
sce <- sce[, !(sce$level2 %in% bad_celltypes)]

#remove FACS sorted samples - doesn't make sense to do comparisons of cell numbers from these.
sce <- sce[, !(sce$before %in% c('Organoid LP sorted', 'Supernatant live-sorted'))]

#This should only really be visible for the outlier samples which may not appear in all tests.

block_var_category <- 'parity'
print(block_var_category)

#subset sce as required
sce_sub <- sce[, sce$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')]
sce_sub <- sce_sub[, sce_sub$parity != 'unknown']

#make milo object
milo <- setup_milo(sce = sce_sub,
                   subset_tf = FALSE,
                   subset_category = 'None',
                   subset_selection = 'None',
                   k_param = 50,
                   d_param = 20,
                   prop = 0.3,
                   reduced_dimension = 'scVI',
                   samples_var = 'sampleID',
                   paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_setup/nghd_size_hist_epi_', block_var_category, '.pdf'),
                   random_seed = 422)

#save
dir.create(paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/'), showWarnings = FALSE, recursive = TRUE)
saveRDS(milo, file=paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/epi_nghd_', block_var_category, '.rds'))

print(paste0('Dimensions: ', dim(milo)))
























