#R

#Do milo analysis on HBCA subsets run without splitting sample_types (organoid/supernatant).

#Use the same neighborhoods for all AR vs HR tests
#scVI embedding


#libraries
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(miloR))
suppressMessages(library(Matrix))

suppressMessages(library(optparse))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

# #example options (testing)

# #Epi
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/epi/sce_epi_annotated.rds'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/epi/output'

# #Str
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/str/sce_str_annotated.rds'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/str/output'

# #Imm
# opt <- list()
# opt$input_path <- '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/initial/output/sce/imm/sce_imm_annotated.rds'
# opt$output_path <- '/nfs/research/marioni/areed/projects/hbca/milo/2023-06-21/scvi/imm/output'

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

sce <- readRDS(opt$input_path)

#remove stripped nulcei/doublets and donor specific clusters we do not want to include in analysis
bad_celltypes <- c('Doublet', 'stripped_nuclei', 'DDC1', 'DDC2')
sce <- sce[, !(sce$level2 %in% bad_celltypes)]

#remove FACS sorted samples - it doesn't make sense to do comparisons of cell numbers from these.
sce <- sce[, !(sce$before %in% c('Organoid LP sorted', 'Supernatant live-sorted'))]




#For testing AR donors (age and parity comparisons)
print('AR')
block_var_category <- 'parity' #naming convention - remnant of old pipeline
print(block_var_category)

#subset sce as required
sce_sub <- sce[, sce$tissue_condition %in% c('Mammoplasty WT')]
sce_sub <- sce_sub[, sce_sub$parity != 'unknown']

#make milo object
milo <- setup_milo(sce = sce_sub,
                   subset_tf = FALSE,
                   subset_category = 'None',
                   subset_selection = 'None',
                   k_param = 50, #optimised previously
                   d_param = 20,
                   prop = 0.3, #optimised previously
                   reduced_dimension = 'scVI',
                   samples_var = 'sampleID',
                   paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_setup/nghd_size_hist_epi_', block_var_category, '.pdf'),
                   random_seed = 422)

#save
dir.create(paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/'), showWarnings = FALSE, recursive = TRUE)
saveRDS(milo, file=paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/AR_nghd_', block_var_category, '.rds'))






#For testing HR donors (AR/HR-BR1 and AR/HR-BR2 comparisons)
#Ensure we create the milo object such that neighbourhoods are shared across 
#each of the comparisons (AR/HR-BR1 and AR/HR-BR2)

print('HR')
block_var_category <- 'parity' #naming convention - remnant of old pipeline
print(block_var_category)

#subset sce as required
sce_sub <- sce[, sce$tissue_condition %in% c('Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2')]
sce_sub <- sce_sub[, sce_sub$parity != 'unknown']

#make milo object
milo <- setup_milo(sce = sce_sub,
                   subset_tf = FALSE,
                   subset_category = 'None',
                   subset_selection = 'None',
                   k_param = 50, #optimised previously
                   d_param = 20,
                   prop = 0.3, #optimised previously
                   reduced_dimension = 'scVI',
                   samples_var = 'sampleID',
                   paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_setup/nghd_size_hist_epi_', block_var_category, '.pdf'),
                   random_seed = 422)

#save
dir.create(paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/'), showWarnings = FALSE, recursive = TRUE)
saveRDS(milo, file=paste0(opt$output_path, '/prop', 0.3, 'k', 50, '/milo_objects/HR_nghd_', block_var_category, '.rds'))
