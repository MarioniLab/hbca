#Regular NB differential abundance testing in the integrated dataset


####A few things to note:

### I need to remove anything strongly FACs sorted (Blocking term if only live sorted)
### There is no kind of normalisation per sample - not sure if the model can account for this...

#model

#Age:       ~ dataset + live_sorted_boolean + sample_type + parous_boolean + age
#parous:    ~ dataset + live_sorted_boolean + sample_type + age + parous_boolean

#librarys
library(plyr)
library(dplyr)
library(scran)
library(edgeR)
library(scater)
library(ggplot2)
library(EnhancedVolcano)
library(stringr)

#load data
metadata_all <- read.csv('/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/output_final/colData/integrated_labelled.csv') # read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/integration/initial/output/outs/all/colData.csv')

#collect extra colData

#Reed
reed_extra_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2023-06-19.csv') #read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
reed_extra_meta$live_sorted_boolean <- mapvalues(reed_extra_meta$before, 
                                                 from=c('Organoid LP sorted', 'Organoid unsorted', 'Supernatant live-sorted', 'Supernatant unsorted'), 
                                                 to = c('remove', 'no_sort', 'live_sorted', 'no_sort'))
reed_extra_meta$sample_type <- reed_extra_meta$sample_type_coarse
reed_extra_meta$risk_status <- mapvalues(reed_extra_meta$tissue_condition, 
                                         from=c('Contralateral BRCA1', 'Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2',
                                               'Mastectomy unknown', 'Mastectomy WT'),
                                         to=c('HR-cBR1', 'AR', 'HR-BR1', 'HR-BR2', 'HR-Unk', 'HR-Unk'))
reed_extra_meta$tissue_origin <- 'frozen'
reed_extra_meta <- reed_extra_meta[, c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')]

#Kumar
kumar_patient_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/original/extra_metadata.csv') ##/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_full_phenodata.csv
kumar_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_phenodata.csv')
kumar_patient_meta$risk_status <- mapvalues(kumar_patient_meta$tissue_source,
                                            from = c('Cancer Mastectomy', 'Prophylatic Mastectomy', 'Reduction Mammoplasty'),
                                            to = c('HR-cUnk', 'HR-Unk', 'AR'))
kumar_extra_meta <- merge(kumar_meta, unique(kumar_patient_meta[,c('patientID', 'age', 'parity', 'risk_status')]), by='patientID', all.x=T)
kumar_extra_meta$live_sorted_boolean <- 'live_sorted'
kumar_extra_meta$sample_type <- 'mixed'
kumar_extra_meta <- kumar_extra_meta[, c('cellID', 'live_sorted_boolean', 'sample_type', 'age', 'parity', 'risk_status')]
kumar_extra_meta$tissue_origin <- 'fresh'
names(kumar_extra_meta) <- c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')

#Nee
nee_extra_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_full_phenodata.csv')
nee_extra_meta$live_sorted_boolean <- 'remove' #all cells FACs sorted
nee_extra_meta$sample_type <- 'mixed'
nee_extra_meta$risk_status <- mapvalues(nee_extra_meta$tissue, 
                                        from = c('Contralateral', 'Prophylactic Mastectomy', 'Prophylatctic Mastectomy',
                                                 'Reduction Mammoplasty'),
                                        to = c('HR-cBR1', 'HR-BR1', 'HR-BR1', 'AR'))
nee_extra_meta <- nee_extra_meta[, c('X', 'live_sorted_boolean', 'sample_type', 'age', 'parity', 'risk_status')]
nee_extra_meta$tissue_origin <- 'frozen'
names(nee_extra_meta) <- c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')

#Gray 
gray_patient_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/original/patient_meta.csv') 
gray_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/gray_phenodata.csv')
gray_patient_meta$tissue_condition <- paste0(gray_patient_meta$Surgery, ' ', word(gray_patient_meta$Genotype..mutation., 1))
gray_patient_meta$risk_status <- mapvalues(gray_patient_meta$tissue_condition, 
                                           from=c('Contralateral prophylactic mastectomy BRCA1',
                                                  'Contralateral prophylactic mastectomy BRCA2',
                                                  'Prophylactic mastectomy BRCA1',
                                                  'Prophylactic mastectomy BRCA2',
                                                  'Reductive mammoplasty RAD51C',
                                                  'Reductive mammoplasty WT'),
                                           to=c('HR-cBR1', 'HR-cBR2', 'HR-BR1', 'HR-BR2', 'HR-RAD', 'AR'))
gray_patient_meta_sub <- gray_patient_meta[,c('scRNA.seq', 'Age', 'Births', 'risk_status')]
names(gray_patient_meta_sub) <- c('patientID', 'patient_age', 'parity', 'risk_status')
gray_extra_meta <- merge(gray_meta, unique(gray_patient_meta_sub), by='patientID', all.x=T, sort=F)
gray_extra_meta$live_sorted_boolean <- 'no_sort'
gray_extra_meta$sample_type <- 'mixed'
gray_extra_meta$tissue_origin <- 'fresh'
gray_extra_meta <- gray_extra_meta[, c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')]


#Twigger
twigger_patient_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/original/patient_meta.csv') 
twigger_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/twigger_phenodata.csv')
twigger_meta$patientID <- mapvalues(twigger_meta$patientID, 
                                    from=c('HMC1', 'HMC2', 'HMC2B', 'HMC3', 
                                           'HMC4', 'HMC5', 'HMC6', 'HMC7', 
                                           'HMC8',  'HMC9', 
                                           'RB1', 'RB2', 'RB3', 'RB4', 'RB5', 
                                           'RB6', 'RB7', 'RB8'),
                                    to=c('LMC1', 'LMC2', 'LMC2B', 'LMC3', 
                                         'LMC4', 'LMC5', 'LMC6', 'LMC7', 
                                         'LMC8',  'LMC9', 
                                         'NMC1', 'NMC2', 'NMC3', 'NMC4', 'NMC5', 
                                         'NMC6', 'NMC7', 'NMC1B'))
twigger_patient_meta_sub <- twigger_patient_meta[,c('patientID', 'patient_age', 'parity')]
twigger_extra_meta <- merge(twigger_meta, unique(twigger_patient_meta_sub), by='patientID', all.x=T, sort=F)
twigger_extra_meta$live_sorted_boolean <- 'no_sort'
twigger_extra_meta$sample_type <- 'mixed'
twigger_extra_meta$risk_status <- 'AR'
twigger_extra_meta$tissue_origin <- 'frozen'
twigger_extra_meta <- twigger_extra_meta[, c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')]

#Murrow
murrow_extra_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_full_phenodata.csv')
murrow_extra_meta$live_sorted_boolean <- mapvalues(murrow_extra_meta$Sort, 
                                                   from=c('Basal', 'Epithelial', 'Live_singlet', 'Luminal'), 
                                                   to = c('remove', 'remove', 'live_sorted', 'remove'))
murrow_extra_meta$sample_type <- 'mixed'
murrow_extra_meta$risk_status <- 'AR'
murrow_extra_meta$tissue_origin <- 'frozen'
murrow_extra_meta <- murrow_extra_meta[, c('X', 'live_sorted_boolean', 'sample_type', 'Age', 'Parity', 'risk_status', 'tissue_origin')]
names(murrow_extra_meta) <- c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')


#Pal
pal_patient_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/original/patient_meta.csv') 
pal_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_full_phenodata.csv')
pal_meta$patientID <- sapply(strsplit(pal_meta$GEO_sampleID, '_'), function(x) x[1])
pal_patient_meta$group <- gsub('-', '_', pal_patient_meta$Sample.Name)
pal_patient_meta$patient_age <- 'Unknown'
pal_patient_meta_sub <- pal_patient_meta[,c('group', 'patient_age', 'Parity', 'Condition')]
pal_extra_meta <- merge(pal_meta, unique(pal_patient_meta_sub), by='group', all.x=T, sort=F)
pal_extra_meta$live_sorted_boolean <- 'Unknown'
pal_extra_meta$sample_type <- 'Unknown'
pal_extra_meta$risk_status <- mapvalues(pal_extra_meta$Condition,
                                        from=c('Normal', 'Normal BRCA1+/- pre-neoplastic'),
                                        to=c('AR', 'HR-BR1'))
pal_extra_meta$tissue_origin <- 'fresh'
pal_extra_meta <- pal_extra_meta[, c('GEO_cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'Parity', 'risk_status', 'tissue_origin')]
names(pal_extra_meta) <- c('cellID', 'live_sorted_boolean', 'sample_type', 'patient_age', 'parity', 'risk_status', 'tissue_origin')



###Join these all together
extra_meta_list <- list('reed'=reed_extra_meta, 'kumar'=kumar_extra_meta,
                        'nee'=nee_extra_meta, 'gray'=gray_extra_meta,
                        'twigger'=twigger_extra_meta, 'pal'=pal_extra_meta,
                        'murrow'=murrow_extra_meta)

#explore and compare the extrametadata looking for any inconsistencies
for (dataset in names(extra_meta_list)) {
  print(dataset)
  print(head(extra_meta_list[[dataset]]))
}

#merge extrameta
all_extra_meta <- do.call(rbind, extra_meta_list)

# names(metadata_all)[2] <- 'X'
metadata_all_extra <- merge(metadata_all, all_extra_meta, by='cellID', all.X=T, all.Y=F, sort=F)

#One Pal sample has both AR and HR-BR1 cells - change to all HR-BR1
metadata_all_extra$risk_status[metadata_all_extra$patientID == 'MH0023'] <- 'HR-BR1'

#save
dir.create('/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/colData/', showWarnings = F, recursive = T)
write.csv(metadata_all_extra, '/nfs/research/marioni/areed/projects/hbca/integrated_da/2023-06-21/scvi/output_final/colData/integrated_datasets_with_extra_meta.csv')
































