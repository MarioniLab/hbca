#Make summary bar plots of the metadata we have.

library(plyr)
library(dplyr)
library(ggplot2)




### Load data

#metadf <- read.csv('C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/metadata/2022-05-17_all_metadata.csv')
metadf <- read.csv('/nfs/research/marioni/areed/projects/hbca/metadata/2022-05-17_all_metadata.csv')
sce_meta <- read.csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
# metadf <- read.csv('/nfs/research/marioni/areed/projects/hbca/metadata/2022-05-17_all_metadata.csv')


#THIS IS WRONG!!!!! CHECK. IT LEADS TO A HRBR1 and HRBR2 being removed unnecessarily (ie these were in the sce object)
#I think at least
#metadf2 <- metadf[!(metadf$library_status %in% c('failed','very_low_cell_count')), ] #these samples were removed from downstream analysis for failing sequencing - NOT ALL!

#fixed using this
metadf <- metadf[metadf$sample_id %in% sce_meta$sampleID, ]


### Setup metadata by patient
columns_to_use <- c("patient_id", "parity", "patient_age", "tissue_condition", 
                    "Ethnicity", "Age_at_first_menstruation", 
                    "Menopause_status.Pre.peri.post.", 
                    "OCP_use.past.present.no.never.NA.", 
                    "HRT_use_.past.present.never.NA.", "Body_mass_index", 
                    "Smoking_status_.past.present.no.never.NA.", 
                    "Alcohol_consumption_levels")
patient_metadf <- metadf[,columns_to_use]
names(patient_metadf) <- c("patientID", "Parity", "Patient_age", "Tissue_condition", 
                           "Ethnicity", "Age_at_first_menstruation", 
                           "Menopause_status", 
                           "OCP_use", 
                           "HRT_use", "BMI", 
                           "Smoking_status", 
                           "Alcohol_consumption")
patient_metadf <- unique(patient_metadf) #make unique per patient (not sample)
#no longer needed
#patient_metadf <- patient_metadf[2:dim(patient_metadf)[1], ] #remove empty row (spike-in?)

#Map values to readable versions
patient_metadf$Ethnicity <- mapvalues(patient_metadf$Ethnicity, from = unique(patient_metadf$Ethnicity), to = c("White", 'NA', "Other", "Other", "Black", "White", "Asian", "Asian", "White"))
patient_metadf$Menopause_status <- mapvalues(patient_metadf$Menopause_status, from = unique(patient_metadf$Menopause_status), to = c("Pre", "Post", "Post (SI)", "NA"))
patient_metadf$OCP_use <- mapvalues(patient_metadf$OCP_use, from = unique(patient_metadf$OCP_use), to = c("Past", 'NA', 'Never', 'Past', 'Never'))
patient_metadf$HRT_use <- mapvalues(patient_metadf$HRT_use, from = unique(patient_metadf$HRT_use), to = c('Never', 'Present', 'Present', 'NA', 'Past'))

#Set up levels for plotting order
patient_metadf[is.na(patient_metadf)] <- 'NA'
patient_metadf[patient_metadf == 'unknown'] <- 'NA'
patient_metadf$Tissue_condition <- factor(patient_metadf$Tissue_condition, levels = c("Mammoplasty WT", "Mastectomy BRCA1", "Mastectomy BRCA2", "Mastectomy WT", "Mastectomy unknown", "Contralateral BRCA1"))
patient_metadf$Ethnicity <- factor(patient_metadf$Ethnicity, levels = c('White', 'Asian', 'Black', 'Other', 'NA'))
patient_metadf$Menopause_status <- factor(patient_metadf$Menopause_status, levels = c('Pre', 'Post', 'Post (SI)', 'NA'))
patient_metadf$OCP_use <- factor(patient_metadf$OCP_use, levels = c('Past', 'Never', 'NA'))
patient_metadf$HRT_use <- factor(patient_metadf$HRT_use, levels = c('Present', 'Past', 'Never', 'NA'))
patient_metadf$Smoking_status <- factor(patient_metadf$Smoking_status, levels = c('Past', 'Never', 'NA'))
patient_metadf$Alcohol_consumption <- factor(patient_metadf$Alcohol_consumption, levels = c('None', 'Low', 'Moderate', 'High', 'Socially', 'NA'))

#Check these
head(patient_metadf) 
dim(patient_metadf) #55, 12 #was 53, 12 but this was missing two patients


### Make barplots
#fig_dir <- 'C:/Users/44756/OneDrive - University of Cambridge/WTKLAB/Projects/hbca/codon/figures/figure1/barplot/'

fig_dir <- '/nfs/research/marioni/areed/projects/hbca/figures/figure1/barplot/'
dir.create(paste0(fig_dir), showWarnings = F, recursive = T)

#for categorical data
for (meta_col in c("Parity", "Tissue_condition", 
                   "Ethnicity", 
                   "Menopause_status", 
                   "OCP_use", 
                   "HRT_use",
                   "Smoking_status", 
                   "Alcohol_consumption")) {
  pdf(file = paste0(fig_dir, 'f1_metaboxplot_', meta_col, '.pdf'))
  print(ggplot(data=patient_metadf, aes_string(x=meta_col)) +
          geom_bar(fill='black') +
          xlab(gsub('_', ' ', meta_col)) +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
} 

#for continuous data
for (meta_col in c("Patient_age", "Age_at_first_menstruation", "BMI")) {
  patient_metadf[, meta_col] <- as.numeric(patient_metadf[, meta_col])
  #patient_metadf[is.na(patient_metadf[, meta_col]), meta_col] <- 'NA'
  
  pdf(file = paste0(fig_dir, 'f1_metaboxplot_', meta_col, '.pdf'))
  print(ggplot(data=patient_metadf, aes_string(x=meta_col)) +
          geom_histogram(color="White", fill='black', bins = 10) +
          xlab(gsub('_', ' ', meta_col)) +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()
} 







