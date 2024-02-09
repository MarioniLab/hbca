#python

#Prepare the Kumar et al dataset for joining and integration in iHBCA.

#conda env milo



#packages
import scanpy as sc
import pandas as pd
import os
import scipy.io as sio

#load data
adata = sc.read('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/original/kumar_adata.h5ad')

#make seperate counts, rowdata and coldata of set formatting
os.makedirs('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/', exist_ok=True)

#rowData
kumar_rowData = pd.DataFrame(data={'symbol': adata.var.feature_name.values,
                                   'gene_id': adata.var.index},
                             index=adata.var.feature_name.values)

print(kumar_rowData.head())
kumar_rowData.to_csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_features.csv')

#colData
kumar_meta = adata.obs

level0_dict = {'Fibro-major': 'Stroma', 'basal': 'Epithelial', 'Fibro-prematrix': 'Stroma', 'LummHR-major': 'Epithelial',
                  'Vas-venous': 'Stroma', 'Fibro-matrix': 'Stroma', 'Lumsec-basal': 'Epithelial',
                  'Vas-capillary': 'Stroma', 'vsmc': 'Stroma', 'pericytes': 'Stroma', 'Lumsec-major': 'Epithelial',
                  'CD8-Trm': 'Immune', 'LummHR-SCGB': 'Epithelial', 'CD4-Th-like': 'Immune', 'Vas-arterial': 'Stroma',
                  'LummHR-active': 'Epithelial', 'CD4-Th': 'Immune', 'Lymph-major': 'Stroma', 'Macro-m2': 'Immune',
                  'CD8-Tem': 'Immune', 'Macro-m1': 'Immune', 'plasma_IgA': 'Immune', 'CD4-naive': 'Immune', 'CD4-Tem': 'Immune',
                  'NK': 'Immune', 'Mono-classical': 'Immune', 'Lumsec-KIT': 'Epithelial', 'bmem_switched': 'Immune',
                  'cDC2': 'Immune', 'NKT': 'Immune', 'Fibro-SFRP4': 'Stroma', 'CD4-Treg': 'Immune', 'b_naive': 'Immune',
                  'CD4-activated': 'Immune', 'Lumsec-myo': 'Epithelial', 'NK-ILCs': 'Immune', 'Macro-m1-CCL': 'Immune',
                  'Macro-m2-CXCL': 'Immune', 'Lumsec-HLA': 'Epithelial', 'Macro-lipo': 'Immune', 'Lumsec-prol': 'Epithelial',
                  'plasma_IgG': 'Immune', 'Lumsec-lac': 'Epithelial', 'CD8-activated': 'Immune',
                  'Mono-non-classical': 'Immune', 'Mast': 'Immune', 'cDC1': 'Immune', 'GD': 'Immune', 'mDC': 'Immune',
                  'Macro-IFN': 'Immune', 'Neutrophil': 'Immune', 'bmem_unswitched': 'Immune', 'mye-prol': 'Immune',
                  'Lymph-valve1': 'Stroma', 'T_prol': 'Immune', 'Lymph-valve2': 'Stroma', 'Lymph-immune': 'Stroma',
                  'pDC': 'Immune'}
level1_dict = {'Fibro-major': 'Fibroblast', 'basal': 'Basal', 'Fibro-prematrix': 'Fibroblast', 'LummHR-major': 'Luminal Hormone Responsive',
                  'Vas-venous': 'Vascular Endothelial', 'Fibro-matrix': 'Fibroblast', 'Lumsec-basal': 'Luminal Secretory',
                  'Vas-capillary': 'Vascular Endothelial', 'vsmc': 'Perivascular', 'pericytes': 'Perivascular', 'Lumsec-major': 'Luminal Secretory',
                  'CD8-Trm': 'Lymphoid', 'LummHR-SCGB': 'Luminal Hormone Responsive', 'CD4-Th-like': 'Lymphoid', 'Vas-arterial': 'Vascular Endothelial',
                  'LummHR-active': 'Luminal Hormone Responsive', 'CD4-Th': 'Lymphoid', 'Lymph-major': 'Lymphatic Endothelial', 'Macro-m2': 'Myeloid',
                  'CD8-Tem': 'Lymphoid', 'Macro-m1': 'Myeloid', 'plasma_IgA': 'Lymphoid', 'CD4-naive': 'Lymphoid', 'CD4-Tem': 'Lymphoid',
                  'NK': 'Lymphoid', 'Mono-classical': 'Myeloid', 'Lumsec-KIT': 'Luminal Secretory', 'bmem_switched': 'Lymphoid',
                  'cDC2': 'Myeloid', 'NKT': 'Lymphoid', 'Fibro-SFRP4': 'Fibroblast', 'CD4-Treg': 'Lymphoid', 'b_naive': 'Lymphoid',
                  'CD4-activated': 'Lymphoid', 'Lumsec-myo': 'Luminal Secretory', 'NK-ILCs': 'Lymphoid', 'Macro-m1-CCL': 'Myeloid',
                  'Macro-m2-CXCL': 'Myeloid', 'Lumsec-HLA': 'Luminal Secretory', 'Macro-lipo': 'Myeloid', 'Lumsec-prol': 'Luminal Secretory',
                  'plasma_IgG': 'Lymphoid', 'Lumsec-lac': 'Luminal Secretory', 'CD8-activated': 'Lymphoid',
                  'Mono-non-classical': 'Myeloid', 'Mast': 'Myeloid', 'cDC1': 'Myeloid', 'GD': 'Lymphoid', 'mDC': 'Myeloid',
                  'Macro-IFN': 'Myeloid', 'Neutrophil': 'Myeloid', 'bmem_unswitched': 'Lymphoid', 'mye-prol': 'Myeloid',
                  'Lymph-valve1': 'Lymphatic Endothelial', 'T_prol': 'Lymphoid', 'Lymph-valve2': 'Lymphatic Endothelial', 'Lymph-immune': 'Lymphatic Endothelial',
                  'pDC': 'Myeloid'}

kumar_colData = pd.DataFrame(data={'cellID': kumar_meta.index,
                                   'patientID': kumar_meta.donor_id.values,
                                   'batch': kumar_meta.donor_id.values,
                                   'level0': kumar_meta.author_cell_type.values.map(level0_dict),
                                   'level1': kumar_meta.author_cell_type.values.map(level1_dict),
                                   'level2': kumar_meta.author_cell_type.values},
                             index=kumar_meta.index)

print(kumar_colData.head())

kumar_meta.to_csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_full_phenodata.csv')
kumar_colData.to_csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_phenodata.csv')

#counts
adata_raw = adata.raw.to_adata()

kumar_counts = adata_raw.X.T #transpose to match format of other saved datasets
sio.mmwrite("/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_counts.mtx", kumar_counts)
