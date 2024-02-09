#Use celltypist logistic regression for classification of celltypes across datasets:

#conda env celltypist

#libraries
import scanpy as sc
import celltypist
from celltypist import models

import pandas as pd
import numpy as np
import os
import time

#load data
adata = sc.read('/nfs/research/marioni/kunz/HBCA/data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')
adata.obs['cellID'] = adata.obs.index

#get metadata
reed_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2023-06-19.csv')
reed_colData.index = reed_colData.cellID.values
gray_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/gray_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)
murrow_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/murrow_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)
nee_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/nee_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)
nee_extra_coldata = pd.read_table('/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/original/epimetadata.txt', sep=' ')
pal_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/pal_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)
twigger_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/twigger_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)
kumar_colData = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/kumar_phenodata.csv').set_index('Unnamed: 0').rename_axis(None)

#add metadata
adata.obs['celltype'] = 'not_mapped'

#reed
temp = adata.obs.loc[adata.obs.dataset=='HBCA',].merge(pd.DataFrame({'cellID': reed_colData.cellID.values,
                                                                 'celltype_new': reed_colData.level2.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='HBCA', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='HBCA'] = temp3.copy()

#nee
df = nee_colData.merge(pd.DataFrame({'cellID': nee_extra_coldata.index,
                                     'epithelial_state': nee_extra_coldata['epithelial_states']}),
                       on='cellID', how='left', sort=False)
df['celltype'] = df.level1
df.celltype[df.level1.isin(['Luminal1','Luminal2','Basal'])] = df.epithelial_state[df.level1.isin(['Luminal1','Luminal2','Basal'])].copy()
df.celltype[df.celltype.isna()] = 'Doublet'
temp = adata.obs.loc[adata.obs.dataset=='Nee',].merge(pd.DataFrame({'cellID': df.cellID.values,
                                                                    'celltype_new': df.celltype.values}),
                                                      how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Nee', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Nee'] = temp3.copy()

#kumar
temp = adata.obs.loc[adata.obs.dataset=='Kumar',].merge(pd.DataFrame({'cellID': kumar_colData.cellID.values,
                                                                 'celltype_new': kumar_colData.level2.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Kumar', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Kumar'] = temp3.copy()

#murrow
temp = adata.obs.loc[adata.obs.dataset=='Murrow',].merge(pd.DataFrame({'cellID': murrow_colData.cellID.values,
                                                                 'celltype_new': murrow_colData.level1.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Murrow', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Murrow'] = temp3.copy()

#pal
temp = adata.obs.loc[adata.obs.dataset=='Pal',].merge(pd.DataFrame({'cellID': pal_colData.cellID.values,
                                                                 'celltype_new': pal_colData.level1.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Pal', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Pal'] = temp3.copy()

#twigger
temp = adata.obs.loc[adata.obs.dataset=='Twigger',].merge(pd.DataFrame({'cellID': twigger_colData.cellID.values,
                                                                 'celltype_new': twigger_colData.level1.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Twigger', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Twigger'] = temp3.copy()

#gray
temp = adata.obs.loc[adata.obs.dataset=='Gray',].merge(pd.DataFrame({'cellID': gray_colData.cellID.values,
                                                                 'celltype_new': gray_colData.level2.values}),
                                                   how='left', on='cellID', sort=False).set_index('cellID')
temp2 = adata.obs.loc[adata.obs.dataset=='Gray', 'cellID'].values
temp3 = temp.loc[temp2, 'celltype_new'].values
adata.obs['celltype'].loc[adata.obs.dataset=='Gray'] = temp3.copy()

os.makedirs('/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/output_final/h5ad/')
os.makedirs('/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/output_final/colData/')
adata.write_h5ad('/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/output_final/h5ad/integrated_unlabelled.h5ad')
adata.obs.to_csv('/nfs/research/marioni/areed/projects/hbca/integrated_celltypes_compare/2023-06-21/scvi/celltypist/output_final/colData/integrated_unlabelled.csv')

##Function

def map_dataset_labels(adata, dataset, dataset_label):
    print(dataset_label)

    # preprocessing
    adata_ref = adata[adata.obs.dataset == dataset, ]
    adata_query = adata[adata.obs.dataset != dataset, ]

    # Sample 500 cells from each cell type for `adata_Elmentaite`.
    # All cells from a given cell type will be selected if the cell type size is < 2000.

    #This fails (?) so do manually
    # sampled_cell_index = celltypist.samples.downsample_adata(adata_ref, mode='each', n_cells=2000, by='celltype',
    #                                                          return_index=True)

    categories = adata_ref.obs['celltype'].unique()
    sampled_cell_index = {}
    num_rows = 2000
    for category in categories:
        subset = adata_ref.obs[adata_ref.obs['celltype'] == category]
        available_rows = len(subset)
        if available_rows <= num_rows:
            sampled_cell_index[category] = list(subset['cellID'])
        else:
            sampled_rows = subset.sample(n=num_rows, random_state=42)
            sampled_cell_index[category] = list(sampled_rows['cellID'])
    sampled_cell_index = np.concatenate(list(sampled_cell_index.values()))

    print(f"Number of downsampled cells for training: {len(sampled_cell_index)}")

    # Use `celltypist.train` to quickly train a rough CellTypist model.
    # You can also set `mini_batch = True` to enable mini-batch training.
    t_start = time.time()
    model_fs = celltypist.train(adata_ref[sampled_cell_index,], 'celltype', n_jobs=10, max_iter=10, use_SGD=True)
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")

    # select genes base on first quick model
    gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -300, axis=1)[:, -300:]
    gene_index = np.unique(gene_index)
    print(f"Number of genes selected: {len(gene_index)}")

    # properly train model
    # Add `check_expression = False` to bypass expression check with only a subset of genes.
    t_start = time.time()
    model = celltypist.train(adata_ref[sampled_cell_index, gene_index], 'celltype', check_expression=False, n_jobs=10,
                             max_iter=400)
    t_end = time.time()
    print(f"Time elapsed: {(t_end - t_start) / 60} minutes")

    # Save the model.
    model.write('output_final/celltypist_demo_folder/model_from_' + dataset_label + '.pkl')

    # CellTypist prediction without over-clustering and majority-voting.
    t_start = time.time()
    predictions = celltypist.annotate(adata_query, model='output_final/celltypist_demo_folder/model_from_' + dataset_label + '.pkl')
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")

    adata_pred = predictions.to_adata()
    adata_pred.write_h5ad('output_final/celltypist_demo_folder/adata_predictions_from_' + dataset_label + '.h5ad')

    # save predictions
    #predictions.predicted_labels.to_csv('output_final/celltypist_demo_folder/predictions_meta_from_' + dataset_label + '.csv')
    adata_pred.obs.to_csv('output_final/celltypist_demo_folder/predictions_meta_from_' + dataset_label + '.csv')


#perform label transfer
map_dataset_labels(adata, 'HBCA', 'reed')
map_dataset_labels(adata, 'Kumar', 'kumar')
map_dataset_labels(adata, 'Nee', 'nee')
map_dataset_labels(adata, 'Gray', 'gray')
map_dataset_labels(adata, 'Twigger', 'twigger')
map_dataset_labels(adata, 'Pal', 'pal')
map_dataset_labels(adata, 'Murrow', 'murrow')


