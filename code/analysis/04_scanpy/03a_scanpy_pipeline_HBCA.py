#!/usr/bin/env python3

#####################
## Scanpy Pipeline ##
#####################

import numpy as np
import pandas as pd
import scanpy as sc
import scipy

import glob
import os
from collections import Counter

import scvi
import torch
device = torch.device("cuda")


def get_adata_10x(sampleIDs, path_dict):
    list_adata = []
    
    for sampleID in sampleIDs:
        print(sampleID)
        adata_temp = sc.read_10x_h5(path_dict['10x_dir'] + sampleID + path_dict['10x_fstruct'])
        adata_temp.var_names_make_unique()
        list_adata.append(adata_temp)
    
    # update to concat() in the future - https://anndata.readthedocs.io/en/latest/concatenation.html
    adata_merged = list_adata[0].concatenate(list_adata[1:],
                                             join='outer',
                                             fill_value=0,
                                             batch_key = 'sampleID',
                                             batch_categories = sampleIDs)
    
    adata_merged.obs['cellID'] = [cellID.replace('-1-', '-') for cellID in adata_merged.obs.index.tolist()]
    adata_merged.obs = adata_merged.obs.set_index('cellID')
    
    return(adata_merged)


def get_adata(path_dict):
    
    files = glob.glob(path_dict['10x_dir']+'*'+path_dict['10x_fstruct'])
    sampleIDs = [os.path.basename(os.path.dirname(samplepath.replace('/outs', ''))) for samplepath in files]
    
    adata = get_adata_10x(sampleIDs, path_dict)
    
    # import metadata
    if(path_dict['scrublet'] is not None):
        # filter scrublet doublet calling if passed along
        scrublet_scores = pd.read_csv(path_dict['scrublet'], index_col=0)
        cellIDs_scrublet = scrublet_scores[~scrublet_scores['is_doublet']].index.values
        
        if 'retain_scrublet' in path_dict.keys():
            cellIDs_scrublet = [*cellIDs_scrublet, *path_dict['retain_scrublet']]
        
        cellIDs_pass = list(set.intersection(set(cellIDs_scrublet), set(adata.obs.index.values)))
        adata = adata[cellIDs_pass, :].copy()
    
    # filter cells based on QC
    metadata_qc = pd.read_csv(path_dict['metadata'], low_memory=False)
    metadata_qc = metadata_qc.set_index('cellID')
    cellIDs_pass = list(set.intersection(set(metadata_qc.index.values), set(adata.obs.index.values)))
    adata = adata[cellIDs_pass, :].copy()
    
    obs_merged = pd.merge(adata.obs, metadata_qc.drop(columns=['sampleID']), how='left', on='cellID')
    adata.obs = obs_merged
    
    # not really needed, but generates n_genes slot
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)
    
    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    # add the total counts per cell as observations-annotation to adata
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    adata.layers['raw'] = adata.X.copy()
    
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    
    adata.layers['lognorm'] = adata.X.copy()
    
    # HVGs
    sc.pp.highly_variable_genes(adata, layer='lognorm', n_top_genes=5000)
    
    return adata


def run_scVI(adata, batchID='batch', n_dims=20):
    # subset to HVGs
    train_adata = adata[:, adata.var['highly_variable']].copy()
    scvi.model.SCVI.setup_anndata(train_adata, layer='raw', batch_key=batchID)
    
    arches_params = {
        "use_layer_norm": "both",
        "use_batch_norm": "none",
        "encode_covariates": True,
        "dropout_rate": 0.2,
        "n_layers": 2,
    }
    
    vae = scvi.model.SCVI(train_adata, n_latent=n_dims, gene_likelihood="nb", **arches_params)
    # vae.train() # 
    vae.train(early_stopping=True,
        train_size=0.9,
        early_stopping_patience=45,
        max_epochs=400, 
        batch_size=1024, 
        limit_train_batches=20)
    
    adata.obsm['X_scVI'] = vae.get_latent_representation()
    
    return {'adata': adata, 'vae': vae}


def post_dimred_processing(adata, dimred='X_pca', run_clustering=True):
    # neighbourhood graph
    print("knn graph")
    sc.pp.neighbors(adata, use_rep=dimred, n_neighbors=15)
    # UMAP
    print("UMAP")
    sc.tl.umap(adata)
    
    # Leiden clustering
    if run_clustering:
        list_leiden_res = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        print("Leiden clustering")
        for leiden_res in list_leiden_res:
            print(leiden_res)
            sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))
    
    print("Diffusion map")
    sc.tl.diffmap(adata)
    
    return adata



def subset_adata_batch(adata, cellIDs_sub=[], batchID="batchID"):
    adata_sub = adata[cellIDs_sub,:].copy()
    print("Remove previous cluster labels")
    # remove previous leiden clustering
    filter_meta_cols = [not x.startswith(('leiden_')) for x in adata_sub.obs.columns]
    adata_sub.obs = adata_sub.obs.loc[:, filter_meta_cols]
    print("Remove batches with less than 3 cells")
    # batch correction needs at least 3 cells per batch:
    dict_donorID_count = Counter(adata_sub.obs[batchID])
    newDict = dict()
    # identify donors with less than 3 cells.
    for (key, value) in dict_donorID_count.items():
        if value < 3:
            newDict[key] = value
    
    filter_sub = [x not in newDict.keys() for x in adata_sub.obs[batchID]]
    adata_sub = adata_sub[filter_sub,:].copy()
    print("Rerun HVG detection")
    sc.pp.highly_variable_genes(adata_sub, layer='lognorm', n_top_genes=5000)
    
    return adata_sub


######################
## From raw samples ##
######################

# srun -t 7-0 -c 10 --mem=300G --gres=gpu:v100:1 --pty singularity shell --nv /hps/software/users/marioni/kunz/scAnalysis_v0.2.3.sif 


# no scrublet

path_dict = {'10x_dir': '/nfs/research/marioni/asteif/projects/hbca/cellranger/count/2021-07-10/',
             '10x_fstruct': '/outs/filtered_feature_bc_matrix.h5',
             'scrublet': None,
             'metadata': 'data/metadata_HBCA_QC_2022-08-24.csv.gz'}

adata = get_adata(path_dict)
adata.write_h5ad('data/HBCA_noscrublet_2023-06-19.h5ad')




adata = sc.read_h5ad('data/HBCA_noscrublet_2023-06-19.h5ad')
out_scVI = run_scVI(adata, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_2022-11-01.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_2022-11-01", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_noscrublet_scVI_processing_date_2023-06-19.h5ad')




######################
# Retain macrophages #
######################

adata_noscrublet = sc.read_h5ad("data/HBCA_noscrublet_scVI_processing_date_2023-06-19.h5ad")
cellIDs_retain_scrublet = adata_noscrublet.obs.index.values[adata_noscrublet.obs['leiden_1_0'] == '24']

path_dict = {'10x_dir': '/nfs/research/marioni/asteif/projects/hbca/cellranger/count/2021-07-10/',
             '10x_fstruct': '/outs/filtered_feature_bc_matrix.h5',
             'scrublet': 'data/scrublet-scores/scrublet_scores_global.csv.gz',
             'retain_scrublet': cellIDs_retain_scrublet,
             'metadata': 'data/metadata_HBCA_QC_2022-08-24.csv.gz'}

adata = get_adata(path_dict)
adata.write_h5ad('data/HBCA_2023-06-19.h5ad')




adata = sc.read_h5ad('data/HBCA_2023-06-19.h5ad')
out_scVI = run_scVI(adata, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_2022-11-01.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_2022-11-01", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_2023-06-19.h5ad')




###################
## Subset Immune ##
###################

adata = sc.read_h5ad('data/HBCA_scVI_processing_date_2023-06-19.h5ad')
# issue potentially related to
# https://github.com/scverse/scanpy/issues/1333
# workaround: (need to check properly and potentially submit bugreport)
adata.uns['log1p']['base'] = None
filter_sub = [x in ['Immune'] for x in adata.obs['level1']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_immune_2022-11-01.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_immune_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_2023-06-19.h5ad')



## remove doublets & debris
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_immune_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None
filter_sub = [x not in ['6', '8', '10'] for x in adata.obs['leiden_0_4']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')



## Subgroup

# # Immune -> T-cells
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['0'] for x in adata.obs['leiden_0_05']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_T-cells_2023-06-19.h5ad')


# Immune -> B-cells & plasma cells
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['2', '7'] for x in adata.obs['leiden_0_5']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_B-cells_plasma_2023-06-19.h5ad')


# Immune -> Myeloid
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['3'] for x in adata.obs['leiden_0_5']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_myeloid_2023-06-19.h5ad')



# # Immune -> CD4T and INFG+ T
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['1', '4', '8'] for x in adata.obs['leiden_0_4']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_CD4T_INFGT_2023-06-19.h5ad')




#######################
## Subset Epithelial ##
#######################

adata = sc.read_h5ad('data/HBCA_scVI_processing_date_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['Luminal progenitor', 'Luminal hormone sensing', 'Basal'] for x in adata.obs['level1']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_epithelial_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_epithelial_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_2023-06-19.h5ad')



# ## remove doublets & debris
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None
filter_sub = [x not in ['12'] for x in adata.obs['leiden_1_0']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.h5ad')




## Subgroup

# Epithelial -> Basal
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['2'] for x in adata.obs['leiden_0_3']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_BSL_2023-06-19.h5ad')



# Epithelial -> LP
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['0', '4', '6'] for x in adata.obs['leiden_0_3']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LP_2023-06-19.h5ad')



# Epithelial -> LHS
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['1', '3', '5'] for x in adata.obs['leiden_0_3']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LHS_2023-06-19.h5ad')






###################
## Subset Stroma ##
###################

adata = sc.read_h5ad('data/HBCA_scVI_processing_date_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['Fibroblast', 'Vascular mural', 'Endothelial [vascular]', 'Endothelial [lymphatic]'] for x in adata.obs['level1']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_stroma_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_stroma_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_2023-06-19.h5ad')



# ## remove doublets & debris
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_stroma_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x not in ['5'] for x in adata.obs['leiden_0_2']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.h5ad')



# Stroma -> vascular mural
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['2'] for x in adata.obs['leiden_0_05']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_VM_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_VM_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_VM_2023-06-19.h5ad')


# Stroma -> Fibroblast
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['0'] for x in adata.obs['leiden_0_05']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_FB_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_FB_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_FB_2023-06-19.h5ad')



# Stroma -> Endothelial
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x in ['1', '3'] for x in adata.obs['leiden_0_05']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_2023-06-19.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_2023-06-19", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_2023-06-19.h5ad')



## remove doublets & debris
adata = sc.read_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_2023-06-19.h5ad')
adata.uns['log1p']['base'] = None

filter_sub = [x not in ['4'] for x in adata.obs['leiden_0_2']]
cellIDs_sub = adata[filter_sub,:].obs_names.values
adata_sub = subset_adata_batch(adata, cellIDs_sub=cellIDs_sub, batchID="processing_date")

out_scVI = run_scVI(adata_sub, batchID="processing_date")
adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI')
adata_scVI.write_h5ad('data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_cleaned_2023-06-19.h5ad')