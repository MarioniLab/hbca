#!/usr/bin/env python3

############################
## adata from raw samples ##
############################

import anndata as ad
import pandas as pd

#  pal_data  twigger_data

def get_adata_Gray():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/gray_data/formatted/"
    # data_dir = "datasets/gray_data/formatted/"
    adata = ad.read_mtx(data_dir+'gray_counts.mtx').T
    metadata = pd.read_csv(data_dir+'gray_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'gray_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)


def get_adata_Murrow():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/murrow_data/formatted/"
    # data_dir = "datasets/murrow_data/formatted/"
    adata = ad.read_mtx(data_dir+'murrow_counts.mtx').T
    metadata = pd.read_csv(data_dir+'murrow_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'murrow_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    # adata.obs['batch'] = adata.obs['Batch']
    # adata.obs['patientID'] = adata.obs['Sample']
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)


def get_adata_Nee():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/nee_data/formatted/"
    # data_dir = "datasets/nee_data/formatted/"
    adata = ad.read_mtx(data_dir+'nee_counts.mtx').T
    metadata = pd.read_csv(data_dir+'nee_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'nee_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)


def get_adata_Pal():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/pal_data/formatted/"
    # data_dir = "datasets/pal_data/formatted/"
    adata = ad.read_mtx(data_dir+'pal_counts.mtx').T
    metadata = pd.read_csv(data_dir+'pal_phenodata.csv', index_col=0, low_memory=False)
    genemap = pd.read_csv(data_dir+'pal_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)



def get_adata_Twigger():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/twigger_data/formatted/"
    # data_dir = "datasets/twigger_data/formatted/"
    adata = ad.read_mtx(data_dir+'twigger_counts.mtx').T
    metadata = pd.read_csv(data_dir+'twigger_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'twigger_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)



def get_adata_Kumar():
    data_dir = "/nfs/research/marioni/areed/projects/hbca/datasets/kumar_data/formatted/"
    # data_dir = "datasets/twigger_data/formatted/"
    adata = ad.read_mtx(data_dir+'kumar_counts.mtx').T
    metadata = pd.read_csv(data_dir+'kumar_phenodata.csv', index_col=0)
    genemap = pd.read_csv(data_dir+'kumar_features.csv', index_col=0)
    adata.obs = metadata
    adata.var = genemap
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    
    adata.layers['raw'] = adata.X.copy()
    
    return(adata)







adata_gray = get_adata_Gray()
adata_murrow = get_adata_Murrow()
adata_nee = get_adata_Nee()
adata_pal = get_adata_Pal()
adata_twigger = get_adata_Twigger()
adata_kumar = get_adata_Kumar()


adata_HBCA = ad.read_h5ad('data/HBCA_scVI_processing_date_2023-06-19.h5ad')
adata_HBCA.obs['batch'] = adata_HBCA.obs['processing_date']
adata_HBCA.var['gene_id'] = adata_HBCA.var['gene_ids']
adata_HBCA.X = adata_HBCA.layers['raw'].copy()



adatas = {
    "Gray": adata_gray,
    "Murrow": adata_murrow,
    "Nee": adata_nee,
    "Pal": adata_pal,
    "Twigger": adata_twigger,
    "Kumar": adata_kumar,
    "HBCA": adata_HBCA
}



##################################
# use ensembl geneIDs for merge
##################################

def filter_ENSG(adata):
    # filter_geneID = [isinstance(geneID, str) for geneID in adata.var['gene_id'].values]
    adata.var['gene_id'] = adata.var['gene_id'].astype(str)
    filter_geneID = [geneID != 'nan' for geneID in adata.var['gene_id'].values]
    adata = adata[:,filter_geneID].copy()
    adata.var_names = adata.var['gene_id'].values
    adata.var_names_make_unique()
    return(adata)




adata_gray = filter_ENSG(adata_gray).copy()
adata_murrow = filter_ENSG(adata_murrow).copy()
adata_nee = filter_ENSG(adata_nee).copy()
adata_pal = filter_ENSG(adata_pal).copy()
adata_twigger = filter_ENSG(adata_twigger).copy()
adata_kumar = filter_ENSG(adata_kumar).copy()
adata_HBCA = filter_ENSG(adata_HBCA).copy()


adatas = {
    "Gray": adata_gray,
    "Murrow": adata_murrow,
    "Nee": adata_nee,
    "Pal": adata_pal,
    "Twigger": adata_twigger,
    "Kumar": adata_kumar,
    "HBCA": adata_HBCA
}


adata_merged_inner = ad.concat(adatas, label="dataset", join="inner")

adata_merged_inner.write_h5ad('data/integration_HBCA_inner_ENSEMBL.h5ad')



#################
## Integration ##
#################


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




# srun -t 7-0 -c 10 --mem=500G --gres=gpu:v100:1 --pty singularity shell --nv /hps/software/users/marioni/kunz/scAnalysis_v0.2.3.sif 



adata = sc.read_h5ad('data/integration_HBCA_inner_ENSEMBL.h5ad')


batchID = "batch"
adata.layers['raw'] = adata.X.copy()

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.layers['lognorm'] = adata.X.copy()


print("Remove batches with less than 3 cells")
# batch correction needs at least 3 cells per batch:
dict_donorID_count = Counter(adata.obs[batchID])
newDict = dict()
# identify donors with less than 3 cells.
for (key, value) in dict_donorID_count.items():
    if value < 3:
        newDict[key] = value

filter_sub = [x not in newDict.keys() for x in adata.obs[batchID]]
adata = adata[filter_sub,:].copy()
print("Rerun HVG detection")
sc.pp.highly_variable_genes(adata, layer='lognorm', n_top_genes=5000)











out_scVI = run_scVI(adata, batchID="batch")
# export model
pd.DataFrame.from_dict({k:v[k] for k,v in out_scVI['vae'].history.items()}).to_csv("data/vae_history_HBCA_scVI_processing_date_2022-11-01.csv")
out_scVI['vae'].save("data/vae_HBCA_scVI_processing_date_2022-11-01", overwrite=True)

adata_scVI = out_scVI['adata'].copy()
adata_scVI = post_dimred_processing(adata_scVI, dimred='X_scVI', run_clustering=False)
adata_scVI.write_h5ad('data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')




####################################
## Add clustering and recalc UMAP ##
####################################

import numpy as np
import pandas as pd
import scanpy as sc


adata = sc.read_h5ad('data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')


# Leiden clustering
list_leiden_res = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
print("Leiden clustering")
for leiden_res in list_leiden_res:
    print(leiden_res)
    sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))

adata.write_h5ad('data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')



adata = sc.read_h5ad('data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')
# sc.tl.umap(adata, random_state=1)
adata.obsm['X_umap_seed_0'] = adata.obsm['X_umap'].copy()
sc.tl.umap(adata, random_state=1)
adata.obsm['X_umap_seed_1'] = adata.obsm['X_umap'].copy()
sc.tl.umap(adata, random_state=2)
adata.obsm['X_umap_seed_2'] = adata.obsm['X_umap'].copy()
adata.obsm['X_umap'] = adata.obsm['X_umap_seed_0'].copy()
adata.write_h5ad('data/integration_HBCA_inner_ENSEMBL_scVI.h5ad')


