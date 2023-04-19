'''Prepare (loompy) files for SCENIC'''

import pandas as pd
import scanpy as sc
import loompy as lp
import os
import numpy as np

output_dir = '/nfs/research/marioni/areed/projects/hbca/scenic/2022-04-05/scvi_new/hs/output'

#Load data
adata = sc.read_h5ad('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2022-11-01.h5ad')
metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_all = metadata_all.set_index('cellID')
adata.obs = metadata_all[metadata_all.index.isin(adata.obs.index)]
sc.pp.filter_genes(adata, min_cells=5)

#get interesting celltypes
celltypes_to_study = ['HS1', 'HS2', 'HS3', 'HS4']
adata = adata[adata.obs.level2.isin(celltypes_to_study), :]

# save full loom
f_loom_path_scenic = output_dir + "/loom/loom_hs_all.loom"
row_attrs = {"Gene": np.array(adata.var_names)}
col_attrs = {
    "CellID": np.array(adata.obs_names),
    "nGene": np.array(np.sum(adata.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata.X.transpose(), axis=0)).flatten(),
}
os.makedirs(output_dir + '/loom/', exist_ok=True)
lp.create(f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)

#save a subsetted version (n=1000)
np.random.seed(seed=42)
rnd_cells = np.random.choice(np.arange(1, adata.n_obs), size=1000, replace=False)
adata_sub = adata[rnd_cells]
# create basic row and column attributes for the loom file:
f_loom_path_scenic = output_dir + "/loom/loom_hs_1k.loom"
row_attrs = {"Gene": np.array(adata_sub.var_names)}
col_attrs = {
    "CellID": np.array(adata_sub.obs_names),
    "nGene": np.array(np.sum(adata_sub.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata_sub.X.transpose(), axis=0)).flatten(),
}
lp.create(f_loom_path_scenic, adata_sub.X.transpose(), row_attrs, col_attrs)


#save another subsetted version (n=4*9000 + #(LP proliferating, Other epithelial (1), Other epithelial (2)))
#this time make sure we have 3000 AR, 3000 HR-BR1 and 3000 HR-BR2
np.random.seed(seed=42)
cells_selected_total = []
for celltype_select in celltypes_to_study:
    print(celltype_select)
    if celltype_select in ['HS1', 'HS2', 'HS3', 'HS4']:
        for tissue_type in ['Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2']:
            cells_of_type = adata.obs.index.values[(adata.obs.level2 == celltype_select) & (adata.obs.tissue_condition == tissue_type)]
            cells_selected = cells_of_type[np.random.choice(np.arange(0, len(cells_of_type)),
                                                            size=min(3000, len(cells_of_type)), replace=False)]
            cells_selected_total = cells_selected_total + cells_selected.tolist()
    else:
        cells_of_type = adata.obs.index.values[adata.obs.level2 == celltype_select]
        cells_selected = cells_of_type[np.random.choice(np.arange(0, len(cells_of_type)),
                                                    size=min(9000, len(cells_of_type)), replace=False)]
        cells_selected_total = cells_selected_total + cells_selected.tolist()

print('Number of subsetted cells: ' + str(len(cells_selected_total)))

adata_sub = adata[cells_selected_total]
# create basic row and column attributes for the loom file:
f_loom_path_scenic = output_dir + "/loom/loom_hs_sub1.loom"
row_attrs = {"Gene": np.array(adata_sub.var_names)}
col_attrs = {
    "CellID": np.array(adata_sub.obs_names),
    #"tissue_condition": np.array(adata_sub.obs.tissue_condition),
    "nGene": np.array(np.sum(adata_sub.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(adata_sub.X.transpose(), axis=0)).flatten(),
}
lp.create(f_loom_path_scenic, adata_sub.X.transpose(), row_attrs, col_attrs)
