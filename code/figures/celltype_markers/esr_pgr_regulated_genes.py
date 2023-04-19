'''Make Suplementary dotplots with estrogen and progesterone regulated genes'''

if True: #for copying purposes
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import os
    from itertools import compress
    import colorcet as cc
    import seaborn as sns


#Load data
adata_hs = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LHS_2022-11-01.h5ad')

#add metadata
if True: #for copying purposes
    metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
    dblt_metadata_all = pd.read_csv('/nfs/research/marioni/kunz/HBCA/data/scrublet-scores/scrublet_scores_global.csv.gz')
    dblt_metadata_all = dblt_metadata_all[['cellID', 'scrublet_score', 'scrublet_cluster_score']]
    metadata_all = metadata_all.merge(dblt_metadata_all, how='left', on='cellID', sort=False)
    metadata_all = metadata_all.set_index('cellID')

adata_hs.obs = metadata_all.loc[adata_hs.obs.index, :]

hs_celltype_order = ["HS1", "HS2", "HS3", "HS4"]

# plot gene markers
fig_dir = '/nfs/research/marioni/areed/projects/hbca/figures/src/marker_plots/output/'
sc.settings.figdir = fig_dir
os.makedirs(fig_dir + '/dotplot_esr_pgr/', exist_ok=True)

esr_pgr_dict = {'Estrogen-regulated': ['ESR1', 'SERPINA1', 'CITED1', 'PDZK1'],
                          'Progesterone-regulated': ['PGR', 'BTG1', 'KCNK1', 'AKAP13', 'WNT4']}

sc.pl.dotplot(adata_hs[adata_hs.obs.level2.isin(hs_celltype_order)],
              var_names=esr_pgr_dict,
              groupby='level2', dendrogram=False,
              standard_scale='var',
              save='esr_pgr/hs_dotplot1.png')



