'''Analyse and visualise output from SCENIC'''

# import dependencies
if True:
    import os
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import loompy as lp
    import anndata as ad
    import umap
    import json
    import base64
    import zlib
    import re
    from pyscenic.binarization import binarize
    from pyscenic.plotting import plot_binarization
    from pyscenic.plotting import plot_rss
    from pyscenic.export import add_scenic_metadata
    from pyscenic.cli.utils import load_signatures
    from pyscenic.rss import regulon_specificity_scores
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
#from adjustText import adjust_text ###module not in container

sc.set_figure_params(vector_friendly=True, dpi_save=500)


### Read in scenic data


### Load Data

#adata
adata_epi = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2022-11-01.h5ad')
metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_all = metadata_all.set_index('cellID')
adata_epi.obs = metadata_all.loc[adata_epi.obs.index, :]

#scenic
subset_number = '1'
subset_depth = 'sub' + subset_number
f_final_loom = './output/scenic_out/SCENIC_OUT_subset' + subset_number +'.loom'

lf = lp.connect(f_final_loom, mode='r', validate=False)
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

#normalize
auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
for col in list(auc_mtx.columns):
    auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

#Get regulons of interest into adata.obs
if True:
    interesting_regulons = ['BRCA1(+)', 'MYBL1(+)', #'BRCA2(+)', 'TP53(+)', #seem to be missing?   ##LP5
                            'EZH2(+)', 'E2F1(+)', 'E2F2(+)', 'E2F3(+)',    ##LP5
                            'E2F4(+)', 'E2F5(+)', 'E2F6(+)', 'E2F7(+)', 'E2F8(+)',    #LP5
                            'CEBPD(+)', 'STAT1(+)', 'ETV7(+)',  ##LP1
                            'ETS1(+)',  'ELK3(+)',#'GLIS3(+)', 'FOXJ3(+)', 'PRKAA1(+)', 'RREB1(+)', 'SMAD3(+)', ##LP2
                            'TEAD1(+)', 'CEBPG(+)',  ##LP3
                            'DLX5(+)', 'SOX8(+)', #'OLIG1(+)',#'ISL1(+)', 'ETV1(+)', 'SOX14(+)', ##LP4
                            'STAT2(+)', 'NFKB1(+)', 'TP63(+)', ##BSL1
                            'SNAI2(+)', 'IRF4(+)', ##BSL2
                            'POU2F3(+)', 'SOX9(+)', 'FOXA1(+)', ##HS1
                            'ESR1(+)', 'HES1(+)', #'MYCN(+)', #'MYB(+)', ##HS2
                            'LEF1(+)', 'RORC(+)', 'EHF(+)', 'FOXC1(+)',   ##HS3
                            'ETV5(+)', 'FOXI1(+)', 'HAND2(+)', 'CUX1(+)', 'RB1(+)', #'TEAP2B(+)',  ##HS4
                            'ZEB1(+)',
                            'TFEB(+)', #'PRA(+)', 'PRB(+)'
                            ]
    adata_epi.obs[interesting_regulons] = 'NA'#pd.NA
    adata_epi.obs.loc[adata_epi.obs_names.isin(auc_mtx_Z.index), interesting_regulons] = auc_mtx_Z[interesting_regulons]

#make umap plots for the regulons of interest
fig_dir = '/nfs/research/marioni/areed/projects/hbca/scenic/2022-04-05/scvi_new/epi/output/plots/'
sc.settings.figdir = fig_dir
os.makedirs(fig_dir + 'umap_rna/', exist_ok=True)

adata_epi_sub = adata_epi[adata_epi.obs_names.isin(auc_mtx_Z.index), ].copy()
for regulon in interesting_regulons:
    sc.pl.umap(adata_epi_sub,
               color=[regulon],
               legend_loc='right margin',
               legend_fontsize=6,
               legend_fontoutline=2,
               vmin=-2, vmax=4,
               size=500000 / adata_epi_sub.obs.shape[0],
               save='_rna/epi_umap_' + regulon + '.pdf')


