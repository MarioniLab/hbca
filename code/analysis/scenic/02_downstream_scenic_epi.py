'''Analyse and visualise output from SCENIC'''

# import dependencies
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



### Read in scenic data


#load data
subset_number = '1'
subset_depth = 'sub' + subset_number
f_final_loom = './output/scenic_out/SCENIC_OUT_subset' + subset_number +'.loom'

# scenic output
lf = lp.connect(f_final_loom, mode='r', validate=False)
meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
exprMat = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i, r in pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene).iteritems():
    regulons[i] = list(r[r == 1].index.values)

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
    [pd.DataFrame(lf.ca.nGene, index=lf.ca.CellID),
     pd.DataFrame(lf.ca.nUMI, index=lf.ca.CellID),
     ],
    axis=1
)
cellAnnot.columns = ['nGene', 'nUMI']

lf.close()

adata = sc.read(f_final_loom, validate=False)

pdat = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
pdat = pdat.set_index("cellID")
pdat = pdat.loc[adata.obs_names, :] #pdat[pdat.cellID.isin(adata.obs_names)]
adata.obs = pdat #.reindex(adata.obs_names)



###Compute DimRed on AUC Matrix

#We can use the Regulon AUC Matrix to compute reduced dimensions. Batch effect remains here as this is on raw data, doesn't impact the questions we're trying ot address.

figdir = '/nfs/research/marioni/areed/projects/hbca/scenic/2022-04-05/scvi_new/epi/output/plots/'
sc.settings.figdir = figdir

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap(auc_mtx)
adata.obsm["scenic_UMAP"] = dr_umap

os.makedirs(figdir + 'umap/', exist_ok=True)
sc.settings.figdir = figdir + 'umap/'
sc.pl.embedding(adata,
           basis='scenic_UMAP',
           color=['level1', 'level2', 'patientID', 'before', 'processing_date', 'tissue_condition'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           #palette='Set1',
           save='_various_colours_' + subset_depth + '.png')
sc.pl.embedding(adata,
           basis='scenic_UMAP',
           color=['level1', 'level2', 'patientID', 'before', 'processing_date', 'tissue_condition'],
           legend_loc='on data',
           legend_fontsize=6,
           legend_fontoutline=2,
           #palette='Set1',
           save='_various_colours_' + subset_depth + 'ondata_.png')

#Build adata object with only the AUC_Matrix in for simplicity
aucdat = ad.AnnData(X=auc_mtx, obs=adata.obs)
aucdat.obsm["scenic_UMAP"] = dr_umap



### Compute RSS (regulon specifity score)

#Here we identify regulons that are specific to any of the cell types
rss_cellType = regulon_specificity_scores(auc_mtx, aucdat.obs.level2)
rss_tissueCond = regulon_specificity_scores(auc_mtx, aucdat.obs.tissue_condition)
print(rss_cellType)

cats = sorted(list(set(aucdat.obs.level2)))
cats2 = sorted(list(set(aucdat.obs.tissue_condition)))

##Karstens plot code
fig = plt.figure(figsize=(12, 4))
for c, num in zip(cats, range(1, len(cats) + 1)):
    x = rss_cellType.T[c]
    ax = fig.add_subplot(1, len(cats)+1, num)
    plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05, x.max() + (x.max() - x.min()) * 0.05)
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    # adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-', color='lightgrey'),
    #             precision=0.001) ###module not in container

fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium',
    'font.size': 5.0,
})
plt.tight_layout()
os.makedirs('/nfs/research/marioni/areed/projects/hbca/scenic/2022-04-05/scvi_new/epi/output/plots/regulons/', exist_ok=True)
plt.savefig(figdir + 'regulons/RSS_plots_level2_karsten_' + subset_depth + '.pdf')
plt.clf()

##tutorial plot code
# sns.set()
# sns.set(style='whitegrid', font_scale=0.8)
# fig, ((ax1, ax2, ax3, ax4, ax5, ax6, ax7)) = plt.subplots(1, 7, figsize=(8, 6), dpi=100)
# plot_rss(rss_cellType, 'LP1', ax=ax1)
# plot_rss(rss_cellType, 'LP2', ax=ax2)
# ax2.set_ylabel('')
# plot_rss(rss_cellType, 'LP3', ax=ax3)
# ax3.set_ylabel('')
# plot_rss(rss_cellType, 'LP4', ax=ax4)
# ax4.set_ylabel('')
# plot_rss(rss_cellType, 'LP proliferating', ax=ax5)
# ax5.set_ylabel('')
# plot_rss(rss_cellType, 'Other epithelial (1)', ax=ax6)
# ax6.set_ylabel('')
# plot_rss(rss_cellType, 'Other epithelial (2)', ax=ax7)
# plt.tight_layout()
# fig.savefig(figdir + 'regulons/RSS_plots_level2_tutorial_' + subset_depth + '.pdf')
# plt.clf()


###Heatmap of the top specific regulons

#Select top 10 regulons per cell type
topreg = []
for i, c in enumerate(cats):
    topreg.extend(list(rss_cellType.T[c].sort_values(ascending=False)[:10].index))
topreg = list(set(topreg))

#try for tissue_type
topreg2 = []
for i, c in enumerate(cats2):
    topreg2.extend(list(rss_tissueCond.T[c].sort_values(ascending=False)[:10].index))
topreg2 = list(set(topreg2))

#Convert to z-scors for plotting
auc_mtx_Z = pd.DataFrame(index=auc_mtx.index)
for col in list(auc_mtx.columns):
    auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

#make heatmap
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0 + idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

colors = sns.color_palette(["#EB675E", "#A23E36",
                            "#540F54", "#53407F",
                            "#EDABB9", "#EB5C79", "#A06A75", "#C00028",
                            "#DA80DA", "#815481", "#C040C0", "#E1AFE1", "#3F0034"])
colors2 = sns.color_palette(['#ff9d9a', '#f1ce63', '#4e79a7', '#a0cbe8', '#f28e2b', '#bab0ac'])
colorsd = dict(zip(cats, colors))
colorsd2 = dict(zip(cats2, colors2))
colormap = pd.DataFrame({'level2': [colorsd[x] for x in aucdat.obs.level2],
                         'tissue_condition': [colorsd2[x] for x in aucdat.obs.tissue_condition]})

sns.set()
sns.set(font_scale=0.8)
fig = palplot(colors, cats, size=1.0)
os.makedirs(figdir + 'heatmap/', exist_ok=True)
fig.savefig(figdir + 'heatmap/level2_colour_guide_' + subset_depth + '.pdf')
plt.clf()

sns.set()
sns.set(font_scale=0.8)
fig = palplot(colors2, cats2, size=1.0)
os.makedirs(figdir + 'heatmap/', exist_ok=True)
fig.savefig(figdir + 'heatmap/tissue_condition_colour_guide_' + subset_depth + '.pdf')
plt.clf()

sns.set(font_scale=1.2)
#auc_mtx_Z = auc_mtx_Z[~aucdat.obs.level2.isin(['Doublet'])]
g = sns.clustermap(auc_mtx_Z[topreg], annot=False, square=False, linecolor='gray',
                   yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
                   cmap="YlGnBu", figsize=(21, 16), rasterized=True)
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
g.savefig(figdir + 'heatmap/heatmap_' + subset_depth + '.png')
g.savefig(figdir + 'heatmap/heatmap_' + subset_depth + '.pdf')



#get ordered topreg
topreg = []
for i, c in enumerate(['LP1', 'LP2', 'LP3', 'LP4', 'LP5',
                       "HS1", "HS2", "HS3", "HS4",
                       "BSL1", "BSL2",
                       'DDC1', 'DDC2']):
    topreg.extend(list(rss_cellType.T[c].sort_values(ascending=False)[:10].index))
topreg = list(set(topreg))

#make heatmap without clustering rows/cols to see the columns identities more clearly.
if True:
    g_nocluster = sns.clustermap(auc_mtx_Z[topreg].reset_index(drop=True).T, annot=False, square=False, linecolor='gray',
                                 xticklabels=False, yticklabels=True, vmin=-2, vmax=4, col_colors=colormap,
                                 cmap="YlGnBu",
                                 figsize=(21, 16), col_cluster=False)
    g_nocluster.cax.set_visible(True)
    g_nocluster.ax_heatmap.set_ylabel('')
    g_nocluster.ax_heatmap.set_xlabel('')
    g_nocluster.savefig(figdir + 'heatmap/heatmap_nocluster' + subset_depth + '.png')

    g_nocluster.savefig(figdir + 'heatmap/heatmap_nocluster' + subset_depth + '.pdf')

g_nocluster_rast = sns.clustermap(auc_mtx_Z[topreg], annot=False, square=False, linecolor='gray',
                             yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
                             cmap="YlGnBu", figsize=(21, 16), row_cluster=False, rasterized=True)
g_nocluster_rast.cax.set_visible(True)
g_nocluster_rast.ax_heatmap.set_ylabel('')
g_nocluster_rast.ax_heatmap.set_xlabel('')
g_nocluster_rast.savefig(figdir + 'heatmap/heatmap_nocluster_rast' + subset_depth + '.pdf')

