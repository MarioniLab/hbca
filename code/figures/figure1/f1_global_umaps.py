'''Here we make the global UMAP plots for figure1 (and supplementary)'''

#Use conda environment clustering

import scanpy as sc
import pandas as pd
import os

#rasterise plots
sc.set_figure_params(vector_friendly=True, dpi_save=500)


### Load Data

#adata
adata = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_2022-11-01.h5ad')
metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_all = metadata_all.set_index('cellID')
adata.obs = metadata_all[metadata_all.index.isin(adata.obs.index)]

### Colour Palletes

# celltype_colours = {'Lp': '#DDA0DD', 'Hs': '#EE3A8C', 'Bsl': '#FF6347',
#                     'Fb_1': '#AB6B4C', 'Fb_2': '#AB814C', 'Fb_1a': '#AB754C', 'Fb_2a': '#AB7B4C',
#                     'T-Lymphocyte': '#99A800',
#                     'BCell': '#A6A100', 'PlasmaCells': '#A8A500', 'DC': '#D0D111', 'Mo': '#D3D423',
#                     'EC_arterial': '#FFEC8B', 'EC_capillary': '#FFFC96', 'EC_venous': '#FFF69E', 'EC_lymphatic': '#CFBF7F',
#                     'Pericytes_1': '#F4A460', 'Pericytes_2': '#F4B264', 'Pericytes_3': '#F4C26E'}
# level0_colour_dictionary = {'Epithelial': '#DDA0DD', 'Stroma': '#804E1A', 'Immune': '#9FC5E8', 'Doublet': '#AAAAAA'} #'Stroma': #AB6B4C
level0_colour_dictionary = {'Epithelial': '#DEA2CB', 'Stroma': '#EEB449', 'Immune': '#95BE42', 'Doublet': '#AAAAAA'}

# level1_colour_dictionary = {'Luminal progenitor': '#DDA0DD', 'Luminal hormone sensing': '#EE3A8C', 'Basal': '#FF6347', # 'Epithelial Other': '#a67fbe',
#                             'Fibroblast': '#804E1A', 'Vascular mural': '#F89440', 'Endothelial': '#FFF08D', # VM: F4A460, 'EC vas': E0D461, 'Endothelial [lymphatic]': '#A59865',
#                              "Lymphoid": '#9FC5E8', "Myeloid": '#AAB256', 'Doublet': '#AAAAAA'}
# level1_colour_dictionary = {'Luminal progenitor': '#FFFF00', 'Luminal hormone sensing': '#006FA6', 'Basal': '#0000A6',
#                             'Fibroblast': '#8FB0FF', 'Vascular mural': '#1B4400', 'Endothelial': '#61615A',
#                              "Lymphoid": '#B903AA', "Myeloid": '#FFB500', 'Doublet': '#AAAAAA'}

# level1_colour_dictionary = {'Luminal progenitor': '#DDA0DD', 'Luminal hormone sensing': '#EE8298', 'Basal': '#E6554A',
#                             'Fibroblast': '#994E2E', 'Vascular mural': '#F89440', 'Endothelial': '#F1C232',  'Lymphatic Endothelial cell': '#4169E1',
#                             "T-cell": '#99A800', "Macrophage": '#64C6A6', 'B-cell': '#9FC5E8, 'Doublet': '#AAAAAA'}
level1_colour_dictionary = {'Luminal progenitor': '#DDA0DD', 'Luminal hormone sensing': '#EE8298', 'Basal': '#E6554A',
                            'Fibroblast': '#994E2E', 'Vascular mural': '#F89440', 'Endothelial': '#F1C232',  'Lymphatic Endothelial cell': '#4169E1',
                            "Lymphoid": '#99A800', "Myeloid": '#64C6A6', 'Doublet': '#AAAAAA'}

patientID_colour_set = ["#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                        "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                        "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                        "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                        "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                        "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                        "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                        "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                        "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                        "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                        "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                        "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                        "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]
sample_type_coarse_colour_dictionary = {'Organoid': '#fddd9d', 'Organoid LP': '#ff8d1f', 'Supernatant': '#76e1f8'}
sample_type_colour_dictionary = {'Organoid unsorted': '#fddd9d', 'Organoid LP sorted': '#ff8d1f',
                                 'Supernatant unsorted': '#76e1f8', 'Supernatant live-sorted': '#176db9'}
tissue_condition_colour_dictionary = {'Mammoplasty WT': '#f1ce63', 'Mastectomy BRCA1': '#4e79a7',
                                      'Mastectomy BRCA2': '#a0cbe8', 'Mastectomy WT': '#f28e2b',
                                      'Mastectomy unknown': '#bab0ac', 'Contralateral BRCA1':  '#ff9d9a'}
tissue_condition_colour_dictionary = {'AR': '#f1ce63', 'HR-BR1': '#4e79a7',
                                      'HR-BR2': '#a0cbe8', 'HR-unk': '#bab0ac', 'HR-other':  '#ff9d9a'}
parity_colour_dictionary = {'0': '#3158AE', '1': '#FBF349', '2': '#FBCC49', '3': '#FB8549', '4': '#FB4949', 'unknown': '#AAAAAA'}


# tissue_condition_colour_dictionary = {'Mammoplasty WT': '#49c918', 'Mastectomy WT': '#faf400', 'Mastectomy BRCA1': '#ff7529',
#                                       'Mastectomy BRCA2': '#c256ff', 'Mastectomy unknown': '#7e86b4', 'Contralateral BRCA1': '#ff4040'} #Old colour choices


### Make UMAPs

## Directories
fig1_dir = '/nfs/research/marioni/areed/projects/hbca/figures/figure1/'
supfig5_dir = '/nfs/research/marioni/areed/projects/hbca/figures/supfigure5/'

## Set up
adata.obs['tissue_condition'] = adata.obs.tissue_condition.astype('category').cat.reorder_categories(
    ['Mammoplasty WT', 'Mastectomy WT', 'Mastectomy BRCA1',
     'Mastectomy BRCA2', 'Mastectomy unknown', 'Contralateral BRCA1']
)


## Cell type (level) annotations
sc.settings.figdir = fig1_dir
os.makedirs(fig1_dir + 'umap/', exist_ok=True)

sc.pl.umap(adata,
           color=['level0'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level0_colour_dictionary,
           save='/fig1_level0.pdf')
sc.pl.umap(adata,
           color=['level0'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level0_colour_dictionary,
           save='/fig1_level0.png')
sc.pl.umap(adata,
           color=['level1'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level1_colour_dictionary,
           save='/fig1_level1.pdf')
sc.pl.umap(adata,
           color=['level1'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level1_colour_dictionary,
           save='/fig1_level1.png')
#also make one with the fixed level1 labels distinguishing lymphatic endothelial cells
adata.obs['level1_new'] = adata.obs.level1.copy()
adata.obs.level1_new[adata.obs.level2.isin(['LEC1', 'LEC2'])] = 'Lymphatic Endothelial cell'
sc.pl.umap(adata,
           color=['level1_new'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level1_colour_dictionary,
           save='/fig1_level1_new.pdf')
sc.pl.umap(adata,
           color=['level1_new'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level1_colour_dictionary,
           save='/fig1_level1_new.png')


sc.settings.figdir = supfig5_dir
os.makedirs(supfig5_dir + 'umap/', exist_ok=True)
## Experimental metadata plots
sc.pl.umap(adata,
           color=['patientID'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           # palette=patientID_colour_dictionary,
           palette=patientID_colour_set,
           save='/supfig5_patientID.pdf')
sc.pl.umap(adata,
           color=['patientID'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           # palette=patientID_colour_dictionary,
           palette=patientID_colour_set,
           save='/supfig5_patientID.png')
sc.pl.umap(adata,
           color=['before'], #sample_type
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=sample_type_colour_dictionary,
           save='/supfig5_sample_type.pdf')
sc.pl.umap(adata,
           color=['before'], #sample_type
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=sample_type_colour_dictionary,
           save='/supfig5_sample_type.png')
sc.pl.umap(adata,
           color=['sample_type_coarse'], #sample_type
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=sample_type_coarse_colour_dictionary,
           save='/supfig5_sample_type_coarse.pdf')
sc.pl.umap(adata,
           color=['sample_type_coarse'], #sample_type
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=sample_type_coarse_colour_dictionary,
           save='/supfig5_sample_type_coarse.png')
sc.pl.umap(adata,
           color=['tissue_condition'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=tissue_condition_colour_dictionary,
           save='/supfig5_tissue_condition.pdf')
sc.pl.umap(adata,
           color=['tissue_condition'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=tissue_condition_colour_dictionary,
           save='/supfig5_tissue_condition.png')
sc.pl.umap(adata,
           color=['patient_age'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           #palette=tissue_condition_colour_dictionary,
           save='/supfig5_patient_age.pdf')
sc.pl.umap(adata,
           color=['patient_age'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           #palette=tissue_condition_colour_dictionary,
           save='/supfig5_patient_age.png')
sc.pl.umap(adata,
           color=['parity'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=parity_colour_dictionary,
           save='/supfig5_parity.pdf')
sc.pl.umap(adata,
           color=['parity'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=parity_colour_dictionary,
           save='/supfig5_parity.png')










