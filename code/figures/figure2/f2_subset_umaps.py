'''Here we make the subsetted UMAP plots for figure2 (and supplementary)'''

#Use conda environment clustering

import scanpy as sc
import pandas as pd
import os

#rasterise plots
sc.set_figure_params(vector_friendly=True, dpi_save=500)

### Load Data

#adata
adata_epi = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2022-11-01.h5ad')
adata_str = sc.read('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/extra_preparation/output/data/HBCA_scVI_processing_date_sub_stroma_cleaned_newUMAP_2022-11-01.h5ad')
adata_imm = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_immune_ADDED_MP_2022-11-01.h5ad')

adata_lp = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LP_2022-11-01.h5ad')
adata_hs = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LHS_2022-11-01.h5ad')
adata_bsl = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_BSL_2022-11-01.h5ad')
adata_fb = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_FB_2022-11-01.h5ad')
adata_ec = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_endothelial_cleaned_2022-11-01.h5ad')
adata_vm = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_VM_2022-11-01.h5ad')

#metadata
metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
dblt_metadata_all = pd.read_csv('/nfs/research/marioni/kunz/HBCA/data/scrublet-scores/scrublet_scores_global.csv.gz')
dblt_metadata_all = dblt_metadata_all[['cellID', 'scrublet_score', 'scrublet_cluster_score']]
metadata_all = metadata_all.merge(dblt_metadata_all, how='left', on='cellID', sort=False)
metadata_all = metadata_all.set_index('cellID')

adata_epi.obs = metadata_all.loc[adata_epi.obs.index, :]
adata_str.obs = metadata_all.loc[adata_str.obs.index, :]
adata_imm.obs = metadata_all.loc[adata_imm.obs.index, :]

adata_lp.obs = metadata_all.loc[adata_lp.obs.index, :]
adata_hs.obs = metadata_all.loc[adata_hs.obs.index, :]
adata_bsl.obs = metadata_all.loc[adata_bsl.obs.index, :]
adata_fb.obs = metadata_all.loc[adata_fb.obs.index, :]
adata_ec.obs = metadata_all.loc[adata_ec.obs.index, :]
adata_vm.obs = metadata_all.loc[adata_vm.obs.index, :]

### Colour Palletes

# celltype_colours = {'Lp': '#DDA0DD', 'Hs': '#EE3A8C', 'Bsl': '#FF6347',
#                     'Fb_1': '#AB6B4C', 'Fb_2': '#AB814C', 'Fb_1a': '#AB754C', 'Fb_2a': '#AB7B4C',
#                     'T-Lymphocyte': '#99A800',
#                     'BCell': '#A6A100', 'PlasmaCells': '#A8A500', 'DC': '#D0D111', 'Mo': '#D3D423',
#                     'EC_arterial': '#FFEC8B', 'EC_capillary': '#FFFC96', 'EC_venous': '#FFF69E', 'EC_lymphatic': '#CFBF7F',
#                     'Pericytes_1': '#F4A460', 'Pericytes_2': '#F4B264', 'Pericytes_3': '#F4C26E'}
level0_colour_dictionary = {'Epithelial': '#DDA0DD', 'Stroma': '#804E1A', 'Immune': '#9FC5E8', 'Doublet': '#AAAAAA'} #'Stroma': #AB6B4C
level1_colour_dictionary = {'Luminal progenitor': '#DDA0DD', 'Luminal hormone sensing': '#EE3A8C', 'Basal': '#FF6347', # 'Epithelial Other': '#a67fbe',
                            'Fibroblast': '#804E1A', 'Vascular mural': '#F89440', 'Endothelial': '#FFF08D', # VM: F4A460, 'EC vas': E0D461, 'Endothelial [lymphatic]': '#A59865',
                             "Lymphoid": '#9FC5E8', "Myeloid": '#AAB256', 'Doublet': '#AAAAAA'}
# level2_colour_dictionary = {"LP1": '#DDA0DD', "LP2": '#F4C4F4', "LP3": '#C98EE8', "LP4": '#B578B5', #"LP5": '#C967C9', #LP5=Doublet now labelled so
#                             "LP proliferating": '#FC8FFC', "HS1": '#EE3A5D', "HS2": '#FF5039', "HS3": '#D06A87',
#                             "HS4": '#A53468', "BSL1": '#FF5039', "BSL2": '#931B06',
#                             "Other epithelial (1)": '#932AD3', "Other epithelial (2)": '#846598',
#                             "FB1": '#804E1A', "FB2": '#994E2E', "FB3": '#C49781', "FB4": '#601C00', "FB5": '#A36F56',
#                             "VM1": '#F89440', "VM2": '#E6A36A', "VM3": '#CB861E', "VM4": '#F7BF8F', "VM5": '#F4B24E',
#                             "EC venous": '#FFF08D', "EC capillary": '#F1C232', "EC arterial": '#CD9B01', #B9AF4A
#                             "EC angiogenic tip": '#F9E518', "LEC1": '#A59865', "LEC2": '#DDD998',
#                             'CD8T 1': '#99A800', 'CD8T 2': '#B7C24C', 'CD8T 3': '#CCD481', 'CD4T naive': '#E0E5B3',
#                             'IFNG+ T': '#86896B', 'NK1': '#64C6A6', 'NK2': '#D0EDE4', 'NK3': '#32635A',
#                             'ILC3': '#00BB6E', 'B cell': '#9FC5E8', 'Plasma cell': '#23D9F1', 'Macrophage': '#AAB256',
#                             'Doublet': '#AAAAAA'}
# level2_colour_dictionary = {
#     'LP1': "#FFFF00", 'LP2': "#1CE6FF", 'LP3': "#FF34FF", 'LP4': "#FF4A46", 'LP5': "#008941",
#     'HS1': "#006FA6", 'HS2': "#A30059", 'HS3': "#FFDBE5", 'HS4': "#7A4900",
#     'BSL1': "#0000A6", 'BSL2': "#63FFAC",
#     'DDC1': "#B79762", 'DDC2': "#004D43",
#     'FB1': "#8FB0FF", 'FB2': "#997D87", 'FB3': "#5A0007", 'FB4': "#809693", 'FB5': "#FEFFE6",
#     'VM1': "#1B4400", 'VM2': "#4FC601", 'VM3': "#3B5DFF", 'VM4': "#4A3B53", 'VM5': "#FF2F80",
#     'EC venous': "#61615A", 'EC_capillary': "#BA0900", 'EC_arterial': "#6B7900", 'EC_angiogenic_tip': "#00C2A0",
#     'LEC1': "#FFAA92", 'LEC2': "#FF90C9",
#     'CD8T 1': "#B903AA", 'CD8T 2': "#D16100", 'CD8T 3': "#DDEFFF", 'CD4T': "#000035", 'IFNG+ T': "#7B4F4B",
#     'NK1': "#A1C299", 'NK2': "#300018", 'NK3': "#0AA6D8",
#     'ILC3': "#013349",
#     'B cell': "#00846F", 'Plasma cell': "#372101",
#     'Macrophage': "#FFB500",
#     'Doublet': '#AAAAAA'}
level2_colour_dictionary = {
    'LP1': "#DA80DA", 'LP2': "#815481", 'LP3': "#C040C0", 'LP4': "#E1AFE1", 'LP5': "#3F0034",
    'HS1': "#EDABB9", 'HS2': "#EB5C79", 'HS3': "#A06A75", 'HS4': "#C00028",
    'BSL1': "#EB675E", 'BSL2': "#A23E36",
    'DDC1': "#540F54", 'DDC2': "#53407F",
    'FB1': "#DFA38A", 'FB2': "#8C3612", 'FB3': "#623623", 'FB4': "#916350", 'FB5': "#DAC3C3",
    'VM1': "#F8770B", 'VM2': "#E09E3A", 'VM3': "#CD7225", 'VM4': "#FFC990", 'VM5': "#AC5812",
    'EC venous': "#FEE083", 'EC capillary': "#897538", 'EC arterial': "#E7B419", 'EC angiogenic tip': "#BCA048",
    'LEC1': "#6F8BE2", 'LEC2': "#3053BC",
    'CD8T 1': "#6D9F58", 'CD8T 2': "#9EB766", 'CD8T 3': "#BDCB10", 'CD4T': "#3A6527", 'IFNG+ T': "#9EA743",
    'NK1': "#E2E8A7", 'NK2': "#5A6209", 'NK3': "#8FE36B",
    'ILC3': "#818A31",
    'B cell': "#9FC5E8", 'Plasma cell': "#23D9F1",
    'Macrophage': "#64C6A6",
    'Doublet': '#AAAAAA'}

#tissue_condition risk_status dictionary
tissue_condition_2_risk_status_dict = {
    'Mammoplasty WT': 'AR',
    'Mastectomy BRCA1': 'HR-BR1',
    'Mastectomy BRCA2': 'HR-BR2',
    'Mastectomy unknown': 'HR-unk',
    'Mastectomy WT': 'HR-unk',
    'Contralateral BRCA1': 'HR-other'
}

### Make UMAPs

## Directories
fig2_dir = '/nfs/research/marioni/areed/projects/hbca/figures/figure2/'
supfig6_dir = '/nfs/research/marioni/areed/projects/hbca/figures/supfigure6/'

## Set up

## Cell type (level2) annotations
sc.settings.figdir = fig2_dir
os.makedirs(fig2_dir + 'umap/', exist_ok=True)

sc.pl.umap(adata_epi,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_epi_level2.pdf')
sc.pl.umap(adata_epi,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_epi_level2.png')
sc.pl.umap(adata_str,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_str_level2.pdf')
sc.pl.umap(adata_str,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_str_level2.png')
sc.pl.umap(adata_imm,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_imm_level2.pdf')
sc.pl.umap(adata_imm,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_imm_level2.png')


## Experimental metadata plots
sc.settings.figdir = supfig6_dir
os.makedirs(supfig6_dir + 'umap/', exist_ok=True)

sc.pl.umap(adata_lp,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_lp_level2.pdf')
sc.pl.umap(adata_lp,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_lp_level2.png')
sc.pl.umap(adata_hs,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_hs_level2.pdf')
sc.pl.umap(adata_hs,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_hs_level2.png')
sc.pl.umap(adata_bsl,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_bsl_level2.pdf')
sc.pl.umap(adata_bsl,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_bsl_level2.png')
sc.pl.umap(adata_fb,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_fb_level2.pdf')
sc.pl.umap(adata_fb,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_fb_level2.png')
sc.pl.umap(adata_ec,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_ec_level2.pdf')
sc.pl.umap(adata_ec,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_ec_level2.png')
sc.pl.umap(adata_vm,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_vm_level2.pdf')
sc.pl.umap(adata_vm,
           color=['level2'],
           legend_loc='right margin',
           legend_fontsize=6,
           legend_fontoutline=2,
           palette=level2_colour_dictionary,
           save='/fig2_vm_level2.png')


