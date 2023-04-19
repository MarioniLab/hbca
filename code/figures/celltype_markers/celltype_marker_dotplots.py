'''Make Suplementary Figure 7(?) dotplots with celltype gene markers'''

if True: #for copying purposes
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import os
    from itertools import compress
    import colorcet as cc
    import seaborn as sns


#Load adata
adata_epi = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2022-11-01.h5ad')
adata_str = sc.read('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/extra_preparation/output/data/HBCA_scVI_processing_date_sub_stroma_cleaned_newUMAP_2022-11-01.h5ad')
adata_imm = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_immune_ADDED_MP_2022-11-01.h5ad')

#add metadata
if True: #for copying purposes
    metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
    dblt_metadata_all = pd.read_csv('/nfs/research/marioni/kunz/HBCA/data/scrublet-scores/scrublet_scores_global.csv.gz')
    dblt_metadata_all = dblt_metadata_all[['cellID', 'scrublet_score', 'scrublet_cluster_score']]
    metadata_all = metadata_all.merge(dblt_metadata_all, how='left', on='cellID', sort=False)
    metadata_all = metadata_all.set_index('cellID')

adata_epi.obs = metadata_all.loc[adata_epi.obs.index, :]
adata_str.obs = metadata_all.loc[adata_str.obs.index, :]
adata_imm.obs = metadata_all.loc[adata_imm.obs.index, :]

#premade lists
epi_celltype_order = ["LP1", "LP2", "LP3", "LP4", "LP5",
                      "HS1", "HS2", "HS3", "HS4",
                      "BSL1", "BSL2", 'DDC1', 'DDC2']
str_celltype_order = ["FB1", "FB2", "FB3", "FB4", "FB5",
                      "VM1", "VM2", "VM3", "VM4", "VM5",
                      "EC venous", "EC capillary", "EC arterial",
                      "EC angiogenic tip", "LEC1", "LEC2"]
imm_celltype_order = ['CD4T', 'CD8T 1', 'CD8T 2', 'CD8T 3',
                      'IFNG+ T', 'NK1', 'NK2', 'NK3',
                      'ILC3', 'B cell', 'Plasma cell', 'Macrophage']


# plot gene markers organised by leiden clusters for cluster annotation.
fig_dir = '/nfs/research/marioni/areed/projects/hbca/figures/src/marker_plots/output/'
sc.settings.figdir = fig_dir
os.makedirs(fig_dir + '/dotplot_1/', exist_ok=True)

#Epi
if True:
    gene_markers_epi = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/misc/markers_paper/epithelia_markers.csv')
    gene_marker_dict_epi = {}
    for key in gene_markers_epi['cell_type'].unique():
        gene_marker_dict_epi[key] = gene_markers_epi['gene_name'][gene_markers_epi['cell_type'] == key]

sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ], var_names=gene_marker_dict_epi,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=epi_celltype_order,
              save='1/epi_markergenes_dotplot.png')
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ], var_names=gene_marker_dict_epi,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=epi_celltype_order,
              save='1/epi_markergenes_dotplot.pdf')
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ], var_names=gene_marker_dict_epi,
              groupby='level2', dendrogram=False, swap_axes=True,
              standard_scale='var', categories_order=epi_celltype_order,
              save='1/epi_markergenes_dotplot_rotated.pdf')
epi_celltype_order.reverse() #reverse order
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ], var_names=gene_marker_dict_epi,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=epi_celltype_order,
              save='1/epi_markergenes_dotplot_reverse.pdf')
epi_celltype_order.reverse() #fix order

#Str
if True:
    gene_markers_str = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/misc/markers_paper/stroma_markers.csv')
    gene_marker_dict_str = {}
    for key in gene_markers_str['cell_type'].unique():
        gene_marker_dict_str[key] = gene_markers_str['gene_name'][gene_markers_str['cell_type'] == key]

sc.pl.dotplot(adata_str[adata_str.obs.level2.isin(str_celltype_order), ], var_names=gene_marker_dict_str,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=str_celltype_order,
              save='1/str_markergenes_dotplot.png')
sc.pl.dotplot(adata_str[adata_str.obs.level2.isin(str_celltype_order), ], var_names=gene_marker_dict_str,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=str_celltype_order,
              save='1/str_markergenes_dotplot.pdf')
sc.pl.dotplot(adata_str[adata_str.obs.level2.isin(str_celltype_order), ], var_names=gene_marker_dict_str,
              groupby='level2', dendrogram=False, swap_axes=True,
              standard_scale='var', categories_order=str_celltype_order,
              save='1/str_markergenes_dotplot_rotated.pdf')
str_celltype_order.reverse() #switch
sc.pl.dotplot(adata_str[adata_str.obs.level2.isin(str_celltype_order), ], var_names=gene_marker_dict_str,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=str_celltype_order,
              save='1/str_markergenes_dotplot_reverse.pdf')
str_celltype_order.reverse() #fix

#Imm
if True:
    gene_markers_imm = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/misc/markers_paper/immune_markers.csv')
    gene_marker_dict_imm = {}
    for key in gene_markers_imm['cell_type'].unique():
        gene_marker_dict_imm[key] = gene_markers_imm['gene_name'][gene_markers_imm['cell_type'] == key]

sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(imm_celltype_order), ], var_names=gene_marker_dict_imm,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=imm_celltype_order,
              save='1/imm_markergenes_dotplot.png')
sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(imm_celltype_order), ], var_names=gene_marker_dict_imm,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=imm_celltype_order,
              save='1/imm_markergenes_dotplot.pdf')
sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(imm_celltype_order), ], var_names=gene_marker_dict_imm,
              groupby='level2', dendrogram=False, swap_axes=True,
              standard_scale='var', categories_order=imm_celltype_order,
              save='1/imm_markergenes_dotplot_rotated.pdf')
imm_celltype_order.reverse() #switch
sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(imm_celltype_order), ], var_names=gene_marker_dict_imm,
              groupby='level2', dendrogram=False,
              standard_scale='var', categories_order=imm_celltype_order,
              save='1/imm_markergenes_dotplot_reverse.pdf')
imm_celltype_order.reverse() #fix



## Plot some specific immune markers (exhaustion etc.)

cytotoxic_celltypes = ['CD8T 1', 'CD8T 2', 'CD8T 3', 'IFNG+ T', 'NK1', 'NK2', 'NK3']
#tissue_types = ['Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2']
tc_types = ['AR', 'HR-BR1', 'HR-BR2']
tc_map = {'Mammoplasty WT': 'AR', 'Mastectomy BRCA1': 'HR-BR1', 'Mastectomy BRCA2': 'HR-BR2', 'Mastectomy WT': 'HR-unk',
          'Mastectomy Unknown': 'HR-unk'}
adata_imm.obs['tc_new'] = adata_imm.obs.tissue_condition.map(tc_map)
exhaustion_marker_dict = {#'Control': ['PTPRC'],
                          'ICRs': ['PDCD1', 'CTLA4', 'LAG3', 'TIGIT', 'HAVCR2']} #HAVCR2 = TIM3
exhaustion_marker_dict2 = ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'NFATC2', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'FOXO1', 'TOX', 'CD28']
adata_imm.obs['tc_new_level2'] = adata_imm.obs['level2'] + ' - ' + adata_imm.obs['tc_new']
sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(cytotoxic_celltypes) & adata_imm.obs.tc_new.isin(tc_types), ],
              var_names=exhaustion_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=cytotoxic_celltypes,
              save='1/imm_exhaustionmarkers_dotplot.png')
sc.pl.dotplot(adata_imm[adata_imm.obs.level2.isin(cytotoxic_celltypes) & adata_imm.obs.tc_new.isin(tc_types), ],
              var_names=exhaustion_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=cytotoxic_celltypes,
              save='1/imm_exhaustionmarkers_dotplot.pdf')

#stacked violin plots
os.makedirs(fig_dir + '/stacked_violin_1/', exist_ok=True)
sc.pl.stacked_violin(adata_imm[adata_imm.obs.level2.isin(cytotoxic_celltypes) & adata_imm.obs.tc_new.isin(tc_types), ],
              var_names=exhaustion_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=cytotoxic_celltypes,
              save='1/imm_exhaustionmarkers_violins.pdf')


#immune checkpoint ligand expression
interesting_celltypes = ["LP1", "LP2", "LP3", "LP4", "LP5",
                         "HS1", "HS2", "HS3", "HS4",
                         "BSL1", "BSL2", 'DDC1', 'DDC2']
#tissue_types = ['Mammoplasty WT', 'Mastectomy BRCA1', 'Mastectomy BRCA2']
tc_types = ['AR', 'HR-BR1', 'HR-BR2']
tc_map = {'Mammoplasty WT': 'AR', 'Mastectomy BRCA1': 'HR-BR1', 'Mastectomy BRCA2': 'HR-BR2', 'Mastectomy WT': 'HR-unk',
          'Mastectomy Unknown': 'HR-unk'}
adata_epi.obs['tc_new'] = adata_epi.obs.tissue_condition.map(tc_map)
checkpoint_ligands_marker_dict = {'ICRs': ['CD274', 'PDCD1LG2', 'CD80', 'CD86', 'PVR', 'LGALS9']} #, 'HLA-A', 'HLA-B', 'HLA-C']} #PVR=CD155
adata_epi.obs['tc_new_level1'] = adata_epi.obs['level1'] + ' - ' + adata_epi.obs['tc_new']
adata_epi.obs['tc_new_level2'] = adata_epi.obs['level2'] + ' - ' + adata_epi.obs['tc_new']

sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(interesting_celltypes) & adata_epi.obs.tc_new.isin(tc_types), ],
              var_names=checkpoint_ligands_marker_dict,
              groupby='tc_new_level1', dendrogram=False,
              standard_scale='var', #categories_order=interesting_celltypes,
              save='1/epi_checkpoint_ligands_level1_dotplot.png')
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(interesting_celltypes) & adata_epi.obs.tc_new.isin(tc_types), ],
              var_names=checkpoint_ligands_marker_dict,
              groupby='tc_new_level1', dendrogram=False,
              standard_scale='var', #categories_order=interesting_celltypes,
              save='1/epi_checkpoint_ligands_level1_dotplot.pdf')

sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(interesting_celltypes) & adata_epi.obs.tc_new.isin(tc_types), ],
              var_names=checkpoint_ligands_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=interesting_celltypes,
              save='1/epi_checkpoint_ligands_level2_dotplot.png')
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(interesting_celltypes) & adata_epi.obs.tc_new.isin(tc_types), ],
              var_names=checkpoint_ligands_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=interesting_celltypes,
              save='1/epi_checkpoint_ligands_level2_dotplot.pdf')

#stacked violin plots
os.makedirs(fig_dir + '/stacked_violin_1/', exist_ok=True)
sc.pl.stacked_violin(adata_epi[adata_epi.obs.level2.isin(interesting_celltypes) & adata_epi.obs.tc_new.isin(tc_types), ],
              var_names=checkpoint_ligands_marker_dict,
              groupby='tc_new_level2', dendrogram=False,
              standard_scale='var', #categories_order=interesting_celltypes,
              save='1/epi_checkpoint_ligands_violins.pdf')



#LP4/BSL2/HS2 specific Marker plotting

epi_dict = {'Stress': ['JUN', 'FOS', 'HSPA1A'], 'Other': ['RECQL', 'FBXW7', 'MCL1', 'GADD45B', 'GADD45G']}

sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ],
              var_names=epi_dict,
              groupby='level2', dendrogram=False,
              standard_scale='var',
              categories_order=epi_celltype_order,
              save='1/epi_LP4HS2BSL2_dotplot.png')
sc.pl.dotplot(adata_epi[adata_epi.obs.level2.isin(epi_celltype_order), ],
              var_names=epi_dict,
              groupby='level2', dendrogram=False,
              standard_scale='var',
              categories_order=epi_celltype_order,
              save='1/epi_LP4HS2BSL2_dotplot.pdf')

