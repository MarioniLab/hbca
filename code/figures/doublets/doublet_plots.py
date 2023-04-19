'''Here we make the doublet detection plots.'''

#Use conda environment clustering

if True:
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import anndata as ad
    import seaborn as sns

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
if True:
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

#fixed plotting error #not sure if neccessary
adata_epi.uns['log1p']["base"] = None
adata_str.uns['log1p']["base"] = None
adata_imm.uns['log1p']["base"] = None
adata_lp.uns['log1p']["base"] = None
adata_hs.uns['log1p']["base"] = None
adata_bsl.uns['log1p']["base"] = None
adata_fb.uns['log1p']["base"] = None
adata_ec.uns['log1p']["base"] = None
adata_vm.uns['log1p']["base"] = None

adata_dict = {'epithelial': adata_epi, 'stroma': adata_str, 'immune': adata_imm,
              'lp': adata_lp, 'hs': adata_hs, 'bsl': adata_bsl,
              'fb': adata_fb, 'ec': adata_ec, 'vm': adata_vm}

### Mapping dictionaries

level2_colour_dictionary = {
    'LP1': "#DA80DA", 'LP2': "#815481", 'LP3': "#C040C0", 'LP4': "#E1AFE1", 'LP5': "#3F0034",
    'HS1': "#EDABB9", 'HS2': "#EB5C79", 'HS3': "#A06A75", 'HS4': "#C00028",
    'BSL1': "#EB675E", 'BSL2': "#A23E36",
    'DDC1': "#540F54", 'DDC2': "#53407F",
    'FB1': "#DFA38A", 'FB2': "#8C3612", 'FB3': "#623623", 'FB4': "#916350", 'FB5': "#DAC3C3",
    'VM1': "#F8770B", 'VM2': "#E09E3A", 'VM3': "#CD7225", 'VM4': "#FFC990", 'VM5': "#AC5812",
    #'EC venous': "#FEE083", 'EC capillary': "#897538", 'EC arterial': "#E7B419", 'EC angiogenic tip': "#BCA048",
    'ECV': "#FEE083", 'ECC': "#897538", 'ECA': "#E7B419", 'ECAT': "#BCA048",
    'LEC1': "#6F8BE2", 'LEC2': "#3053BC",
    'CD8T 1': "#6D9F58", 'CD8T 2': "#9EB766", 'CD8T 3': "#BDCB10", 'CD4T': "#3A6527", 'IFNG+ T': "#9EA743",
    'NK1': "#E2E8A7", 'NK2': "#5A6209", 'NK3': "#8FE36B",
    'ILC3': "#818A31",
    'B cell': "#9FC5E8", 'Plasma cell': "#23D9F1",
    'Macrophage': "#64C6A6",
    'Doublet': '#AAAAAA'}

level2_update_dictionary = {'EC venous': "ECV", 'EC capillary': "ECC", 'EC arterial': "ECA", 'EC angiogenic tip': "ECAT"}

celltype_order = np.array(list(level2_colour_dictionary.keys()))


for cell_group in adata_dict.keys():
    print(cell_group)
    figdir = '/nfs/research/marioni/areed/projects/hbca/figures/src/doublets/output/' + cell_group + '/'
    sc.settings.figdir = figdir

    #load data
    adata_group = adata_dict[cell_group]
    adata_group.var_names_make_unique()
    adata_group.obs['cellID'] = adata_group.obs.index

    adata_group.obs.level2 = adata_group.obs.level2.replace(level2_update_dictionary)

    #celltypes order for plots
    binary_celltype_choice = np.array(pd.DataFrame(celltype_order).isin(adata_group.obs.level2.unique())[0])
    celltype_order_use = celltype_order[binary_celltype_choice]


    fig_dir = '/nfs/research/marioni/areed/projects/hbca/figures/src/doublets/output/' + cell_group + '/celltype_proportions/'
    os.makedirs(fig_dir, exist_ok=True)

    # umi counts
    fig1, ax1 = plt.subplots()
    ax1 = sns.boxplot(data=adata_group.obs, x='level2', y='n_counts',
                      order=celltype_order_use,
                      palette=level2_colour_dictionary)
    plt.xticks(rotation=90)
    plt.tight_layout()
    fig1.savefig(fig_dir + 'boxplot_n_counts.pdf')
    fig1.get_axes()[0].set_yscale('log')
    fig1.savefig(fig_dir + 'boxplot_n_counts_log.pdf')

    # genes detected
    fig2, ax2 = plt.subplots()
    ax2 = sns.boxplot(data=adata_group.obs, x='level2', y='n_genes',
                      order=celltype_order_use,
                      palette=level2_colour_dictionary)
    plt.xticks(rotation=90)
    plt.tight_layout()
    fig2.savefig(fig_dir + 'boxplot_n_genes.pdf')
    fig2.get_axes()[0].set_yscale('log')
    fig2.savefig(fig_dir + 'boxplot_n_genes_log.pdf')

    # mito percent
    fig3, ax3 = plt.subplots()
    ax3 = sns.boxplot(data=adata_group.obs, x='level2', y='percent_mito',
                      order=celltype_order_use,
                      palette=level2_colour_dictionary)
    plt.xticks(rotation=90)
    plt.tight_layout()
    fig3.savefig(fig_dir + 'boxplot_percent_mito.pdf')
    fig3.get_axes()[0].set_yscale('log')
    fig3.savefig(fig_dir + 'boxplot_percent_mito_log.pdf')

    # Doublet scores
if True:
    fig4, ax4 = plt.subplots()
    ax4 = sns.boxplot(data=adata_group.obs, x='level2', y='scrublet_score',
                      order=celltype_order_use,
                      palette=level2_colour_dictionary)
    plt.xticks(rotation=90)
    plt.tight_layout()
    fig4.savefig(fig_dir + 'boxplot_scrublet_score.pdf')
    fig4.get_axes()[0].set_yscale('log')
    fig4.savefig(fig_dir + 'boxplot_scrublet_score_log.pdf')

