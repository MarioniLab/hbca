'''Label cells based off scvi embedding clustering.'''

import pandas as pd
import scanpy as sc
import os

# set up directory for saving figures
fig_dir = '/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/labelling/output/level2'


#load data
adata_all = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_2023-06-19.h5ad')
metadata_all = adata_all.obs

adata_epi = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_2023-06-19.h5ad')
adata_str = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_2023-06-19.h5ad')
adata_imm = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_immune_cleaned_2023-06-19.h5ad')
adata_lp = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LP_2023-06-19.h5ad')
adata_hs = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_LHS_2023-06-19.h5ad')
adata_bsl = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_epithelial_cleaned_sub_BSL_2023-06-19.h5ad')
adata_ec = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_endothelial_2023-06-19.h5ad')
adata_vm = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_VM_2023-06-19.h5ad')
adata_fb = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_stroma_cleaned_sub_FB_2023-06-19.h5ad')
adata_BLym = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_B-cells_plasma_2023-06-19.h5ad')
adata_mye = sc.read('/nfs/research/marioni/kunz/HBCA/data/HBCA_scVI_processing_date_sub_immune_cleaned_sub_myeloid_2023-06-19.h5ad')


#dictionaries for mapping levels

#Update names to the new nomenclature
lvl1_dict = {'Luminal progenitor': 'Luminal adaptive secretory precurser',
             'Luminal hormone sensing': 'Luminal hormone sensing',
             'Basal': 'Basal-myoepithelial',
             'Endothelial [vascular]': 'Vascular endothelial',
             'Endothelial [lymphatic]': 'Lymphatic endothelial',
             'Vascular mural': 'Perivascular',
             'Fibroblast': 'Fibroblast',
             'Immune': 'Immune',
             'Doublet': 'Doublet'}

#lvl2
lp_dict_lvl2 = {'0': 'LASP1',
                '1': 'LASP2',
                '2': 'LASP3',
                '3,0': 'LASP4',
                '3,1': 'LASP5',
                '4': 'Doublet',
                '5': 'DDC2'}
hs_dict_lvl2 = {'0': 'LHS1',
                '1': 'LHS2',
                '2,0': 'LHS3',
                '2,1': 'Doublet',
                '3,0': 'stripped_nuclei',
                '3,1': 'stripped_nuclei',
                '4': 'DDC1'}
bsl_dict_lvl2 = {'0': 'BMYO1',
                 '1': 'BMYO1',
                 '2': 'BMYO1',
                 '3': 'BMYO2',
                 '4': 'Doublet'}
fb_dict_lvl2 = {'0': 'FB1',
                '1': 'FB2',
                '2': 'FB3',
                '3': 'FB4',
                '4': 'Doublet'}
ec_dict_lvl2 = {'0': 'VEV',
                '1': 'VEAT',
                '2,0': 'VEV',
                '2,1': 'VEC',
                '2,2': 'VEC',
                '3': 'VEA',
                '4,0': 'LE1',
                '4,1': 'LE2',
                '5': 'VEA'}
vm_dict_lvl2 = {'0': 'PV1',
                '1': 'PV3',
                '2': 'PV4',
                '3': 'PV2',
                '4': 'PV5',
                '5': 'PV1',
                '6': 'Doublet'}
imm_dict_lvl2 = {'0': 'CD8_Trm',
                 '1,0': 'CD4_naive',
                 '1,1': 'CD4_Th',
                 '2': 'Bcell',
                 '3': 'Myeloid',
                 '4': 'CD8_Tc1',
                 '5,0': 'CD8_Tem',
                 '5,1': 'CD8_Trm',
                 '5,2,0': 'NKT',
                 '5,2,1': 'NK',
                 '6': 'stripped_nuclei',
                 '7': 'Plasma',
                 '8': 'CD8_Tc1',
                 '9': 'ILC'}
Blym_dict_lvl2 = {'0': 'B_mem_switched',
                  '1': 'B_mem_switched',
                  '2': 'B_mem_unswitched',
                  '3': 'B_naive',
                  '4': 'B_naive',
                  '5': 'Plasma_cell',
                  '6': 'stripped_nuclei',
                  '7': 'Plasma_cell'}
mye_dict_lvl2 = {'0': 'Macro',
                 '1': 'Macro',
                 '2': 'Macro',
                 '3': 'Macro',
                 '4': 'DC',
                 '5': 'Macro-lipo'}

### Map cell types level2

#LP mapping
sc.tl.leiden(adata_lp, resolution=0.04, restrict_to=('leiden_0_5', ['3']), random_state=0,
             key_added='leiden_0_5_sub_3')
adata_lp.obs['level2'] = adata_lp.obs['leiden_0_5_sub_3'].map(lp_dict_lvl2)

#HS mapping
sc.tl.leiden(adata_hs, resolution=0.1, restrict_to=('leiden_0_2', ['3']), random_state=0,
                 key_added='leiden_0_2_sub_3')
sc.tl.leiden(adata_hs, resolution=0.1, restrict_to=('leiden_0_2_sub_3', ['2']), random_state=0,
             key_added='leiden_0_2_sub_3_sub_2')
adata_hs.obs['level2'] = adata_hs.obs['leiden_0_2_sub_3_sub_2'].map(hs_dict_lvl2)

#BSL mapping
adata_bsl.obs['level2'] = adata_bsl.obs['leiden_0_3'].map(bsl_dict_lvl2)

#FB mapping
adata_fb.obs['level2'] = adata_fb.obs['leiden_0_2'].map(fb_dict_lvl2)

#EC mapping
sc.tl.leiden(adata_ec, resolution=0.2, restrict_to=('leiden_0_4', ['2']), random_state=0,
                 key_added='leiden_0_4_sub_2')
sc.tl.leiden(adata_ec, resolution=0.1, restrict_to=('leiden_0_4_sub_2', ['4']), random_state=0,
             key_added='leiden_0_4_sub_2_sub_4')
adata_ec.obs['level2'] = adata_ec.obs['leiden_0_4_sub_2_sub_4'].map(ec_dict_lvl2)

#VM mapping
adata_vm.obs['level2'] = adata_vm.obs['leiden_0_4'].map(vm_dict_lvl2)

#Imm mapping
sc.tl.leiden(adata_imm, resolution=0.2, restrict_to=('leiden_0_5', ['5']), random_state=0,
                 key_added='leiden_0_5_sub_5')
sc.tl.leiden(adata_imm, resolution=0.2, restrict_to=('leiden_0_5_sub_5', ['5,2']), random_state=0,
             key_added='leiden_0_5_sub_5_sub_5,2')
sc.tl.leiden(adata_imm, resolution=0.2, restrict_to=('leiden_0_5_sub_5_sub_5,2', ['1']), random_state=0,
             key_added='leiden_0_5_sub_5_sub_5,2_sub_1')
adata_imm.obs['level2'] = adata_imm.obs['leiden_0_5_sub_5_sub_5,2_sub_1'].map(imm_dict_lvl2)

#immune_subsets and map these to adata_imm
adata_BLym.obs['level2'] = adata_BLym.obs['leiden_0_5'].map(Blym_dict_lvl2)
adata_mye.obs['level2'] = adata_mye.obs['leiden_0_4'].map(mye_dict_lvl2)

adata_imm.obs.level2[adata_imm.obs.level2.isin(['Plasma', 'Bcell']).values] = adata_BLym.obs.level2.copy()
adata_imm.obs.level2[adata_imm.obs.level2.isin(['Myeloid']).values] = adata_mye.obs.level2.copy()

### Add labels to metadata and save
dfs_lvl1 = [adata_epi.obs.level1.reset_index(level=0), adata_str.obs.level1.reset_index(level=0),
           adata_imm.obs.level1.reset_index(level=0)]
dfs_lvl2 = [adata_lp.obs.level2.reset_index(level=0), adata_hs.obs.level2.reset_index(level=0),
           adata_bsl.obs.level2.reset_index(level=0),
           adata_fb.obs.level2.reset_index(level=0), adata_ec.obs.level2.reset_index(level=0),
           adata_vm.obs.level2.reset_index(level=0),
           adata_imm.obs.level2.reset_index(level=0)]
merged_level1 = pd.concat(dfs_lvl1).groupby('cellID', sort=False, as_index=False).first()
merged_level1 = merged_level1.set_index('cellID')
merged_level2 = pd.concat(dfs_lvl2).groupby('cellID', sort=False, as_index=False).first()
merged_level2 = merged_level2.set_index('cellID')
metadata_all = metadata_all.drop(labels=['level1', 'leiden_0_05', 'leiden_0_1', 'leiden_0_2', 'leiden_0_3',
                                         'leiden_0_4', 'leiden_0_5', 'leiden_0_6', 'leiden_0_7', 'leiden_0_8',
                                         'leiden_0_9', 'leiden_1_0'], axis=1)
metadata_all = metadata_all.join(merged_level1, how='left', sort=False)
metadata_all = metadata_all.join(merged_level2, how='left', sort=False)

#update level1 nomenclature
metadata_all.level1 = metadata_all.level1.map(lvl1_dict)

#make initial cluster levels which we will overwrite with doublet labels etc
metadata_all['level0_global'] = metadata_all.level0
metadata_all['level1_global'] = metadata_all.level1

#map level2 to level1 and level0 for consistency
level2_to_level1_dict = {'LASP1': "Luminal adaptive secretory precurser", 'LASP2': "Luminal adaptive secretory precurser",
                         'LASP3': "Luminal adaptive secretory precurser", 'LASP4': "Luminal adaptive secretory precurser",
                         'LASP5': "Luminal adaptive secretory precurser",
                         'LHS1': "Luminal hormone sensing", 'LHS2': "Luminal hormone sensing",
                         'LHS3': "Luminal hormone sensing",
                         'BMYO1': "Basal-myoepithelial", 'BMYO2': "Basal-myoepithelial",
                         'DDC1': "DDC1", 'DDC2': "DDC2",
                         'FB1': "Fibroblast", 'FB2': "Fibroblast", 'FB3': "Fibroblast", 'FB4': "Fibroblast",
                         'FB5': "Fibroblast",
                         'PV1': "Perivascular", 'PV2': "Perivascular", 'PV3': "Perivascular", 'PV4': "Perivascular",
                         'PV5': "Perivascular",
                         'VEV': "Vascular endothelial", 'VEC': "Vascular endothelial", 'VEA': "Vascular endothelial",
                         'VEAT': "Vascular endothelial",
                         'LE1': "Lymphatic endothelial", 'LE2': "Lymphatic endothelial",
                         'CD4_naive': "Lymphoid", 'CD4_Th': "Lymphoid",
                         'CD8_Tem': "Lymphoid", 'CD8_Trm': "Lymphoid",
                         'CD8_Tc1': "Lymphoid",
                         'NKT': "Lymphoid", 'NK': "Lymphoid",
                         'ILC': "Lymphoid",
                         'B_naive': "Lymphoid", 'B_mem_unswitched': "Lymphoid", 'B_mem_switched': "Lymphoid",
                         'Plasma_cell': "Lymphoid",
                         'Macro': "Myeloid", 'Macro-lipo': 'Myeloid',
                         'DC': 'Myeloid',
                         'Doublet': 'Doublet',
                         'stripped_nuclei': 'stripped_nuclei'}

level2_to_level0_dict = {'LASP1': "Epithelial", 'LASP2': "Epithelial",
                         'LASP3': "Epithelial", 'LASP4': "Epithelial",
                         'LASP5': "Epithelial",
                         'LHS1': "Epithelial", 'LHS2': "Epithelial",
                         'LHS3': "Epithelial",
                         'BMYO1': "Epithelial", 'BMYO2': "Epithelial",
                         'DDC1': "Epithelial", 'DDC2': "Epithelial",
                         'FB1': "Stroma", 'FB2': "Stroma", 'FB3': "Stroma", 'FB4': "Stroma",
                         'FB5': "Stroma",
                         'PV1': "Stroma", 'PV2': "Stroma", 'PV3': "Stroma", 'PV4': "Stroma",
                         'PV5': "Stroma",
                         'VEV': "Stroma", 'VEC': "Stroma", 'VEA': "Stroma",
                         'VEAT': "Stroma",
                         'LE1': "Stroma", 'LE2': "Stroma",
                         'CD4_naive': "Immune", 'CD4_Th': "Immune",
                         'CD8_Tem': "Immune", 'CD8_Trm': "Immune",
                         'CD8_Tc1': "Immune",
                         'NKT': "Immune", 'NK': "Immune",
                         'ILC': "Immune",
                         'B_naive': "Immune", 'B_mem_unswitched': "Immune", 'B_mem_switched': "Immune",
                         'Plasma_cell': "Immune",
                         'Macro': "Immune", 'Macro-lipo': 'Immune',
                         'DC': 'Immune',
                         'Doublet': 'Doublet',
                         'stripped_nuclei': 'stripped_nuclei'}

metadata_all.level1 = metadata_all.level2.map(level2_to_level1_dict)
metadata_all.level0 = metadata_all.level2.map(level2_to_level0_dict)

#cells labelled as doublets by Daniel are labelled so in the level0,1 annotations
metadata_all.level0[metadata_all.level0.isna()] = 'Doublet'
metadata_all.level1[metadata_all.level1.isna()] = 'Doublet'

os.makedirs('/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/labelling/output/data/', exist_ok=True)
metadata_all.to_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2023-06-21/scvi/labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2023-06-19.csv')
