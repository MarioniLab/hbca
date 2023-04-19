'''Make the stacked barchat for celltypes over the patients'''

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

import os

plt.rcParams['pdf.fonttype'] = 42


##load Data
metadata_all = pd.read_csv('/nfs/research/marioni/areed/projects/hbca/clustering/2022-04-05/scvi_new/round1_labelling/output/data/metadata_scanpy_HBCA_scVI_processing_date_2022-11-18.csv')
metadata_all = metadata_all.set_index('cellID')


#figure params
sns.set(context='talk',
        style='ticks',
        font='Helvetica',
        rc={'figure.titlesize': 15,
            'axes.titlesize': 15,
            'axes.labelsize': 15,
            'xtick.labelsize': 15,
            'ytick.labelsize': 15,
            'legend.fontsize': 15,
            'legend.title_fontsize': 15,
            'axes.linewidth': 1.3,
            'xtick.major.width': 1.3,
            'ytick.major.width': 1.3})


### Make figure
fig = plt.figure(figsize=(18, 9))

ax = fig.gca()

metadata_sub = metadata_all[['patientID', 'level1']]
metadata_all.plot.bar(x='patientID', stacked=True, ax=ax, edgecolor='none', width=1, legend=False, color=colour_dict)

cell_type_labels = colour_dict.keys()
labels = [x.replace('_', ' ').capitalize() for x in cell_type_labels]
cols = [colour_dict[x] for x in cell_type_labels]

sample_labels = df_metadata.loc[df_metadata['sample_id'].isin(df_pivot.index), 'label'].tolist()

ax.set_ylabel('Number of cells')

if label_metadata:
    ax.set_xticklabels(sample_labels, rotation=45, ha='right', fontdict={'fontsize': 10})
else:
    ax.set_xlabel('Sample')

lg = add_legend(ax, labels, cols, 3, type='square', location='upper left', markersize=12, frameon=True)

sns.despine()
plt.tight_layout()

pdf.savefig()
plt.close()