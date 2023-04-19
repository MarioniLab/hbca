'''Plot Vireo sample assignment UMAP and cell genotype heatmap.'''
import gzip
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.io import mmread
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import add_legend


def parse_args():
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser()

    parser.add_argument('--sample_ids', nargs='+',
                        help='''List of sample identifiers.''')

    parser.add_argument('--summary_files', nargs='+',
                        help='''Path to Vireo summary file for each sample in .csv format.''')

    parser.add_argument('--genotype_files', nargs='+',
                        help='''Path to Vireo genotype file for each sample in .vcf format.''')

    parser.add_argument('--assign_files', nargs='+',
                        help='''Path to Vireo assignment file for each sample in .csv format.''')

    parser.add_argument('--out_dir',
                        help='''Path to output directory.''')

    args = parser.parse_args()

    return args


# Vireo colour palettes


def get_genotype_colours(sample_id):
    '''Get colour dictionary for the sample, spike-in, doublet and unassigned cells.'''
    tableau_miller_stone = ['#4f6980', '#849db1', '#a2ceaa', '#638b66', '#bfbb60',
                            '#f47942', '#fbb04e', '#b66353', '#d7ce9f', '#b9aa97', '#7e756d']

    colour_dict = {sample_id: tableau_miller_stone[1],
                   'Spike-in': tableau_miller_stone[6],
                   'Doublet': tableau_miller_stone[5],
                   'Unassigned': tableau_miller_stone[10]}

    return colour_dict


def get_patient_colours(num_cols):
    '''Get list of colours corresponding to the patient IDs.'''
    new_tableau_20 = ['#4e79a7', '#a0cbe8', '#f28e2b', '#ffbe7d', '#59a14f',
                      '#8cd17d', '#b6992d', '#f1ce63', '#499894', '#86bcb6',
                      '#e15759', '#ff9d9a', '#79706e', '#bab0ac', '#d37295',
                      '#fabfd2', '#b07aa1', '#d4a6c8', '#9d7660', '#d7b5a6']

    bardot_sunny_day_hike = ['#005dc6', '#00b2d7', '#7ca001', '#009b7c', '#d33a02',
                             '#027d59', '#bacc47', '#5e89ca', '#cddaf0', '#fb9361',
                             '#938e7e', '#6e582f', '#e8909e', '#c45762', '#fe9e88']

    bardot_midnight_produce = ['#030303', '#3edcb6', '#ff9509', '#d9871d', '#febb2c',
                               '#1a8dbf', '#0b656c', '#c20011', '#7b0406', '#50062b',
                               '#9fbebe', '#b1ab6b', '#7c6619', '#cecb1d']

    bardot_rosy_days = ['#f7dfd7', '#e7b6b0', '#ed8a8e', '#b33a6d', '#707a42',
                        '#5da4b9', '#88dad2', '#b5e5df', '#e5c86e', '#cd9d3c']

    cols = new_tableau_20 + bardot_sunny_day_hike + bardot_midnight_produce + bardot_rosy_days

    return cols[:num_cols]


# Vireo scatter plots


def plot_vireo_sample_scatter(sample_id, df, pdf, spike_in='SLX-20005-20446_SIGAE10'):
    '''Plot Vireo cell genotypes on a UMAP scatter plot.'''
    donors_sample = pd.Series(df['donor_id'].unique()).tolist()

    donors_all = [sample_id, 'doublet', spike_in, 'unassigned']

    donor_all_labels = pd.Series(donors_all).replace(spike_in, 'spike-in'.capitalize()).replace(
                                                 'doublet', 'doublet'.capitalize()).replace(
                                                 'unassigned', 'unassigned'.capitalize()).tolist()

    donor_sample_labels = [y for x, y in zip(donors_all, donor_all_labels) if x in donors_sample]

    colour_dict = get_genotype_colours(sample_id)

    donor_cols = [colour_dict[x] for x in donor_all_labels]

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

    fig = plt.figure(figsize=(11.25, 10))

    ax = fig.gca()

    labels = []

    for donor, col, label in zip(donors_all, donor_cols, donor_all_labels):
        if donor in donors_sample:
            df_donor = df[df['donor_id'] == donor]
            ax.scatter(df_donor['umap1'], df_donor['umap2'], color=col, edgecolor='none', alpha=0.6, s=16)
            label_counts = '{} (n={})'.format(label, df['donor_id'].value_counts()[donor])
            labels.append(label_counts)
        else:
            label_counts = '{} (n={})'.format(label, 0)
            labels.append(label_counts)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.tick_params(axis='both', which='both', length=0)

    add_legend(ax, labels, donor_cols, 1,
               type='circle', location='upper right', markersize=12, frameon=True)

    sns.despine()
    plt.tight_layout()

    pdf.savefig()
    plt.close()


def plot_vireo_scatter(sample_ids, vireo_assign_dict, cellranger_umap_dict, out_file,
                       spike_in='SLX-20005-20446_SIGAE10'):
    '''Plot Vireo cell genotypes on UMAP using Cellranger dimensionality reduction coordinates.'''
    with PdfPages(out_file) as pdf:
        for sample_id in sample_ids:
            df_assign = pd.read_csv(vireo_assign_dict[sample_id])
            df_umap = pd.read_csv(cellranger_umap_dict[sample_id])
            df_umap.columns = ['cell_barcode', 'umap1', 'umap2']
            df = df_assign.merge(df_umap)
            plot_vireo_sample_scatter(sample_id, df, pdf, spike_in=spike_in)


# Vireo bar plots


def combine_vireo_summary_files(vireo_summary_dict, out_file, spike_in='SLX-20005-20446_SIGAE10'):
    '''Combine the Vireo summary files across all samples.'''
    df = pd.DataFrame()

    for sample_id in vireo_summary_dict.keys():
        df_summary = pd.read_csv(vireo_summary_dict[sample_id])
        df_summary['sample_id'] = sample_id
        df_summary['assignment'] = df_summary['assignment'].replace(sample_id, 'Donor')
        df = pd.concat([df, df_summary])

    df['assignment'] = df['assignment'].replace(spike_in, 'Spike-in')
    df['assignment'] = df['assignment'].replace('doublet', 'Doublet')
    df['assignment'] = df['assignment'].replace('unassigned', 'Unassigned')

    df_pivot = df.pivot_table(columns='assignment', index='sample_id', values='assignment')
    df_pivot = df_pivot.fillna(0).astype(int).rename_axis(None, axis=1).reset_index()

    df_pivot.to_csv(out_file, index=False)


def plot_vireo_stacked_bar_summary(df_pivot, pdf):
    '''Stacked bar plot showing number of cells assigned to each category.'''
    sns.set(context='talk',
            style='ticks',
            font='Helvetica',
            rc={'figure.titlesize': 15,
                'axes.titlesize': 15,
                'axes.labelsize': 15,
                'xtick.labelsize': 8,
                'ytick.labelsize': 15,
                'legend.fontsize': 15,
                'legend.title_fontsize': 15,
                'axes.linewidth': 1.3,
                'xtick.major.width': 1.3,
                'ytick.major.width': 1.3})

    cols = get_genotype_colours('Donor')

    assignment_labels = df_pivot.columns.values.tolist()
    assignment_cols = [cols[x] for x in assignment_labels]

    assignment_cmap = colors.ListedColormap(assignment_cols)

    index_total = df_pivot.sum(axis=1).sort_values(ascending=False).index

    df_pivot_reorder = df_pivot.loc[index_total]

    fig = plt.figure(figsize=(24, 9))

    ax = fig.gca()

    df_pivot_reorder.plot.bar(stacked=True, ax=ax, edgecolor='none', width=1, legend=False, colormap=assignment_cmap)

    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Number of cells')

    ax = add_legend(ax, assignment_labels, assignment_cols, 1,
                    type='square', location='upper right', markersize=12, frameon=True)

    sns.despine()
    plt.tight_layout()

    pdf.savefig()
    plt.close()


def plot_vireo_single_bar_summary(df_pivot, pdf, assign_column='Spike-in'):
    '''Bar plot showing number of cells of a particular genotype category identified in each sample.'''
    sns.set(context='talk',
            style='ticks',
            font='Helvetica',
            rc={'figure.titlesize': 15,
                'axes.titlesize': 15,
                'axes.labelsize': 15,
                'xtick.labelsize': 8,
                'ytick.labelsize': 15,
                'legend.fontsize': 15,
                'legend.title_fontsize': 15,
                'axes.linewidth': 1.3,
                'xtick.major.width': 1.3,
                'ytick.major.width': 1.3})

    cols = get_genotype_colours('Donor')

    assignment_col = cols[assign_column]

    assignment_cmap = colors.ListedColormap(assignment_col)

    df_bar = df_pivot[assign_column].sort_values(ascending=False)

    fig = plt.figure(figsize=(24, 9))

    ax = fig.gca()

    df_bar.plot.bar(ax=ax, edgecolor='none', width=1, legend=False, colormap=assignment_cmap)

    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Number of {} cells'.format(assign_column.lower()))

    sns.despine()
    plt.tight_layout()

    pdf.savefig()
    plt.close()


def plot_vireo_bar(combined_summary_file, out_file):
    '''Plot Vireo bar plots showing distribution of genotype assignments across samples.'''
    df_pivot = pd.read_csv(combined_summary_file)

    df_pivot.set_index('sample_id', inplace=True)

    with PdfPages(out_file) as pdf:
        plot_vireo_stacked_bar_summary(df_pivot, pdf)
        plot_vireo_single_bar_summary(df_pivot, pdf, assign_column='Donor')
        plot_vireo_single_bar_summary(df_pivot, pdf, assign_column='Spike-in')
        plot_vireo_single_bar_summary(df_pivot, pdf, assign_column='Doublet')
        plot_vireo_single_bar_summary(df_pivot, pdf, assign_column='Unassigned')


# Vireo heatmap plots


def read_vireo_vcf_file(donor_gt_file):
    '''Read the Vireo donor genotype file and extract the genotype state for each donor.'''
    header = None

    with gzip.open(donor_gt_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                pass
            elif line.startswith('#'):
                header = line.strip('\n').lstrip('#').split('\t')

    donors = header[9:]

    header = [x.lower() for x in header[0:9]] + donors

    df_donor_gt = pd.read_csv(donor_gt_file, compression='gzip', sep='\t', comment='#', header=None, names=header)

    for donor in donors:
        donor_columns = ['{}_GT'.format(donor), '{}_AD'.format(donor), '{}_DP'.format(donor), '{}_PL'.format(donor)]
        df_donor_gt[donor_columns] = df_donor_gt[donor].str.split(':', expand=True)

    df_donor_gt.insert(0, 'chrom_pos', df_donor_gt['chrom'] + '_' + df_donor_gt['pos'].astype('str'))

    return donors, df_donor_gt


def combine_vireo_vcf_files(donor_gt_dict, out_file, spike_in='SLX-20005-20446_SIGAE10'):
    '''Read the Vireo VCF output for multiple donors and merge the dataframes.'''
    df = pd.DataFrame()
    i = 1

    for sample_id in donor_gt_dict.keys():
        donors, df_donor_gt = read_vireo_vcf_file(donor_gt_dict[sample_id])
        if len(df) == 0:
            # drop the spike-in columns because these only represent reads from the spike-in cells found in this patient
            drop_columns = [spike_in] + ['{}_{}'.format(spike_in, col_str) for col_str in ['GT', 'AD', 'DP', 'PL']]
            df_donor_gt.drop(drop_columns, axis=1, inplace=True)
            df = df_donor_gt
        else:
            subset_columns = ['chrom_pos', sample_id] + \
                             ['{}_{}'.format(sample_id, col_str) for col_str in ['GT', 'AD', 'DP', 'PL']]
            df = df.merge(df_donor_gt[subset_columns], how='outer', on='chrom_pos')

        print('{}/{}: {}'.format(i, len(donor_gt_dict.keys()), sample_id))
        i += 1

    df.to_csv(out_file)


def format_combined_vireo_vcf(combined_vcf_file):
    '''Format the combined Vireo VCF file (saved in CSV format without a the VCF header).'''
    df = pd.read_csv(combined_vcf_file)

    df.set_index('chrom_pos', inplace=True)

    donor_list = [x.replace('_GT', '') for x in df.columns.values.tolist() if x.endswith('_GT')]

    gt_columns = ['{}_GT'.format(x) for x in donor_list]
    dp_columns = ['{}_DP'.format(x) for x in donor_list]
    ad_columns = ['{}_AD'.format(x) for x in donor_list]

    df_gt = df[gt_columns].copy()
    df_gt.replace('0/0', 0, inplace=True)
    df_gt.replace('1/0', 1, inplace=True)
    df_gt.replace('1/1', 2, inplace=True)

    df_dp = df[dp_columns].copy()
    df_ad = df[ad_columns].copy()

    df_dp = df_dp.astype(int)
    df_ad = df_ad.astype(int)

    df_gt.columns = donor_list
    df_dp.columns = donor_list
    df_ad.columns = donor_list

    return df, df_gt, df_dp, df_ad


def informative_row(df_row):
    '''Return true if at least one sample differs in the row genotype.'''
    return len(df_row.unique()) > 1


def select_informative_variants(donors, df_gt, df_dp, df_ad, min_depth=10):
    '''Subset genotypes that differ between samples and have the minimum specified coverage per sample.'''
    df_index = pd.DataFrame(index=df_dp.index)

    df_index['informative'] = df_gt.apply(informative_row, axis=1)

    # each sample must have coverage greater than min_depth
    df_index['min_depth'] = (df_dp >= min_depth).sum(axis=1) == len(donors)

    df_index['both'] = df_index['informative'] & df_index['min_depth']

    df_gt_inform = df_gt[df_index['both']]
    df_dp_inform = df_dp[df_index['both']]
    df_ad_inform = df_ad[df_index['both']]

    return df_gt_inform, df_dp_inform, df_ad_inform


def plot_vireo_genotype_heatmap(df_gt_inform, master_table, pdf):
    '''Plot heatmap for the donor genotypes, clustering samples by donor.'''
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

    genotype_cols = ['#97cfd0', '#f7c480', '#f3a546']

    genotype_cmap = colors.ListedColormap(genotype_cols)

    # assign colours based on patient

    metadata = pd.read_csv(master_table)

    patient_ids = [x for x in metadata['patient_id'].unique() if x is not np.nan]

    patient_cols = get_patient_colours(len(patient_ids))

    df_map = metadata[['sample_id', 'patient_id']].set_index('sample_id')

    df_map_ordered = df_map.loc[df_gt_inform.columns.values.tolist()]

    df_map_ordered['Patient'] = df_map_ordered['patient_id']

    df_map_ordered['Patient'] = df_map_ordered['Patient'].replace(patient_ids, patient_cols)

    patient_colours = df_map_ordered['Patient']

    # plot heatmap

    plt.rcParams['pdf.fonttype'] = 42

    cm = sns.clustermap(df_gt_inform,
                        cmap=genotype_cmap,
                        figsize=(60, 20),
                        dendrogram_ratio=0.1, col_colors=patient_colours)

    ax = cm.ax_heatmap
    ax.tick_params(right=False)
    ax.set_ylabel('SNP (n={})'.format(len(df_gt_inform)), rotation=-90)
    ax.set_yticklabels([])

    cm.cax.set_visible(False)
    #cm.ax_row_dendrogram.set_visible(False)

    lg = add_legend(ax, ['1/1', '0/1', '0/0'],
                    [x for x in reversed(genotype_cols)],
                    1, type='square', location='upper left', markersize=12,
                    frameon=True, markerfirst=True, title=None)

    plt.tight_layout()

    pdf.savefig()
    plt.close()


def plot_vireo_heatmap(combined_vcf_file, master_table, out_file):
    '''Plot genotype heatmap for all patient samples.'''
    df, df_gt, df_dp, df_ad = format_combined_vireo_vcf(combined_vcf_file)

    donors = df_gt.columns.values.tolist()

    df_gt_inform, df_dp_inform, df_ad_inform = select_informative_variants(donors, df_gt, df_dp, df_ad, min_depth=10)

    with PdfPages(out_file) as pdf:
        plot_vireo_genotype_heatmap(df_gt_inform, master_table, pdf)


'''
# Not used

def read_vartrix_allele_counts(alt_mat_file, cov_mat_file, variants_file):
    #Read the vartrix per cell allele count matrices.
    alt_counts = mmread(alt_mat_file)
    cov_counts = mmread(cov_mat_file)
    variants = pd.read_csv(variants_file, sep='\t', header=None, names=['chrom_pos'])

    return alt_counts, cov_counts, variants


def subset_vartrix_frequency_matrix(df, alt_counts, ref_counts, variant_subset):
    #Subset the Vartrix matrices with the specified variants in Pandas Series format.
    tableau_miller_stone = ['#4f6980', '#849db1', '#a2ceaa', '#638b66', '#bfbb60',
                            '#f47942', '#fbb04e', '#b66353', '#d7ce9f', '#b9aa97', '#7e756d']

    cols = [tableau_miller_stone[i] for i in [1, 6, 5, 10]]

    variant_indices = variant_subset.index
    alt_counts_subset = alt_counts.tocsr()[variant_indices,:].toarray()
    ref_counts_subset = ref_counts.tocsr()[variant_indices,:].toarray()

    freq_mat = alt_counts_subset / (alt_counts_subset + ref_counts_subset)

    freq_mat_no_nan = np.nan_to_num(freq_mat, nan=-1)

    sns.clustermap(freq_mat_no_nan, method='average', mask=np.isnan(freq_mat), vmin=0, vmax=1)

    gt_mat = freq_mat
    gt_mat[(gt_mat > 0) & (gt_mat < 1)] = 3
    gt_mat[gt_mat == 1] = 2
    gt_mat[gt_mat == 3] = 1

    gt_mat_no_nan = np.nan_to_num(gt_mat, nan=-1)

    donors = df['donor_id'].unique()

    donor_col_dict = dict(zip(donors, cols))

    column_cols = df['donor_id'].map(donor_col_dict)

    freq_cmap = sns.color_palette('RdYlBu_r', as_cmap=True)

    sns.clustermap(pd.DataFrame(gt_mat_no_nan), col_colors=column_cols, mask=np.isnan(freq_mat), vmin=0, vmax=2, cmap=freq_cmap)
'''


def main():
    args = parse_args()
    pass


if __name__ == '__main__':
    main()
