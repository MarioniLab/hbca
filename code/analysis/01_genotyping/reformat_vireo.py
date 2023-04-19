'''Read, reformat, and output Vireo files.

import os

vireo_dir = '/Users/steif01/Documents/projects/hbca/tmp_analysis/vireo/SLX-20005-SIGAE11/original'
vireo_out = '/Users/steif01/Documents/projects/hbca/tmp_analysis/vireo/SLX-20005-SIGAE11/reformat'

known_id = 'SLX-20005-SIGAE10'
sample_id = 'SLX-20005-SIGAE11'

donor_gt_file = os.path.join(vireo_dir, 'GT_donors.vireo.vcf.gz')
donor_id_file = os.path.join(vireo_dir, 'donor_ids.tsv')
singlet_prob_file = os.path.join(vireo_dir, 'prob_singlet.tsv.gz')
doublet_prob_file = os.path.join(vireo_dir, 'prob_doublet.tsv.gz')
summary_file = os.path.join(vireo_dir, 'summary.tsv')
log_file = os.path.join(vireo_dir, '_log.txt')

out_donor_gt = os.path.join(vireo_out, '{}.vireo_genotypes.vcf.gz'.format(sample_id))
out_donor_id = os.path.join(vireo_out, '{}.vireo_assign.csv'.format(sample_id))
out_singlet_prob = os.path.join(vireo_out, '{}.vireo_singlet_prob.csv'.format(sample_id))
out_doublet_prob = os.path.join(vireo_out, '{}.vireo_doublet_prob.csv'.format(sample_id))
out_summary = os.path.join(vireo_out, '{}.vireo_summary.csv'.format(sample_id))
out_log = os.path.join(vireo_out, '{}.vireo_log.txt'.format(sample_id))

'''
import gzip
import argparse
import shutil
import pandas as pd


def parse_args():
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser()

    parser.add_argument('--known_id',
                        help='''Identifier for the genotyped spike-in sample.''')

    parser.add_argument('--sample_id',
                        help='''Identifier for the donor sample.''')

    parser.add_argument('--donor_gt_file',
                        help='''Path to Vireo donor genotype file in .vcf.gz format.''')

    parser.add_argument('--donor_id_file',
                        help='''Path to Vireo cell assignment file in .tsv format.''')

    parser.add_argument('--singlet_prob_file',
                        help='''Path to Vireo singlet probability file in .tsv.gz format.''')

    parser.add_argument('--doublet_prob_file',
                        help='''Path to Vireo doublet probability file in .tsv.gz format.''')

    parser.add_argument('--summary_file',
                        help='''Path to Vireo summary file in .tsv format.''')

    parser.add_argument('--log_file',
                        help='''Path to Vireo log file in .txt format.''')

    parser.add_argument('--out_donor_gt',
                        help='''Path to reformatted donor genotype output file in .vcf.gz format.''')

    parser.add_argument('--out_donor_id',
                        help='''Path to reformatted cell assignment output file in .csv format.''')

    parser.add_argument('--out_singlet_prob',
                        help='''Path to reformatted singlet probability output file in .csv format.''')

    parser.add_argument('--out_doublet_prob',
                        help='''Path to reformatted doublet probability output file in .csv format.''')

    parser.add_argument('--out_summary',
                        help='''Path to reformatted summary output file in .csv format.''')

    parser.add_argument('--out_log',
                        help='''Path to renamed log file in .txt format.''')

    args = parser.parse_args()

    return args


def parse_vireo_vcf(known_id, sample_id, vcf_file):
    '''Parse and reformat Vireo VCF output file with genotypes for each sample.'''
    df = pd.read_csv(vcf_file, compression='gzip', sep='\t', comment='#', header=None,
                     names=['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', known_id, sample_id])

    new_header = [x.upper() for x in df.columns.values.tolist()]

    header_line = '#{}\n'.format('\t'.join(new_header))

    header = []

    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                header.append(line)
            elif line.startswith('#'):
                header.append(header_line)

    return header, df


def write_vireo_vcf(header, df_vcf, out_file):
    '''Write the reformatted VCF file, including modified header.'''
    with gzip.open(out_file, 'wb') as f:
        for line in header:
            f.write(line.encode())

    df_vcf.to_csv(out_file, compression='gzip', mode='a', sep='\t', header=False)


def parse_vireo_cell_assignment_file(sample_id, donor_id_file):
    '''Parse and reformat Vireo output file with cell donor assignment.'''
    df = pd.read_csv(donor_id_file, sep='\t')

    df.rename(columns={'cell': 'cell_barcode'}, inplace=True)

    df = df.replace('donor1', sample_id)

    df['best_doublet'] = df['best_doublet'].str.replace('donor1', sample_id).str.replace(',', ';')

    return df


def parse_vireo_probability_files(sample_id, singlet_prob_file, doublet_prob_file):
    '''Parse and reformat Vireo singlet and doublet probability files.'''
    df_singlet_prob = pd.read_csv(singlet_prob_file, compression='gzip', sep='\t')
    df_doublet_prob = pd.read_csv(doublet_prob_file, compression='gzip', sep='\t')

    df_singlet_prob.rename(columns={'cell': 'cell_barcode',
                                    'donor1': sample_id}, inplace=True)

    df_doublet_prob.rename(columns={'cell': 'cell_barcode'}, inplace=True)

    df_doublet_prob.columns = [x.replace('donor1', sample_id).replace(',', ';') for x in df_doublet_prob.columns]

    return df_singlet_prob, df_doublet_prob


def parse_vireo_summary(sample_id, summary_file):
    '''Parse and reformat Vireo summary output file.'''
    df = pd.read_csv(summary_file, sep='\t')

    df.columns = ['assignment', 'num_cells']

    df['assignment'] = df['assignment'].replace('donor1', sample_id)

    return df


def reformat_vireo_files(known_id, sample_id, donor_gt_file, donor_id_file,
                         singlet_prob_file, doublet_prob_file, summary_file):
    '''Reformat the Vireo output files to include the sample ID instead of 'donor1'.'''
    header, df_vcf = parse_vireo_vcf(known_id, sample_id, donor_gt_file)

    df_assign = parse_vireo_cell_assignment_file(sample_id, donor_id_file)

    df_singlet, df_doublet = parse_vireo_probability_files(sample_id, singlet_prob_file, doublet_prob_file)

    df_summary = parse_vireo_summary(sample_id, summary_file)

    return header, df_vcf, df_assign, df_singlet, df_doublet, df_summary


def write_reformatted_vireo_files(known_id, sample_id, donor_gt_file, donor_id_file,
                                  singlet_prob_file, doublet_prob_file, summary_file, log_file,
                                  out_donor_gt, out_donor_id, out_singlet_prob, out_doublet_prob, out_summary, out_log):
    '''Parse and reformat the Vireo output files with sample IDs, then write to new files (csv instead of tsv).'''
    header, df_vcf, df_assign, \
    df_singlet, df_doublet, df_summary = reformat_vireo_files(known_id, sample_id, donor_gt_file, donor_id_file,
                                                              singlet_prob_file, doublet_prob_file, summary_file)

    write_vireo_vcf(header, df_vcf, out_donor_gt)

    df_assign.to_csv(out_donor_id, index=False)

    df_singlet.to_csv(out_singlet_prob, index=False)

    df_doublet.to_csv(out_doublet_prob, index=False)

    df_summary.to_csv(out_summary, index=False)

    shutil.copy(log_file, out_log)


def main():
    args = parse_args()

    write_reformatted_vireo_files(args.known_id, args.sample_id,
                                  args.donor_gt_file, args.donor_id_file, args.singlet_prob_file,
                                  args.doublet_prob_file, args.summary_file, args.log_file,
                                  args.out_donor_gt, args.out_donor_id, args.out_singlet_prob,
                                  args.out_doublet_prob, args.out_summary, args.out_log)


if __name__ == '__main__':
    main()
