'''
Pypeliner workflow for demultiplexing cells based on genotype.
Author: Adi Steif <as2886@cam.ac.uk>
'''
import argparse
import os
import yaml
import pypeliner
import pypeliner.managed as mgd
from run_vireo import run_vireo
from plot_vireo import combine_vireo_summary_files, combine_vireo_vcf_files, \
     plot_vireo_scatter, plot_vireo_bar, plot_vireo_heatmap


def parse_args():
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser()

    pypeliner.app.add_arguments(parser)

    parser.add_argument('--config_file',
                        help='''Path to YAML config file.''')

    args = vars(parser.parse_args())

    return args


def parse_config(config_file):
    '''Parse the configuration file.'''
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    return config


def create_input_file_dict(sample_ids, extensions, input_dir, sub_dir='outs'):
    '''Create input dictionaries for each file type and sample.'''
    file_dict = {}

    for file_category in extensions.keys():
        file_dict[file_category] = {}

        for file_type, file_extension in extensions[file_category].items():
            file_dict[file_category][file_type] = {}

            for sample_id in sample_ids:
                input_path = os.path.join(input_dir, '{}'.format(sample_id), sub_dir)
                file_dict[file_category][file_type][sample_id] = os.path.join(input_path, file_extension)

    return file_dict


def create_output_file_dict(sample_ids, extensions, output_dir, sub_dir='output'):
    '''Create output dictionaries for each file type and sample.'''
    file_dict = {}

    for file_category in extensions.keys():
        file_dict[file_category] = {}

        for file_type, file_extension in extensions[file_category].items():
            file_dict[file_category][file_type] = {}

            for sample_id in sample_ids:
                file_dict[file_category][file_type][sample_id] = os.path.join(output_dir,
                                                                           sub_dir,
                                                                           file_category,
                                                                           '{}{}'.format(sample_id, file_extension))

    return file_dict


def create_genotyping_workflow(config):
    '''Create Pypeliner workflow for demultiplexing cells based on genotype. '''
    input_extensions = {'cellranger': {'bam': 'possorted_genome_bam.bam',
                                       'barcodes': os.path.join('filtered_feature_bc_matrix', 'barcodes.tsv.gz'),
                                       'umap': os.path.join('analysis', 'umap', '2_components', 'projection.csv')}}

    output_extensions = {'vartrix': {'alt_mat': '.vartrix_alt.mtx',
                                     'ref_mat': '.vartrix_ref.mtx',
                                     'variants': '.vartrix_variants.tsv'},
                         'vireo': {'donor_gt': '.vireo_genotypes.vcf.gz',
                                   'donor_id': '.vireo_assign.csv',
                                   'singlet_prob': '.vireo_singlet_prob.csv',
                                   'doublet_prob': '.vireo_doublet_prob.csv',
                                   'summary': '.vireo_summary.csv',
                                   'log': '.vireo_log.txt'}}

    input_dict = create_input_file_dict(config['sample_ids'], input_extensions, config['input_dir'])

    output_dict = create_output_file_dict(config['sample_ids'], output_extensions, config['output_dir'])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=config['sample_ids']
    )

    # run vartrix to extract coverage for each cell

    workflow.commandline(
        name='vartrix',
        ctx={'mem': 32000, 'ncpus': 16, 'num_retry': 0},
        axes=('sample_id',),
        args=(
            'vartrix',
            '--log-level', 'info',
            '--bam', mgd.InputFile('bam', 'sample_id', fnames=input_dict['cellranger']['bam']),
            '--cell-barcodes', mgd.InputFile('barcodes', 'sample_id', fnames=input_dict['cellranger']['barcodes']),
            '--fasta', config['reference_genome'],
            '--vcf', config['known_genotype_vcf'],
            '--mapq', config['min_mapq'],
            '--umi',
            '--scoring-method', 'coverage',
            '--out-matrix', mgd.OutputFile('vartrix_alt_mat', 'sample_id', fnames=output_dict['vartrix']['alt_mat']),
            '--ref-matrix', mgd.OutputFile('vartrix_ref_mat', 'sample_id', fnames=output_dict['vartrix']['ref_mat']),
            '--out-variants', mgd.OutputFile('vartrix_variants', 'sample_id', fnames=output_dict['vartrix']['variants']),
            '--threads', config['num_threads']
        )
    )

    # run vireo to infer genotypes and assign cells to a donor sample

    workflow.transform(
        name='vireo',
        func=run_vireo,
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        axes=('sample_id',),
        args=(
            mgd.InputFile('vartrix_alt_mat', 'sample_id', fnames=output_dict['vartrix']['alt_mat']),
            mgd.InputFile('vartrix_ref_mat', 'sample_id', fnames=output_dict['vartrix']['ref_mat']),
            mgd.InputFile('barcodes', 'sample_id', fnames=input_dict['cellranger']['barcodes']),
            config['known_genotype_vcf'],
            mgd.TempSpace('vireo', 'sample_id'),
            config['num_donors'],
            config['random_seed'],
            config['known_genotype_sample'],
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('vireo_genotypes', 'sample_id', fnames=output_dict['vireo']['donor_gt']),
            mgd.OutputFile('vireo_assign', 'sample_id', fnames=output_dict['vireo']['donor_id']),
            mgd.OutputFile('vireo_singlet', 'sample_id', fnames=output_dict['vireo']['singlet_prob']),
            mgd.OutputFile('vireo_doublet', 'sample_id', fnames=output_dict['vireo']['doublet_prob']),
            mgd.OutputFile('vireo_summary', 'sample_id', fnames=output_dict['vireo']['summary']),
            mgd.OutputFile('vireo_log', 'sample_id', fnames=output_dict['vireo']['log'])
        )
    )

    # combine Vireo summary files

    workflow.transform(
        name='combine_vireo_summary_files',
        func=combine_vireo_summary_files,
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        args=(
            mgd.InputFile('vireo_summary', 'sample_id', fnames=output_dict['vireo']['summary']),
            mgd.OutputFile(os.path.join(config['output_dir'], 'output', 'summary', 'vireo_counts_summary.csv')),
            config['known_genotype_sample']
        )
    )

    # combine Vireo VCF files

    workflow.transform(
        name='combine_vireo_vcf_files',
        func=combine_vireo_vcf_files,
        ctx={'mem': 32000, 'ncpus': 1, 'num_retry': 0},
        args=(
            mgd.InputFile('vireo_genotypes', 'sample_id', fnames=output_dict['vireo']['donor_gt']),
            mgd.OutputFile(os.path.join(config['output_dir'], 'output', 'summary', 'vireo_vcf_summary.csv')),
            config['known_genotype_sample']
        )
    )

    # plot Vireo scatter plots on CellRanger UMAP coordinates

    workflow.transform(
        name='plot_vireo_scatter',
        func=plot_vireo_scatter,
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        args=(
            config['sample_ids'],
            mgd.InputFile('vireo_assign', 'sample_id', fnames=output_dict['vireo']['donor_id']),
            mgd.InputFile('umap', 'sample_id', fnames=input_dict['cellranger']['umap']),
            mgd.OutputFile(os.path.join(config['output_dir'], 'output', 'plots', 'vireo_scatter.pdf')),
            config['known_genotype_sample']
        )
    )

    # plot Vireo bar plots

    workflow.transform(
        name='plot_vireo_bar',
        func=plot_vireo_bar,
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        args=(
            mgd.InputFile(os.path.join(config['output_dir'], 'output', 'summary', 'vireo_counts_summary.csv')),
            mgd.OutputFile(os.path.join(config['output_dir'], 'output', 'plots', 'vireo_bar.pdf'))
        )
    )

    # plot Vireo heatmap

    workflow.transform(
        name='plot_vireo_heatmap',
        func=plot_vireo_heatmap,
        ctx={'mem': 16000, 'ncpus': 1, 'num_retry': 0},
        args=(
            mgd.InputFile(os.path.join(config['output_dir'], 'output', 'summary', 'vireo_vcf_summary.csv')),
            config['master_table'],
            mgd.OutputFile(os.path.join(config['output_dir'], 'output', 'plots', 'vireo_heatmap.pdf'))
        )
    )

    return workflow


def main():
    args = parse_args()

    pyp = pypeliner.app.Pypeline(modules=[], config=args)

    config = parse_config(args['config_file'])

    workflow = create_genotyping_workflow(config)

    pyp.run(workflow)


if __name__ == '__main__':
    main()
