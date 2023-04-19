'''
Pypeliner single-sample quality control workflow for 10x scRNA-seq data.
Will create sce and pD objects for plotting QC and later normalization.
Author: Austin Reed <adr44@cam.ac.uk>
'''

import argparse
import os
import yaml
import pypeliner
import pypeliner.managed as mgd


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


def create_input_file_dict(sample_ids, extensions, input_dir):
    '''Create input dictionaries for each file type and sample.'''
    file_dict = {}

    for file_category in extensions.keys():
        file_dict[file_category] = {}

        for file_type, file_extension in extensions[file_category].items():
            file_dict[file_category][file_type] = {}

            for sample_id in sample_ids:
                file_dict[file_category][file_type][sample_id] = os.path.join(input_dir,
                                                                              '{}{}'.format(sample_id, file_extension))

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


def create_qualitycontrol_workflow(config):
    '''Create Pypeliner workflow for completing droplet quality control.'''
    input_extensions = {'downsample': {'barcodes': '.downsample_barcodes_raw.tsv.gz',
                                       'features': '.downsample_features.tsv.gz',
                                       'counts': '.downsample_matrix_raw.mtx.gz'}}

    # final_qc_data extensions are made to match those used in downsampling so that the normalisation scripts
    # work immediately.
    output_extensions = {'sce': {'sce': '.sce.rds'},
                         'sce.alldrops': {'sce.alldrops': '.sce.alldrops.rds'},
                         'pD': {'pD': '.pD.csv'},
                         'pD.alldrops': {'pD.alldrops': '.pD.alldrops.csv'},
                         'final_qc_data': {'barcodes_out': '.downsample_barcodes_filtered.tsv.gz',
                                           'features_out': '.downsample_features.tsv.gz',
                                           'counts_out': '.downsample_matrix_filtered.mtx.gz'},
                         'plots': {'cellQC': '.qualitycontrol.pdf',
                                   'dropsQC': '.alldrops.pdf'}}

    input_dict = create_input_file_dict(config['sample_ids'], input_extensions, config['input_dir'])

    output_dict = create_output_file_dict(config['sample_ids'], output_extensions, config['output_dir'])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=config['sample_ids']
    )

    qualitycontrol_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'qualitycontrol.R')

    workflow.commandline(
        name='qualitycontrol',
        ctx={'mem': 10000, 'ncpus': 1, 'num_retry': 0},
        axes=('sample_id',),
        args=(
            'Rscript', qualitycontrol_script,
            '--sample_id', mgd.InputInstance('sample_id'),
            '--barcode_file', mgd.InputFile('barcodes', 'sample_id', fnames=input_dict['downsample']['barcodes']),
            '--feature_file', mgd.InputFile('features', 'sample_id', fnames=input_dict['downsample']['features']),
            '--count_matrix', mgd.InputFile('counts', 'sample_id', fnames=input_dict['downsample']['counts']),
            '--mito_file', config['mito_file'],
            '--tf_file', config['tf_file'],
            '--metadata_file', config['master_table'],
            '--lower', config['lower'],
            '--niters', config['niters'],
            '--sig', config['sig'],
            '--umi_nmads', config['umi_nmads'],
            '--detected_nmads', config['detected_nmads'],
            '--mt_nmads', config['mt_nmads'],
            '--min_count_depth', config['min_count_depth'],
            '--min_num_genes', config['min_num_genes'],
            '--max_mito', config['max_mito'],
            '--random_seed', config['random_seed'],
            '--out_sce', mgd.OutputFile('sce', 'sample_id', fnames=output_dict['sce']['sce']),
            '--out_pD', mgd.OutputFile('pD', 'sample_id', fnames=output_dict['pD']['pD']),
            '--out_sce.alldrops', mgd.OutputFile('sce.alldrops', 'sample_id', fnames=output_dict['sce.alldrops']['sce.alldrops']),
            '--out_pD.alldrops', mgd.OutputFile('pD.alldrops', 'sample_id', fnames=output_dict['pD.alldrops']['pD.alldrops']),
            '--features_out_pw', mgd.OutputFile('features_out', 'sample_id', fnames=output_dict['final_qc_data']['features_out']),
            '--barcodes_out_pw', mgd.OutputFile('barcodes_out', 'sample_id', fnames=output_dict['final_qc_data']['barcodes_out']),
            '--counts_out_pw', mgd.OutputFile('counts_out', 'sample_id', fnames=output_dict['final_qc_data']['counts_out'])
        )
    )

    return workflow


def main():
    args = parse_args()

    pyp = pypeliner.app.Pypeline(modules=[], config=args)

    config = parse_config(args['config_file'])

    workflow = create_qualitycontrol_workflow(config)

    pyp.run(workflow)


if __name__ == '__main__':
    main()
