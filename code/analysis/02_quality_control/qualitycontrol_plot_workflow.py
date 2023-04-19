'''
Pypeliner single-sample quality control workflow for 10x scRNA-seq data. Will create summary graphs and prepare data for normalization.
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
    '''Create Pypeliner workflow for plotting droplet quality control.'''

    output_extensions = {'sce': {'sce': '.sce.rds'},
                         'sce.alldrops': {'sce.alldrops': '.sce.alldrops.rds'},
                         'pD': {'pD': '.pD.csv'},
                         'pD.alldrops': {'pD.alldrops': '.pD.alldrops.csv'},
                         'plots': {'cellQC': '.qualitycontrol.pdf',
                                   'dropsQC': '.alldrops.pdf'}}

    output_dict = create_output_file_dict(config['sample_ids'], output_extensions, config['output_dir'])

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=config['sample_ids']
    )

    plot_qualitycontrol_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'plot_qualitycontrol.R')

    workflow.commandline(
        name='plot_qualitycontrol',
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        axes=('sample_id',),
        args=(
            'Rscript', plot_qualitycontrol_script,
            '--sample_id', mgd.InputInstance('sample_id'),
            '--pD_file', mgd.InputFile('pD', 'sample_id', fnames=output_dict['pD']['pD']),
            '--sce.alldrops_file',
            mgd.InputFile('sce.alldrops', 'sample_id', fnames=output_dict['sce.alldrops']['sce.alldrops']),
            '--random_seed', config['random_seed'],
            '--bins', config['bins'],
            '--lower', config['lower'],
            '--out_QC_plots', mgd.OutputFile('cellQC', 'sample_id', fnames=output_dict['plots']['cellQC']),
            '--out_drops_plots', mgd.OutputFile('dropsQC', 'sample_id', fnames=output_dict['plots']['dropsQC'])
        )
    )

    plot_qualitycontrol_joint_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                    'plot_qualitycontrol_joint.R')

    pD_pwd = os.path.join(config['output_dir'], 'output', 'pD')
    out_QC_joint_tissuecondition_pw = os.path.join(config['output_dir'], 'output', 'plots', 'allsamples.qualitycontrol.tissuecondition.pdf')
    out_QC_joint_tissuecondition_clean_pw = os.path.join(config['output_dir'], 'output', 'plots', 'allsamples.qualitycontrol.tissuecondition.clean.pdf')
    out_QC_joint_sampletype_pw = os.path.join(config['output_dir'], 'output', 'plots', 'allsamples.qualitycontrol.sampletype.pdf')
    out_QC_joint_sampletype_clean_pw = os.path.join(config['output_dir'], 'output', 'plots', 'allsamples.qualitycontrol.sampletype.clean.pdf')

    workflow.commandline(
        name='plot_qualitycontrol_joint',
        ctx={'mem': 8000, 'ncpus': 1, 'num_retry': 0},
        args=(
            'Rscript', plot_qualitycontrol_joint_script,
            '--pD_pwd', mgd.InputFile(pD_pwd),
            '--mastertable_pw', config['master_table'],
            '--random_seed', config['random_seed'],
            '--out_QC_plots_joint_tissuecondition', mgd.OutputFile(out_QC_joint_tissuecondition_pw),
            '--out_QC_plots_joint_tissuecondition_clean', mgd.OutputFile(out_QC_joint_tissuecondition_clean_pw),
            '--out_QC_plots_joint_sampletype', mgd.OutputFile(out_QC_joint_sampletype_pw),
            '--out_QC_plots_joint_sampletype_clean', mgd.OutputFile(out_QC_joint_sampletype_clean_pw)
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