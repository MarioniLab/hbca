'''
Run Vireo directing output to a Pypeliner temporary directory, then reformat and rename the output files.
'''
import os
import pypeliner
from reformat_vireo import write_reformatted_vireo_files


def run_vireo(alt_mat_file, ref_mat_file, barcodes, known_genotype_vcf, vireo_tmp_dir, num_donors, random_seed,
              known_id, sample_id, out_donor_gt, out_donor_id, out_singlet_prob, out_doublet_prob, out_summary, out_log,
              **kwargs):
    '''Run Vireo directing output to a Pypeliner temporary directory, then reformat and rename the output files.'''

    cmd = ['vireo',
           '--vartrixData={},{},{},{}'.format(alt_mat_file, ref_mat_file, barcodes, known_genotype_vcf),
           '--donorFile={}'.format(known_genotype_vcf),
           '--outDir={}'.format(vireo_tmp_dir),
           '--nDonor={}'.format(num_donors),
           '--randSeed={}'.format(random_seed)]

    pypeliner.commandline.execute(*cmd, **kwargs)

    donor_gt_file = os.path.join(vireo_tmp_dir, 'GT_donors.vireo.vcf.gz')
    donor_id_file = os.path.join(vireo_tmp_dir, 'donor_ids.tsv')
    singlet_prob_file = os.path.join(vireo_tmp_dir, 'prob_singlet.tsv.gz')
    doublet_prob_file = os.path.join(vireo_tmp_dir, 'prob_doublet.tsv.gz')
    summary_file = os.path.join(vireo_tmp_dir, 'summary.tsv')
    log_file = os.path.join(vireo_tmp_dir, '_log.txt')

    write_reformatted_vireo_files(known_id, sample_id, donor_gt_file, donor_id_file,
                                  singlet_prob_file, doublet_prob_file, summary_file, log_file,
                                  out_donor_gt, out_donor_id, out_singlet_prob, out_doublet_prob, out_summary, out_log)
