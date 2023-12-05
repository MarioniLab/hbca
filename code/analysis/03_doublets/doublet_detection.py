#######################
## Doublet detection ##
#######################


import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns
import scipy.stats
import anndata
import os



def bh(pvalues):
	# Computes the Benjamini-Hochberg FDR correction
	n = int(pvalues.shape[0])
	new_pvalues = np.empty(n)
	values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
	values.sort()
	values.reverse()
	new_values = []
	for i, vals in enumerate(values):
		rank = n - i
		pvalue, index = vals
		new_values.append((n/rank) * pvalue)
	for i in range(0, int(n)-1):
		if new_values[i] < new_values[i+1]:
			new_values[i+1] = new_values[i]
	for i, vals in enumerate(values):
		pvalue, index = vals
		new_pvalues[index] = new_values[i]
	return new_pvalues
    

import re
import fnmatch

files_processed = fnmatch.filter(os.listdir('data/scrublet-scores'), "*.csv")
sampleIDs_processed = [re.sub('\.csv$', '', sample) for sample in files_processed]
sampleIDs = os.listdir('/nfs/research/marioni/asteif/projects/hbca/cellranger/count/2021-07-10')

sample_list = [x for x in sampleIDs if x not in sampleIDs_processed]


#there's loads of clustering going on, so set verbosity low unless you enjoy walls of text
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)

scorenames = ['scrublet_score','scrublet_cluster_score','bh_pval']
os.makedirs('data/scrublet-scores', exist_ok=True)
#loop over the subfolders of the rawdata folder
for sample in [i for i in sample_list if os.path.isfile('/nfs/research/marioni/asteif/projects/hbca/cellranger/count/2021-07-10/'+i+'/outs/filtered_feature_bc_matrix.h5')]:
    #import data
    adata_sample = sc.read_10x_h5('/nfs/research/marioni/asteif/projects/hbca/cellranger/count/2021-07-10/'+sample+'/outs/filtered_feature_bc_matrix.h5')
    adata_sample.var_names_make_unique()
    #rename cells to SAMPLE_BARCODE, cleaving the trailing -1
    #adata_sample.obs_names = [sample+'_'+i.split('-')[0] for i in adata_sample.obs_names]
    adata_sample.obs_names = [i.split('-')[0]+'-'+sample for i in adata_sample.obs_names]
    #set up and run Scrublet
    scrub = scr.Scrublet(adata_sample.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata_sample.obs['scrublet_score'] = doublet_scores
    #overcluster prep. run turbo basic scanpy pipeline
    sc.pp.filter_genes(adata_sample, min_cells=3)
    sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_sample)
    sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
    sc.pp.scale(adata_sample, max_value=10)
    sc.tl.pca(adata_sample, svd_solver='arpack')
    sc.pp.neighbors(adata_sample)
    #eoverclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(adata_sample)
    for clus in np.unique(adata_sample.obs['leiden']):
        sc.tl.leiden(adata_sample, restrict_to=('leiden',[clus]))
        adata_sample.obs['leiden'] = adata_sample.obs['leiden_R']
    #compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(adata_sample.obs['leiden']):
        adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_cluster_score'] = \
            np.median(adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_score'])
    #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(adata_sample.obs['scrublet_cluster_score'])
    mask = adata_sample.obs['scrublet_cluster_score']>med
    mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
    #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(adata_sample.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    adata_sample.obs['bh_pval'] = bh(pvals)
    #create results data frame for single sample and copy stuff over from the adata object
    scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
    for meta in scorenames:
        scrublet_sample[meta] = adata_sample.obs[meta]
    #write out complete sample scores
    scrublet_sample.to_csv('data/scrublet-scores/'+sample+'.csv')



##########################
## Merge doublet scores ##
##########################

import os
import glob
import pandas as pd

sample_filenames = [i for i in glob.glob('data/scrublet-scores/*.{}'.format('csv'))]

scrublet_merged = pd.concat([pd.read_csv(f) for f in sample_filenames])
scrublet_merged['is_doublet'] = scrublet_merged['bh_pval'] < 0.1
scrublet_merged.rename(columns={scrublet_merged.columns[0]:'cellID'}, inplace=True)

#export to csv
scrublet_merged.to_csv( "data/scrublet-scores/scrublet_scores.csv.gz", index=False, encoding='utf-8-sig', compression='gzip')


# run on 2022-05-16
# bsub -M 60GB -n 10 -N -q "research" -Is "singularity shell --home /nfs/research/marioni/kunz:/home /hps/software/users/marioni/kunz/scAnalysis_v0.1.6.sif"
# sc.logging.print_header()
# scanpy==1.8.1 anndata==0.7.6 umap==0.5.1 numpy==1.20.2 scipy==1.6.2 pandas==1.2.5 scikit-learn==0.24.2 statsmodels==0.12.2 python-igraph==0.9.1 pynndescent==0.5.4
# pip3 freeze | grep scrublet
# scrublet==0.2.3