#!/usr/bin/env bash

WORKFLOW=/nfs/research/marioni/asteif/projects/hbca/code/single_cell_gene_expression/genotyping/vireo_workflow.py
CONFIG=/nfs/research/marioni/asteif/projects/hbca/genotyping/vireo/2021-07-10/config_vireo.yaml

python $WORKFLOW \
	--config $CONFIG \
	--submit lsf \
	--maxjobs 200 \
	--nativespec ' -n {ncpus} -M {mem} -W 48:00 -R "rusage[mem={mem}]"' \
	--nocleanup

