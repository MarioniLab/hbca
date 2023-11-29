#!/usr/bin/env bash

WORKFLOW=/nfs/research/marioni/areed/projects/hbca/code/single_cell_gene_expression/qualitycontrol/qualitycontrol_workflow.py
CONFIG=/nfs/research/marioni/areed/projects/hbca/qualitycontrol/2022-03-25/raw/config_qualitycontrol.yaml

python $WORKFLOW \
	--config_file $CONFIG \
	--submit lsf \
	--maxjobs 200 \
	--nativespec ' -n {ncpus} -M {mem} -W 48:00 -R "rusage[mem={mem}]"' \
	--nocleanup

