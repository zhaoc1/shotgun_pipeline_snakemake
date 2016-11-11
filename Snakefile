####
# shotgun metagenomic pipeline
# author: Chunyu Zhao
# time: 2016-11-09
####

import glob
import pathlib
import random
import configparser
import os
import sys
import shutil
import yaml
import subprocess
from functions import * 

"""
run on respublica

unset PYTHONHOME
unset PYTHONPATH
source activate shotgun-pipeline

snakemake -j 16 --cluster-config configs/cluster.yaml -c "qsub -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"
snakemake -j 2 --cluster-config configs/cluster.test.yaml -c "qsub -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"
"""

configfile: "configs/localconfig.yaml"
workdir: config["data_dir"]

#shell.prefix("source ~/.bashrc; ")

include:
	"rules/targets.rules"
include:
	"rules/qc.rules"
include:
	"rules/decontam.rules"
include:
	"rules/kraken.rules"
include:
	"rules/reports.rules"

rule all:
	input: 
		TARGET_FPS

onsuccess:
	print("Workflow finished, no error")
	shell("mail -s 'workflow finished' " + config['admins']+" <{log}")
onerror:
	print("An error occurred")
	shell("mail -s 'an error occurred' " + config['admins']+" < {log}")
