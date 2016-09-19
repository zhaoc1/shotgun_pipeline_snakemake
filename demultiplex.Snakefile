import itertools
import glob
import pathlib
import random

from snakemake.utils import update_config

def get_sample_ids_from_barcode(bc_file):
	with open(bc_file) as f:
		lines = f.read().splitlines()
	filenames = []
	for line in lines:
		filenames.append(line.split("\t")[0])
	return(filenames)

default_config = {
	"work_dir": "/Users/zhaoc1/pipeline1/test",
	"LANE_NUM": "001"
}
update_config(default_config, config)
config=default_config

# SampleIDs 
BC = config["work_dir"] + "/barcodes.txt"
SAMPLE_IDS = get_sample_ids_from_barcode(BC)
print(SAMPLE_IDS)

# the list of fastq files you want
TARGET_FPS = expand(
    "{out}/PCMP_{sample_id}_{read}.fastq",
    out=config["work_dir"]+"/dnabc_results", sample_id=SAMPLE_IDS, 
    read=["R1","R2"])
print("\n".join(TARGET_FPS))

rule all:
	input: TARGET_FPS

rule demultiplex:
	input:
		read1 = config["work_dir"] + "/Undetermined_S0_L" + config["LANE_NUM"] + "_R1_001.fastq",
		read2 = config["work_dir"] + "/Undetermined_S0_L" + config["LANE_NUM"] + "_R2_001.fastq"
	output:
	    TARGET_FPS
	params:
		dnabc_output_dir = config["work_dir"] + "/dnabc_results",
		dnabc_summary = config["work_dir"] + "/summary/summary-dnabc.json"
	threads: 8
	shell:
		"""
		mkdir -p {config[work_dir]}/summary
		dnabc.py --forward-reads {input.read1} --reverse-reads {input.read2} \
		--barcode-file {BC} --output-dir {params.dnabc_output_dir} \
		--summary-file {params.dnabc_summary}
		"""