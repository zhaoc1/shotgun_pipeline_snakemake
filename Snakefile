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

snakemake -j 2 --cluster-config configs/cluster.yaml -c "qsub -r n -V -l h_vmem={cluster.g_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"

"""

configfile: "configs/localconfig.yaml"

#shell.prefix("source ~/.bashrc; ")


## SampleIDs 
barcodes = config["data_dir"] + "/barcodes.txt"
SAMPLE_IDS = get_sample_ids_from_barcode(barcodes)

## target fastq for demultiplex
TARGET_FPS0 = expand("{out}/PCMP_{sample_id}_{read}.fastq",
    out=config["data_dir"]+"/dnabc_results", sample_id=SAMPLE_IDS, 
    read=["R1","R2"])

## target fastq files we want
TARGET_FPS1 = expand(
	"{out}/PCMP_{sample_id}.assembled.fastq",
	out=config["data_dir"]+"/decontam_host_assembled_results", 
	sample_id=SAMPLE_IDS)

TARGET_FPS2 = expand(
	"{out}/PCMP_{sample_id}.{type}.fastq",
	out=config["data_dir"]+"/decontam_host_unassembled_results", 
	sample_id=SAMPLE_IDS, type=["unassembled.forward","unassembled.reverse"])

TARGET_FPS = TARGET_FPS1 + TARGET_FPS2


rule all:
	input: TARGET_FPS

## demultiplex 
rule demultiplex:
	input:
		read1 = config["data_dir"] + "/Undetermined_S0_L" + config["lane_num"] + "_R1_001.fastq",
		read2 = config["data_dir"] + "/Undetermined_S0_L" + config["lane_num"] + "_R2_001.fastq"
	output:
	    TARGET_FPS0
	params:
		dnabc_output = config["data_dir"] + config["output"]["dnabc"],
		dnabc_summary = config["data_dir"] + "/summary/dnabc.json"
	log: config["data_dir"] + "/log/dnabc.log"
	threads: 8
	shell:
		"""
		mkdir -p {config[data_dir]}/summary
		mkdir -p {config[data_dir]}/log
		dnabc.py --forward-reads {input.read1} --reverse-reads {input.read2} \
		--barcode-file {barcodes} --output-dir {params.dnabc_output} \
		--summary-file {params.dnabc_summary} 2> {log}
		"""

## quality control
rule illqc:
	input:
		read1 = config["data_dir"] + config["output"]["dnabc"] + "/{sample_id}_R1.fastq",
		read2 = config["data_dir"] + config["output"]["dnabc"] + "/{sample_id}_R2.fastq"
	output:
		config["data_dir"] + config["output"]["illqc"] + "/{sample_id}_R1.fastq",
		config["data_dir"] + config["output"]["illqc"] + "/{sample_id}_R2.fastq"
	params:
		illqc_output = config["data_dir"] + config["output"]["illqc"],
		illqc_report = config["data_dir"] + config["report"]["illqc"],
		illqc_summary = config["data_dir"] + "/summary/{sample_id}.illqc.json"
	shell:
		"illqc.py --forward-reads {input.read1} --reverse-reads {input.read2} "
		"--output-dir {params.illqc_output} --qc-output-dir {params.illqc_report} "
		"--summary-file {params.illqc_summary}"

## stitch paired-end reads
rule stitch_reads:
	input:
		read1 = config["data_dir"] + config["output"]["illqc"] + "/{sample_id}_R1.fastq",
		read2 = config["data_dir"] + config["output"]["illqc"] + "/{sample_id}_R2.fastq"
	output:
		config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.assembled.fastq",
		config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.discarded.fastq",
		config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.unassembled.forward.fastq",
		config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.unassembled.reverse.fastq"
	params:
		prefix = config["data_dir"] + config["output"]["stitch"] + "/{sample_id}",
		stitch_summary = config["data_dir"] + "/summary/{sample_id}.stitch.json",
		log = config["data_dir"] + "/log/{sample_id}.stitch.log"
	threads: 8 
	run:
		shell(
			"""
			pear --forward-fastq {input.read1} --reverse-fastq {input.read2} \
			--output {params.prefix} --threads {threads} --memory 1G \
			--p-value 0.01 --min-overlap 10 --min-assembly-length 50 \
			--quality-threshold 0 --test-method 1 --score-method 2 \
			--phred-base 33  --keep-original 1> {params.log}
			""")
		generate_stitch_summary_json(params.log,params.stitch_summary)

#＃ phix decontamination for unassembled forward and reverse reads pair
rule decontam_phix_unassembled:
	input:
		read1 = config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.unassembled.forward.fastq",
		read2 = config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.unassembled.reverse.fastq"
	output:
		config["data_dir"] + config["output"]["decontam_phix_unassembled"] + "/{sample_id}.unassembled.forward.fastq",
		config["data_dir"] + config["output"]["decontam_phix_unassembled"] + "/{sample_id}.unassembled.reverse.fastq"
	params:
		decontam_phix_output_unassembled = config["data_dir"] + config["output"]["decontam_phix_unassembled"],
		decontam_phix_summary_unassembled = config["data_dir"] + "/summary/{sample_id}.unassembled.decontam.phix.json"
	shell:
		"decontaminate.py --forward-reads {input.read1} --reverse-reads {input.read2} "
		"--output-dir {params.decontam_phix_output_unassembled} "
		"--summary-file {params.decontam_phix_summary_unassembled} --organism phix"

#＃ human decontamination for unassembled forward and reverse reads pair
rule decontam_host_unassembled:
	input:
		read1 = config["data_dir"] + config["output"]["decontam_phix_unassembled"] + "/{sample_id}.unassembled.forward.fastq",
		read2 = config["data_dir"] + config["output"]["decontam_phix_unassembled"] + "/{sample_id}.unassembled.reverse.fastq"
	output:
		config["data_dir"] + config["output"]["decontam_host_unassembled"] + "/{sample_id}.unassembled.forward.fastq",
		config["data_dir"] + config["output"]["decontam_host_unassembled"] + "/{sample_id}.unassembled.reverse.fastq"
	params:
		decontam_host_output_unassembled = config["data_dir"] + config["output"]["decontam_host_unassembled"],
		decontam_host_summary_unassembled = config["data_dir"] + "/summary/{sample_id}.unassembled.decontam.host.json"
	shell:
		"decontaminate.py --forward-reads {input.read1} --reverse-reads {input.read2} "
		"--output-dir {params.decontam_host_output_unassembled} "
		"--summary-file {params.decontam_host_summary_unassembled} --organism human"

## phix decontamination for assembled reads
rule decontam_phix_assembled:
	input:
		read = config["data_dir"] + config["output"]["stitch"] + "/{sample_id}.assembled.fastq"
	output:
		config["data_dir"] + config["output"]["decontam_phix_assembled"] + "/{sample_id}.assembled.fastq"
	params:
		decontam_phix_output_assembled = config["data_dir"] + config["output"]["decontam_phix_assembled"],
		decontam_phix_summary_assembled = config["data_dir"] + "/summary/{sample_id}.assembled.decontam.phix.json"
	shell:
		"decontaminate.py --forward-reads {input.read} --output-dir {params.decontam_phix_output_assembled} "
		"--summary-file {params.decontam_phix_summary_assembled} --organism phix"

#＃ human decontamination for assembled reads
rule decontam_host_assembled:
	input:
		read = config["data_dir"] + config["output"]["decontam_phix_assembled"] + "/{sample_id}.assembled.fastq"
	output:
		config["data_dir"] + config["output"]["decontam_host_assembled"] + "/{sample_id}.assembled.fastq"
	params:
		decontam_host_output_assembled = config["data_dir"] + config["output"]["decontam_host_assembled"],
		decontam_host_summary_assembled = config["data_dir"] + "/summary/{sample_id}.assembled.decontam.host.json"
	shell:
		"decontaminate.py --forward-reads {input.read} --output-dir {params.decontam_host_output_assembled} "
		"--summary-file {params.decontam_host_summary_assembled} --organism human"

onsuccess:
	print("Workflow finished, no error")
	shell("mail -s 'workflow finished' "+config['admins']+" <{log}")
onerror:
	print("An error occurred")
    shell("mail -s 'an error occurred' "+config['admins']+" < {log}")