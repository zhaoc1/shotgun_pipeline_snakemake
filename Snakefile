<<<<<<< HEAD
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
=======
import itertools
import glob
import pathlib
import random

from snakemake.utils import update_config

import os
from functions import *

def get_sample_ids_from_barcode(bc_file):
	with open(bc_file) as f:
		lines = f.read().splitlines()
	filenames = []
	for line in lines:
		filenames.append(line.split("\t")[0])
	return(filenames)

default_config = {
	"work_dir": "/Users/zhaoc1/Work/github/shotgun_pipeline_snakemake_v2/test",
	"LANE_NUM": "001"
}
update_config(default_config, config)
config=default_config


## SampleIDs 
BC = config["work_dir"] + "/barcodes.txt"
SAMPLE_IDS = get_sample_ids_from_barcode(BC)
print(SAMPLE_IDS)

## target fastq files for the demultiplex
TARGET_FPS0 = expand(
    "{out}/PCMP_{sample_id}_{read}.fastq",
    out=config["work_dir"]+"/dnabc_results", sample_id=SAMPLE_IDS, 
    read=["R1","R2"])
#print("\n".join(TARGET_FPS0))
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552

## target fastq files we want
TARGET_FPS1 = expand(
	"{out}/PCMP_{sample_id}.assembled.fastq",
<<<<<<< HEAD
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
=======
	out=config["work_dir"]+"/decontam_host_assembled_results", 
	sample_id=SAMPLE_IDS)
#print("\n".join(TARGET_FPS1))

TARGET_FPS2 = expand(
	"{out}/PCMP_{sample_id}.{type}.fastq",
	out=config["work_dir"]+"/decontam_host_unassembled_results", 
	sample_id=SAMPLE_IDS, type=["unassembled.forward","unassembled.reverse"])
#print("\n".join(TARGET_FPS2))

TARGET_FPS = TARGET_FPS1 + TARGET_FPS2

rule all:
	input: TARGET_FPS

## demultiplex the fastq files
rule demultiplex:
	input:
		read1 = config["work_dir"] + "/Undetermined_S0_L" + config["LANE_NUM"] + "_R1_001.fastq",
		read2 = config["work_dir"] + "/Undetermined_S0_L" + config["LANE_NUM"] + "_R2_001.fastq"
	output:
	    TARGET_FPS0
	params:
		dnabc_output_dir = config["work_dir"] + "/dnabc_results",
		dnabc_summary = config["work_dir"] + "/summary/summary-dnabc.json"
	threads: 1
	shell:
		"""
		mkdir -p {config[work_dir]}/summary
		dnabc.py --forward-reads {input.read1} --reverse-reads {input.read2} \
		--barcode-file {BC} --output-dir {params.dnabc_output_dir} \
		--summary-file {params.dnabc_summary}
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552
		"""

## quality control
rule illqc:
	input:
<<<<<<< HEAD
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
=======
		read1=config["work_dir"] + "/dnabc_results/{sample_id}_R1.fastq",
		read2=config["work_dir"] + "/dnabc_results/{sample_id}_R2.fastq"
	output:
		config["work_dir"] + "/illqc_results/{sample_id}_R1.fastq",
		config["work_dir"] + "/illqc_results/{sample_id}_R2.fastq"
	params:
		illqc_output_dir = config["work_dir"] + "/illqc_results",
		illqc_report_dir = config["work_dir"] + "/illqc_reports",
		illqc_summary = config["work_dir"]+"/summary/summary-illqc_{sample_id}.json"
	shell:
		"illqc.py --forward-reads {input.read1} "
		"--reverse-reads {input.read2} --output-dir {params.illqc_output_dir} "
		"--qc-output-dir {params.illqc_report_dir} --summary-file {params.illqc_summary}"
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552

## stitch paired-end reads
rule stitch_reads:
	input:
<<<<<<< HEAD
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
=======
		read1=config["work_dir"] + "/illqc_results/{sample_id}_R1.fastq",
		read2=config["work_dir"] + "/illqc_results/{sample_id}_R2.fastq"
	output:
		config["work_dir"] + "/stitch_results/{sample_id}.assembled.fastq",
		config["work_dir"] + "/stitch_results/{sample_id}.discarded.fastq",
		config["work_dir"] + "/stitch_results/{sample_id}.unassembled.forward.fastq",
		config["work_dir"] + "/stitch_results/{sample_id}.unassembled.reverse.fastq"
	params:
		prefix = config["work_dir"] + "/stitch_results/{sample_id}",
		log=config["work_dir"] + "/stitch_results/{sample_id}.log",
		stitch_summary = config["work_dir"]+"/summary/summary-stitch_{sample_id}.json"
	threads: 1
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552
	run:
		shell(
			"""
			pear --forward-fastq {input.read1} --reverse-fastq {input.read2} \
<<<<<<< HEAD
			--output {params.prefix} --threads {threads} --memory 1G \
			--p-value 0.01 --min-overlap 10 --min-assembly-length 50 \
			--quality-threshold 0 --test-method 1 --score-method 2 \
			--phred-base 33  --keep-original 1> {params.log}
=======
			--output {params.prefix} --threads {threads} \
			--p-value 0.01 --min-overlap 10 --min-assembly-length 50 \
			--quality-threshold 0 --test-method 1 --score-method 2 \
			--phred-base 33 --memory 200M --keep-original 1> {params.log}
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552
			""")
		generate_stitch_summary_json(params.log,params.stitch_summary)

#＃ phix decontamination for unassembled forward and reverse reads pair
rule decontam_phix_unassembled:
	input:
<<<<<<< HEAD
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
=======
		read1=config["work_dir"] + "/stitch_results/{sample_id}.unassembled.forward.fastq",
		read2=config["work_dir"] + "/stitch_results/{sample_id}.unassembled.reverse.fastq"
	output:
		config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.forward.fastq",
		config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.reverse.fastq"
	params:
		decontam_phix_output_dir_unassembled = config["work_dir"] + 
						"/decontam_phix_unassembled_results",
		decontam_phix_summary_file_unassembled = config["work_dir"] + 
						"/summary/summary-decontam_phix_{sample_id}.unassembled.json"
	shell:
		"decontaminate.py --forward-reads {input.read1} --reverse-reads {input.read2} "
		"--output-dir {params.decontam_phix_output_dir_unassembled} "
		"--summary-file {params.decontam_phix_summary_file_unassembled} --organism phix"

## phix decontamination for assembled reads
rule decontam_phix_assembled:
	input:
		index = "/Users/zhaoc1/biodata/phix174.fasta",
		read = config["work_dir"] + "/stitch_results/{sample_id}.assembled.fastq"
	output:
		config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.fastq"
	params:
		decontam_phix_output_dir_assembled = config["work_dir"] + "/decontam_phix_assembled_results",
		decontam_phix_summary_file_assembled = config["work_dir"] + "/summary/summary-decontam_phix_{sample_id}.assembled.json",
		sam_assembled = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.sam",
		err_assembled = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.err"
	threads: 1
	run:
		## check human genome index exists: if not then make index
		if not os.path.exists(input.index +".amb"):
			shell("bwa index {input.index}")
		## align the assembled reads to the human genome
		shell("bwa mem -M -t {threads} {input.index} {input.read} 1> {params.sam_assembled} 2> {params.err_assembled}")
		## decontaminate: sam -> two fastqs + one json
		annotations = annotate(params.sam_assembled,input.read)
		with FastqSplitter(input.read, params.decontam_phix_output_dir_assembled) as s:
			s.partition(annotations, "phix")
		summary_data = dict(collections.Counter(a for _, a in annotations))
		save_summary(params.decontam_phix_summary_file_assembled, summary_data) #no config info though

#＃ phix decontamination for unassembled forward and reverse reads pair
rule decontam_host_unassembled:
	input:
		read1 = config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.forward.fastq",
		read2 = config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.reverse.fastq"
	output:
		config["work_dir"] + "/decontam_host_unassembled_results/{sample_id}.unassembled.forward.fastq",
		#config["work_dir"] + "/decontam_host_unassembled_results/{sample_id}.unassembled.forward_human.fastq",
		config["work_dir"] + "/decontam_host_unassembled_results/{sample_id}.unassembled.reverse.fastq"
		#config["work_dir"] + "/decontam_host_unassembled_results/{sample_id}.unassembled.reverse_human.fastq"
	params:
		decontam_host_output_dir_unassembled = config["work_dir"] + "/decontam_host_unassembled_results",
		decontam_host_summary_file_unassembled = config["work_dir"] + "/summary/summary-decontam_host_{sample_id}.unassembled.json"
	shell:
		"decontaminate.py --forward-reads {input.read1} --reverse-reads {input.read2} "
		"--output-dir {params.decontam_host_output_dir_unassembled} "
		"--summary-file {params.decontam_host_summary_file_unassembled} --organism human"
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552

#＃ human decontamination for assembled reads
rule decontam_host_assembled:
	input:
<<<<<<< HEAD
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
=======
		index = "/Users/zhaoc1/biodata/human_GRch38.fasta",
		read = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.fastq"
	output:
		assembled = config["work_dir"] + "/decontam_host_assembled_results/{sample_id}.assembled.fastq"
		#assembled_human = config["work_dir"] + "/decontam_host_assembled_results/{sample_id}.assembled_human.fastq"
	params:
		decontam_host_output_dir_assembled = config["work_dir"] + "/decontam_host_assembled_results",
		decontam_host_summary_file_assembled = config["work_dir"] + "/summary/summary-decontam_host_{sample_id}.assembled.json",
		sam_assembled = temp(config["work_dir"] + "/decontam_host_assembled_results/{sample_id}.assembled.sam"),
		err_assembled = temp(config["work_dir"] + "/decontam_host_assembled_results/{sample_id}.assembled.err")
	threads: 1
	run:
		## check human genome index exists: if not then make index
		if not os.path.exists(input.index +".amb"):
			shell("bwa index {input.index}")
		## align the assembled reads to the human genome
		shell("bwa mem -M -t {threads} {input.index} {input.read} 1> {params.sam_assembled} 2> {params.err_assembled}")
		## decontaminate: sam -> two fastqs + one json
		annotations = annotate(params.sam_assembled,input.read)
		with FastqSplitter(input.read, params.decontam_host_output_dir_assembled) as s:
			s.partition(annotations, "human")
		summary_data = dict(collections.Counter(a for _, a in annotations))
		save_summary(params.decontam_host_summary_file_assembled, summary_data) #no config info though

>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552
