import os

from snakemake.utils import update_config

from functions import *


default_config = {
	"work_dir": "/Users/zhaoc1/pipeline1/test"
}
update_config(default_config, config)
config=default_config


SAMPLE_IDS = get_sample_ids_from_filenames(
	config["work_dir"]+"/dnabc_results", "_R1.fastq")
print(SAMPLE_IDS)

TARGET_FPS1 = expand(
	"{out}/{sample_id}.assembled.fastq",
	out=config["work_dir"]+"/decontam_host_assembled_results", 
	sample_id=SAMPLE_IDS)
print("\n".join(TARGET_FPS1))

TARGET_FPS2 = expand(
	"{out}/{sample_id}.{type}.fastq",
	out=config["work_dir"]+"/decontam_host_unassembled_results", 
	sample_id=SAMPLE_IDS, type=["unassembled.forward","unassembled.reverse"])#,"unassembled.forward_human","unassembled.reverse_human"
print("\n".join(TARGET_FPS2))

TARGET_FPS = TARGET_FPS1 + TARGET_FPS2

rule all:
	input: TARGET_FPS

## quality control
rule illqc:
	input:
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

## stitch paired-end reads
#todo: need to generate json summary
rule stitch_reads:
	input:
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
	threads: 4
	run:
		shell(
			"""
			pear --forward-fastq {input.read1} --reverse-fastq {input.read2} \
			--output {params.prefix} --threads {threads} \
			--p-value 0.01 --min-overlap 10 --min-assembly-length 50 \
			--quality-threshold 0 --test-method 1 --score-method 2 \
			--phred-base 33 --memory 200M --keep-original 1> {params.log}
			""")
		generate_stitch_summary_json(params.log,params.stitch_summary)

#＃ phix decontamination for unassembled forward and reverse reads pair
rule decontam_phix_unassembled:
	input:
		read1=config["work_dir"] + "/stitch_results/{sample_id}.unassembled.forward.fastq",
		read2=config["work_dir"] + "/stitch_results/{sample_id}.unassembled.reverse.fastq"
	output:
		config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.forward.fastq",
		#config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.forward_phix.fastq",
		config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.reverse.fastq"
		#config["work_dir"] + "/decontam_phix_unassembled_results/{sample_id}.unassembled.reverse_phix.fastq"
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
		#assembled_phix = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled_phix.fastq"#can be empty and not created..
	params:
		decontam_phix_output_dir_assembled = config["work_dir"] + "/decontam_phix_assembled_results",
		decontam_phix_summary_file_assembled = config["work_dir"] + "/summary/summary-decontam_phix_{sample_id}.assembled.json",
		sam_assembled = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.sam",
		err_assembled = config["work_dir"] + "/decontam_phix_assembled_results/{sample_id}.assembled.err"
	threads: 8
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

#＃ human decontamination for assembled reads
rule decontam_host_assembled:
	input:
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
	threads: 8
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

