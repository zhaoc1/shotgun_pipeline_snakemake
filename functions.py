<<<<<<< HEAD
import json
import re

def get_sample_ids_from_barcode(bc_file):
	with open(bc_file) as f:
		lines = f.read().splitlines()
	filenames = []
	for line in lines:
		filenames.append(line.split("\t")[0])
	return(filenames)
=======
import itertools
import glob
import pathlib
import random
import os
import collections
import json
import re


from decontam.decontamlib import sam
from decontam.decontamlib import utils
from decontam.decontamlib.fastq import FastqSplitter
from decontam.decontamlib.version import __version__

def get_sample_ids_from_filenames(input_dir, suffix):
    dir_path = pathlib.Path(input_dir)
    input_file_paths = [x for x in dir_path.iterdir() if x.is_file()]
    suffix_paths = [x for x in input_file_paths if x.name.endswith(suffix)]
    return [x.name.replace(suffix, "") for x in suffix_paths]

def get_mapped_reads(filename,pct,frac):
	mapped = set()
	for qname, is_read1, rname in sam.get_mapped_reads(filename, pct, frac):
		if rname is not None:
			mapped.add(qname)
	return mapped

def annotate(sam_file,R1_file,pct=0.5,frac=0.6):
	# pct: default percent identity; frac: default fraction of alignment length
	mapped = get_mapped_reads(sam_file, pct, frac)
	ids = utils.parse_read_ids(R1_file)
	return([(id, True if id in mapped else False) for id in ids])

def save_summary(f,data):
	result = {
		"program": "decontam",
		"version": __version__,
		"data": data,
	}
	with open(f, "w") as f:
		f.write(json.dumps(result))
>>>>>>> 52b77d8d95bf1bec29f585765d2e0e7460c2a552

def generate_stitch_summary_json(log_file, json_file):
	with open(log_file) as f:
		lines = f.read().splitlines()

	for line in lines:
		if line.startswith('Assembled reads ...................:'):
			assembled_num = re.findall('[0-9]+',line)[0]
		if line.startswith('Discarded reads ...................: '):
			discarded_num = re.findall('[0-9]+',line)[0]
		if line.startswith('Not assembled reads ...............: '):
			unassembled_num = re.findall('[0-9]+',line)[0]
	
	result = {
		"program": "stitch",
		"version": "0.0.1",
		"data": {'assembled_num':assembled_num,'discarded_num':discarded_num,"unassembled_num":unassembled_num}
	}
	
	with open(json_file, "w") as f:
		f.write(json.dumps(result))