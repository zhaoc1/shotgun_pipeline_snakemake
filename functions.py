import json
import re

def get_sample_ids_from_barcode(bc_file):
	with open(bc_file) as f:
		lines = f.read().splitlines()
	filenames = []
	for line in lines:
		filenames.append(line.split("\t")[0])
	return(filenames)

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