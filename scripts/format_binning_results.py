#!/usr/bin/python

'''
USAGE: python format_binning_results --tool TOOL --assembly PATH_TO_ASSEMBLY_GFA \
							  		 --results PATH_TO_BINNING_RESULTS_FILE \
							  		 --outdir PATH_TO_OUTPUT_FOLDER --outfile OUTPUT_FILENAME
'''

import argparse
import pandas as pd
import os
import io_utils

'''
Formatting functions:
- change the format of binning results into PlasEval input format,
- input: path to binning results file, bins (output) file, assembly graph object,
- output: None, changed format written directly to output file.
- (assembly graph required to access contig lengths)
'''
def format_pbf(RES_FILE, bins, assembly):
	'''
	PlasBin-flow:
	- binning results have one plasmid bin per line,
	- plasmid bin id in first column, list of contigs in last column of a TSV file,
	- list of contigs separated by commas,
	- each contig is of the form ctg_id:mul (type str:float)
	'''
	lines = open(RES_FILE, "r").read().split("\n")
	bins_file = open(bins, "w")
	bins_file.write("plasmid\tcontig\tcontig_len\n")
	for line in lines:
		if line and line[0] != '#':
			line = line.split("\t")
			line = [x for x in line if x]
			pls, ctgs = line[0], line[-1].split(',')
			for ctg in ctgs:
				ctg_id, mul = ctg.split(':')[0], float(ctg.split(':')[1])
				ctg_len = assembly.ctg_dict[ctg_id].length
				copy_num = int(mul)
				if copy_num >= 1:
					bins_file.write('PBF_'+pls + '\t' + ctg_id + '\t' + str(ctg_len) + '\n')
			bins_file.write("\n")

def format_hy(RES_FILE, bins, assembly):
	'''
	HyAsP:
	- binning results have one plasmid bin per line,
	- plasmid bin id and contigs separated by semicolon,
	- list of contigs separated by commas,
	- each contig is of the form ctg_id+ or ctg_id-.
	'''
	lines = open(RES_FILE, "r").read().split("\n")
	bins_file = open(bins, "w")
	bins_file.write("plasmid\tcontig\tcontig_len\n")
	for line in lines:
		if line and line[0] != '#':
			pls, ctgs = line.split(';')[0], line.split(';')[1].split(',')
			curr_ctg_set = set()
			for ctg in ctgs:
				ctg_id = ctg[:-1]	
				ctg_len = assembly.ctg_dict[ctg_id].length
				if ctg_id not in curr_ctg_set:
					bins_file.write('HY_'+pls + '\t' + ctg_id + '\t' + str(ctg_len) + '\n')
			bins_file.write("\n")

def format_mob(RES_FILE, bins, assembly):
	'''
	MOB-recon:
	- binning results in TSV file (one contig per line),
	- only lines with molecule_type 'plasmid' relevant
	- plasmid bin id is obtained from the primary_cluster_id column
	- ctg_id is obtained from the contig_id column 
	'''
	bins_file = open(bins, "w")
	bins_file.write("plasmid\tcontig\tcontig_len\n")
	ctg_df = pd.read_csv(RES_FILE, sep='\t')
	pls_df = ctg_df[ctg_df['molecule_type'] == 'plasmid']
	for index, row in pls_df.iterrows():
		pls, ctg_id = row['primary_cluster_id'], str(row['contig_id'])
		ctg_len = assembly.ctg_dict[ctg_id].length
		bins_file.write('MOB_'+pls + '\t' + ctg_id + '\t' + str(ctg_len) + '\n')

def format_gp(RES_FILE, bins, assembly):
	'''
	gplas2:
	- binning results in a space separated file (one contig per line),
	- ctg_id in first column, plasmid bin id in last column,
	- some contigs might be unbinned, these are treated as individual plasmids,
	'''
	lines = open(RES_FILE, "r").read().split("\n")
	bins_file = open(bins, "w")
	bins_file.write("plasmid\tcontig\tcontig_len\n")
	unbinned_count = 0
	for line in lines[1:]:
		ctg_id, pls = line.split(' ')[0], line.split(' ')[-1]
		if pls == 'Unbinned':
			unbinned_count += 1
			pls = 'unb_'+str(unbinned_count)
		ctg_len = assembly.ctg_dict[ctg_id].length
		bins_file.write('GP_'+pls + '\t' + ctg_id + '\t' + str(ctg_len) + '\n')



argparser = argparse.ArgumentParser()
argparser.add_argument('--tool', help = 'pbf / hy / mob / gp')
argparser.add_argument('--assembly', help = 'path to assembly gfa file')
argparser.add_argument('--results', help = 'path to binning results file')
argparser.add_argument('--outdir', help = 'output directory')
argparser.add_argument('--outfile', help = 'name of output TSV file')

args = argparser.parse_args()

TOOL = args.tool
ASSEMBLY_FILE = args.assembly
RES_FILE = args.results
OUT_DIR = args.outdir
OUT_FILE = args.outfile

if not os.path.exists(OUT_DIR):
	os.makedirs(OUT_DIR)
bins = open(OUT_DIR + '/' + OUT_FILE, "w")

assembly = io_utils.Assembly(ASSEMBLY_FILE)

if TOOL == 'pbf':
	format_pbf(TOOL, RES_FILE, bins, assembly)
elif TOOL == 'hy':
	format_hy(TOOL, RES_FILE, bins, assembly)
elif TOOL == 'mob':
	format_mob(TOOL, RES_FILE, bins, assembly)
elif TOOL == 'gp':
	format_gp(TOOL, RES_FILE, bins, assembly)