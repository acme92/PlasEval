#!/usr/bin/python

'''
USAGE: python ground_truth_mapping.py --plasmids TRUE_PLASMIDS_FASTA --assembly ASSEMBLY_FASTA \
                                  	  --mapping BLAST_OUTPUT_TSV --outdir OUT_DIR --outfile OUT_FILENAME \
									  --min_pid PID_THR --min_cov COV_THR --max_len LEN_THR		
'''

from __future__ import division
import os
import argparse
import pandas as pd
from Bio import SeqIO

def read_blast_output(file):
	'''
	- Reads BLAST output file (outfmt 6) into table
	'''
	col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", \
			     "qstart", "qend", "sstart", "send", "evalue", "bitscore"]  
	return pd.read_csv(file, sep = '\t', names = col_names, dtype = str)

def num_covered_positions(intervals):
	'''
	- Computes number of positions covered by a list of potentially overlapping intervals
	- Input: List of intervals, each interval is a pair of integers (start, end)
	- Output: Length of short read contig mapped to reference sequences (int)
	'''
	intervals.sort(key = lambda x: x[0])  # intervals is now sorted by start position
	num_pos_covered = 0
	last_pos_covered = 0  # last (right-most) position of contig covered so far
	for start, end in intervals:
		if end <= last_pos_covered:
			pass  # contained in previous interval -> no new position covered
		else:
			num_pos_covered += end - max(last_pos_covered + 1, start) + 1
			last_pos_covered = end
	return num_pos_covered

def analyse_from_contigs(hits, contig_seqs, reference_seqs, thresholds, out):
	'''
	- Analyse the matches between the contigs of predicted plasmids and the reference sequences
	- Input:
		- hits: dataframe of blast output
		- contig_seqs: dictionary of short read contigs with contig lengths as values
		- reference_seqs: dictionary of reference sequences with sequence lengths as values
		- thresholds: list of three thresholds
			- pident threshold (pid_thr) (type: float)
			(percent identity of a particular match b/n a contig and a reference sequence)
			- coverage threshold (cov_thr) (type: float)
			(proportion of short read contig cumulatively covered by reference sequences)
			-length threshold (len_thr) (type: int)
			(maximum length of reference sequences to be considered a plasmid;
			 note that this threshold is specific to hybrid sequences)
		- out: path to output file
	- Output: None, Directly written to TSV output file.
		- output format: one contig per line
		- plasmid\tcontig\tcontig_len
	'''
	covered_sections = dict([(pred, []) for pred in contig_seqs])	
	covered_per_ref = {}		#Covered proportion per reference plasmid
	ref_ids = hits.sseqid.unique()	#Set of reference ids in the blast output
	pid_thr = thresholds[0]
	cov_thr = thresholds[1]
	len_thr = thresholds[2]
	for ref in sorted(ref_ids):
		covered_per_ref[ref] = {}
		for ctg in contig_seqs:
			covered_per_ref[ref][ctg] = []
			ctg_ref_hits = hits.loc[hits.qseqid == ctg].loc[hits.sseqid == ref]
			for index, row in ctg_ref_hits.iterrows():
				qstart, qend = int(row[6]), int(row[7])
				pident = float(row[2])/100
				interval = (qstart, qend) if qstart <= qend else (qend, qstart)
				if pident >= pid_thr:
					covered_sections[ctg].append(interval)
					covered_per_ref[ref][ctg].append(interval)
	for ref in covered_per_ref:
		for ctg in covered_per_ref[ref]:
			ctg_len = len(str(contig_seqs[ctg]))
			percent_mapping = num_covered_positions(covered_per_ref[ref][ctg])/ctg_len
			ref_len = reference_seqs[ref]
			if percent_mapping >= cov_thr:
				if ref_len <= len_thr:
					out.write('GT_'+ref + '\t' + ctg + '\t' + str(ctg_len) + "\n")
	   
argparser = argparse.ArgumentParser()
argparser.add_argument('--references', help = 'reference FASTA file')
argparser.add_argument('--assembly', help = 'assembly FASTA file')
argparser.add_argument('--mapping', help = 'blast output TSV file')
argparser.add_argument('--min_pid', const=0.95, type=float, help = 'minimum percent identity threshold (float)')
argparser.add_argument('--min_cov', const=0.8, type=float, help = 'minimum coverage proportion threshold (float)')
argparser.add_argument('--max_len', const=500000, type=int, help = 'maximum length threshold (int)')
argparser.add_argument('--outdir', help = 'output directory')
argparser.add_argument('--outfile', help = 'name of output TSV file')

args = argparser.parse_args()

REF_FASTA = args.references
CTG_FASTA = args.assembly
BLAST_OUT = args.mapping
OUT_DIR = args.outdir
OUT_FILENAME = args.outfile
THRESHOLDS = [args.min_pid, args.min_cov, args.max_len]  

mapping = read_blast_output(BLAST_OUT)

ctg_seqs = dict()
with open(CTG_FASTA, "r") as f:
	for record in SeqIO.parse(f, "fasta"):
		ctg_seqs[record.id] = record.seq

ref_seqs = dict()
with open(REF_FASTA, "r") as f:
	for record in SeqIO.parse(f, "fasta"):
		ref_seqs[record.id] = len(record.seq)			

os.system('mkdir -p '+OUT_DIR)
OUT_FILE = OUT_DIR + '/' + OUT_FILENAME
with open(OUT_FILE, "w") as out:
	analyse_from_contigs(mapping, ctg_seqs, ref_seqs, THRESHOLDS, out)