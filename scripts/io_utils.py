#!/usr/bin/python

from __future__ import division
import scipy.special as sc
import scipy.integrate as integrate
import math



#Initiating, computing, storing and returning contig details.
class Contig():
	def __init__(self, seq, length, rd):
		'''
		Initiates Contig attributes
		sequence -> seq (str), length -> len (int), 
		read depth -> rd (float), GC content -> GC (float)
		'''
		self.seq = seq
		self.length = length
		self.rd = rd
		self.GC = (seq.count('G') + seq.count('C'))/length if length > 0 else 0

	def set_seq(self, seq):
		'''
		Stores the contig sequence (str)
		'''
		self.seq = seq

	def set_len(self, length):
		'''
		Stores the length of the contig sequence (integer)
		'''
		self.length = int(length)

	def set_rd(self, rd):
		'''
		Stores the read depth of the contig sequence (float)
		'''
		self.rd = float(rd)	

	def get_seq(self):
		'''
		Returns the sequence of the contig
		'''
		return self.seq

	def get_len(self):
		'''
		Returns the length of the contig
		'''		
		return self.length

	def get_rd(self):
		'''
		Returns the read depth of the contig
		'''		
		return self.rd

	def get_GC(self):
		'''
		Returns the GC content of the contig
		'''		
		return self.GC
	


#Initiating, computing, storing and returning assembly graph details.
class Assembly:
	def __init__(self, assembly_file):
		'''
		Initiates Assembly graph attributes
		Contigs -> ctg_dict (dictionary of Contigs), 
		Edges -> edge_lst (list of edges)		
		'''
		self.ctg_dict = {}
		self.edge_lst = []
		self.parse_gfa(assembly_file)

	def init_ctg(self, ctg_id, seq, length, rd):
		'''
		Initiates an instance of class Contig
		'''
		self.ctg_dict[ctg_id] = Contig(seq, length, rd)

	def add_edge(self, edge):
		'''
		Adds edge to the edge list
		An edge is a set of a pair of extremities 
		An extremity is a 2-tuple of (ctg_id, ctg_ext)
		'''
		self.edge_lst.append(edge)

	def parse_gfa(self, assembly_file):
		'''
		Reads and processes the assembly graph (gfa) file
		'''
		def parse_vertex(line):
			tmp = line.split('\t')
			# Contig id and sequence
			ctg_id, ctg_seq = tmp[1], tmp[2]
			# Tags
			tags_dict = {x.split(':')[0]:x.split(':')[2] for x in tmp[3:]}
			tags_keys = tags_dict.keys()
			# Reading contig length
			ctg_len = int(tags_dict['LN']) if 'LN' in tags_keys else len(ctg_seq)
			# Reading read depth 
			if 'dp' in tags_keys: ctg_rd = float(tags_dict['dp'])
			elif 'KC' in tags_keys: ctg_rd = int(tags_dict['KC'])/ctg_len
			else: ctg_rd = None
			# Recording values 
			self.init_ctg(ctg_id, ctg_seq, ctg_len, ctg_rd)		
		def parse_edge(line):
			[ctg1, ori1, ctg2, ori2] = line.split('\t')[1:-1]
			ext1 = 'h' if ori1 == '+' else 't'
			ext2 = 't' if ori2 == '+' else 'h'
			self.add_edge(((ctg1, ext1),(ctg2, ext2)))
		
		with open(assembly_file, 'r') as ag:
			line = next(ag)
			while line:
				if line[0].startswith('S'): parse_vertex(line)
				elif line.startswith('L'): parse_edge(line)
				line = next(ag, None)

	def get_ctgs(self):
		'''
		Returns list of contigs
		'''
		return self.ctg_dict.keys()
	
	def get_edges(self):
		'''
		Returns list of edges
		'''
		return self.edge_lst
	
	def get_nctg(self):
		'''
		Returns number of contigs in the graph
		'''
		return len(self.ctg_dict.keys())
	
	def get_nedges(self):
		'''
		Returns number of edges
		'''
		return len(self.edge_lst)

	def get_degree(self, ctg_id):
		'''
		Returns number of edges incident on a contig
		'''
		d = 0
		for edge in self.edge_lst:
			if edge[0][0] == ctg_id or edge[1][0] == ctg_id:
				d += 1
		return d	



#Initiating, computing and storing input specific to PlasBin-Flow.
class PBF_input():
	def __init__(self, assembly_file, gcint_file, mapping_file, seed_gd_thr, seed_len_thr):
		'''
		Initiates attributes specific to PlasBin-Flow
		'''
		self.assembly = Assembly(assembly_file)
		self.set_GC_intervals(gcint_file)       #List of GC interval endpoints: list of floats
		self.set_GC_probs()                 #Dictionary of GC interval probabilities: Key: Contig and Value: List of floats
		self.set_gd(mapping_file)
		self.set_seeds(seed_gd_thr, seed_len_thr)

	def set_GC_intervals(self, gcint_file):
		'''
		Defining the GC intervals
		'''
		self.GC_intervals = sorted([ 
			float(x) 
			for x in open(gcint_file, "r").read().split("\n") 
			if x and x[0] != '#'
		])
		
	#Following functions are required to compute the GC probabilities
	def set_GC_probs(self, m = 10):
		def gprob2(n,g,p,m):
			'''
			Compute probability of observing g GC nucleotides in a contig
			of length n within a molecule of GC content p using pseducount m. 
			Done via logarithm to avoid overflow.
			'''
			alpha = m*p
			beta = m*(1-p)
			combln = sc.gammaln(n + 1) - sc.gammaln(g + 1) - sc.gammaln(n - g + 1)
			resultln = combln + sc.betaln(g + alpha, n - g + beta) - sc.betaln(alpha, beta)
			return math.exp(resultln)
		
		#Computes the GC probabilities for all contigs and stores them in a dictionary
		self.GC_probs = {}
		for ctg in self.assembly.get_ctgs():
			# Contig length (n) and GC content proportion (gc) and content (n_gc)
			n = self.assembly.ctg_dict[ctg].get_len()
			gc = self.assembly.ctg_dict[ctg].get_GC()
			n_gc = n * gc
			# GC intervals
			GC_ints = self.GC_intervals
			# Computing probability for each GC interval
			gp_array, total = [], 0		
			for i in range(0, len(GC_ints)-1):
				gp_aux = integrate.quad(lambda x: gprob2(n,n_gc,x,m), GC_ints[i], GC_ints[i+1])
				gp = gp_aux[0]/(GC_ints[i+1] - GC_ints[i])
				total += gp
				gp_array.append(gp)
			ctg_probs = [gp/total for gp in gp_array]
			self.GC_probs[ctg] = ctg_probs

	#Following functions are required to compute the gene densities
	def set_gd(self, mapping_file):
		'''
		Computes the gene density for each contig
		'''
		def parse_mapping(mapping_file):
			'''
			Computing the gene coverage intervals for each contig

			The mapfile is to be provided in BLAST output fmt 6. It is a tab separated file.
			Each row of the file has the following information in tab separated format
			qseqid      query or gene sequence id (str)
			sseqid      subject or contig sequence id (str)
			pident      percentage of identical positions (float)
			length      alignment or overlap length (int)
			mismatch    number of mismatches (int)
			gapopen     number of gap openings (int)
			qstart      start of alignment in query (int)
			qend        end of alignment in query (int)
			sstart      start of alignment in subject (int)
			send        end of alignment in subject (int)
			evalue      expect value (float)
			bitscore    bit score (int)
			'''

			covg_int = {}
			with open(mapping_file, 'r') as map:
				line = next(map)
				while line:	
					tmp = line.split("\t")	
					ctg_id = tmp[1]
					sstart, send = tmp[8], tmp[9]
					if ctg_id not in covg_int:
						covg_int[ctg_id] = []
					if int(sstart) > int(send):
						covg_int[ctg_id].append((int(send), int(sstart)))
					else:
						covg_int[ctg_id].append((int(sstart), int(send)))
					line = next(map, None)
			return covg_int
		def get_union(intervals):
			'''
			Takes the gene covering intervals for a contig and finds their union
			The length of the union is used to compute gene coverage
			'''
			union = []
			for begin,end in sorted(intervals):
				if union and union[-1][1] >= begin-1:
					union[-1][1] = max(union[-1][1],end)
				else:
					union.append([begin,end])
			return union	
		self.gd = {}
		covg_int = parse_mapping(mapping_file)
		for ctg in self.assembly.get_ctgs():
			union = []
			if ctg in covg_int:
				union = get_union(covg_int[ctg])
			covered = 0
			for interval in union:
				covered += interval[1] - interval[0] + 1
			ctg_len = self.assembly.ctg_dict[ctg].get_len()
			self.gd[ctg] = covered/ctg_len
	
	def set_seeds(self, gd_thr, len_thr):
		'''
		Computes the seed value for each contig
		'''
		self.seeds = {}
		for ctg in self.assembly.get_ctgs():
			ctg_len = self.assembly.ctg_dict[ctg].get_len()	
			self.seeds[ctg] = 1 if (self.gd[ctg] >= gd_thr) and (ctg_len >= len_thr) else 0

	def get_GC_intervals(self):
		'''
		Returns list of GC interval endpoints
		'''
		return self.GC_intervals



#Initiating, computing and storing plasmid bin details.
class Plasmids:
	def __init__(self, source, pls_file, pident_thr = 0.95):
		'''
		Note: pident_thr (float): percent identity threshold for accepting contig-plasmid matches. 
		 				  		Contig belongs to plasmid only if plasmid covers >0.95 length of the contig. 
		plasmids (dict): Key: Plasmid ID, Value: List of pairs [ctg_id, mul]
		pls_lengths (dict): Key: Plasmid ID, Value: Length
		source (str): source of plasmid / plasmid bins, either ground truth or name of tool
		'''
		self.plasmids = {}
		self.pls_lengths = {}
		self.source = source
		self.parse_file(pls_file, pident_thr)

	def add_ctg_to_pls(self, pls, ctg, mul):
		'''
		Adds contig to list of contig representing plasmids. Each contig is represented as a pair [ctg_id, ctg_mul]. 
		'''
		if pls not in self.plasmids:
			self.plasmids[pls] = []
		self.plasmids[pls].append([ctg, mul])

	def set_pls_len(self, pls, pls_len):
		'''
		Sets length of plasmid
		'''
		if pls not in self.pls_lengths:
			self.pls_lengths[pls] = int(pls_len)

	def parse_file(self, pls_file, pident_thr):
		'''
		Reads the plasmid / plasmid bins file and stores the list of contigs with multiplicities
		If source is not PlasBin-Flow, default multiplicity is 1.		
		'''
		with open(pls_file, 'r') as pf:
			line = next(pf)
			while line:
				if line[0] != '#':			
					if self.source == 'ground_truth':
						mul = 1
						[pls, ctg, pident, pls_len, ctg_len] = line.split('\t')
						self.set_pls_len(pls, pls_len)
						if float(pident) >= pident_thr:
							self.add_ctg_to_pls(pls, ctg, mul)
					elif self.source == 'plasbin_flow':
						[bin_id, flow, GC_bin, ctgs] = filter(None,line.split('\t'))
						ctgs = ctgs.split(',')
						for pair in ctgs:
							ctg, mul = pair.split(':')[0], float(pair.split(':')[1])
							self.add_ctg_to_pls('pbf_'+bin_id, ctg, mul)
				line = next(pf, None)

	def get_pls_ids(self):
		'''
		Return list of plasmids
		'''
		return self.plasmids.keys()



#Initiating, computing and storing output specific to PlasBin-Flow.
class PBF_output():
	def __init__(self, source, pls_bins):
		'''
		Initiates attributes specific to PlasBin-Flow output
		flow (dict): Dictionary of flow values associated with each plasmid bin: Key: Bin_id, value: float
		GC_bin (dict): Dictionary of GC bin associated with each plasmid bin: Key: Bin_id, value: int
		GC_range (dict): Dictionary of GC range associated with each plasmid bin: Key: Bin_id, value: pair (list) of floats
		source (str): Method name	
		'''
		self.plasmids = Plasmids(source, pls_bins)
		self.flow = {}			#Dictionary of flow values associated with each plasmid bin: Key: Bin_id, value: float
		#self.GC_bin = {}		#Dictionary of GC bin associated with each plasmid bin: Key: Bin_id, value: int
		self.GC_range = {}		#Dictionary of GC range associated with each plasmid bin: Key: Bin_id, value: pair (list) of floats
		self.source = source	#Tool name: str
		self.parse_bins(pls_bins)

	def parse_bins(self, pls_bins):
		with open(pls_bins, 'r') as bins:
			line = next(bins)
			while line:
				if line[0] != '#':			
					[bin_id, flow, GC_bin, ctgs] = filter(None,line.split('\t'))
					self.flow['pbf_'+bin_id] = float(flow)
					#self.GC_bin['pbf_'+bin_id] = GC_bin
					self.GC_range['pbf_'+bin_id] = [float(GC_bin.split('-')[0]), float(GC_bin.split('-')[1])]
				line = next(bins, None)