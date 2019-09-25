import requests, io, re, os, sys, argparse, random
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as numpy
from Bio import Seq
#import dna_features_viewer as dfv
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)

#TODO: Additional error checking, ie cant plot sequence before it's set

class DNABindingMap:

	def __init__(self, sequence = None, fasta = None, genbank_id = None):
		#TODO: handle multiple optional arguments
		self.sequence = ""
		self.binding_list = []

	def set_sequence(self, sequence):
		if DNABindingMap._checkDNA(sequence):
			self.sequence = sequence.lower()
			self.binding_list = {}
		else:
			print("ERROR: unrecognized characters in input sequence.")
			print(sequence)
			sys.exit(-1)

	def set_sequence_from_fasta(self, fasta_file):
		seq_from_fasta = DNABindingMap._parge_args(fasta_file)
		self.set_sequence(seq_from_fasta)

	def get_sequence(self):
		return self.sequence

	# Check if the input DNA only contains ATGC
	def _checkDNA(dna):
		check = re.search("[^atgcATGC]+", dna)
		print(check)
		if check == None:
			return True
		else:
			return False
		
	def _parge_args(file_name):
		fasta_file = open(os.path.abspath(file_name),'r')
		input_dna = ""
		for line in fasta_file.readlines():
			if line[0] != '>':
				input_dna += line.strip()
		return input_dna

	def findBindingProteins(self, dna_len_threshold = None):
		global df
		dna_seq = df["DNA_SEQS"].str.split(";").values
		match = {} 
		input_dna = self.sequence
		# Record the row number of the protein matched as keys
		#Record the starting and end positions of the matched dna as values
		for i in range(len(dna_seq)):
			for dna in dna_seq[i]:
				if ( dna_len_threshold == None or len(dna) >= dna_len_threshold) and dna in input_dna:
					start = input_dna.find(dna)
					end = start + len(dna) - 1
					match[i] = [start,end]
					break;
		self.binding_data = match
		print(match)
		return match


	def _random_color():
		r = lambda: random.randint(0,255)
		return '#{:02x}{:02x}{:02x}'.format(r(), r(), r())


	def showResults(self):
		global df
		# TODO: make PDB ids clickable, go to rcsb page?

		feats=[]

		#consolidate proteins that bind to same indeces
		bars = {}
		for key, value in self.binding_data.items():
			new_key = str(value[0]) + "," + str(value[1]) 
			if new_key in bars.keys():
				bars[new_key].append(df.at[key, '#PDB_ID'])
			else:
				bars[new_key] = [df.at[key, '#PDB_ID']]

		#generate bars
		for indices in bars.keys():
			start, end = indices.split(',')
			gf = GraphicFeature(start=int(start), end=int(end), strand=+1, color=DNABindingMap._random_color(),
                   label=", ".join(bars[indices]))
			feats.append(gf)

		record = GraphicRecord(sequence=self.sequence,  features=feats)

		fig, ax = plt.subplots()
		#record.plot_sequence(ax)
		ax1, _ = record.plot(ax=ax)
		#record.plot_sequence(ax1, {'size': 11})

		fig.tight_layout()
		plt.show()

def dbRead():
	url = ("http://melolab.org/pdidb/web/download/pdidb.txt")
	test = requests.get(url).content
	df = pd.read_csv(io.StringIO(test.decode('utf-8')), sep='\t', skiprows=5)
	#df1 = df[['#PDB_ID', 'PROT_SEQS', 'DNA_SEQS']].copy()
	lower = lambda x: x.lower()
	df['DNA_SEQS']=df['DNA_SEQS'].apply(lower)
	return df

	# Index(['#PDB_ID', 'ENTRY_ID', 'PUBMED_ID', 'RESOLUTION', 'SPECIES',
      #  'MULTICOMPLEXITY', 'NO_INTERFACES', 'ASYMMETRIC_UNIT',
      #  'BIOLOGICAL_UNITS', 'IS_CRUCIFORM', 'HAS_ZDNA', 'HAS_WATER',
      #  'CLASSIFICATION', 'TYPE', 'SUBTYPE', 'COMPLEX_ID', 'INTERFACE_ID',
      #  'DNA_DOUBLE_STRAND', 'DNA_SINGLE_STRAND', 'STICKY_ENDS', 'FLIPPED_BASE',
      #  'NICKED_DNA', 'GAPPED_DNA', 'OPEN_DNA', 'MODIFIED_DNA',
      #  'PROTEIN_MONOMERS', 'MULTIMERIZATION', 'PROT_PROT_INTERACTION',
      #  'DNA_PROT_INTERACTION', 'PROT_CHAIN_IDS', 'PROT_SEQS', 'DNA_CHAIN_IDS',
      #  'DNA_SEQS', 'SEQ_GROUPS', 'INTERFACE_GROUP', 'NUMBER_EFF_CONTACTS',
      #  'GROOVE_CONTACTS(MA;MI;BB;NA)',
      #  'TYPE_CONTACTS(1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20)'],
      # dtype='object')

df = dbRead()

if __name__ == "__main__":
	#input_dna = DNABindingMap._parge_args(sys.argv[1])
	
	tester_map = DNABindingMap()
	tester_map.set_sequence_from_fasta(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults()


