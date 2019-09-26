import requests, io, re, os, sys, argparse, random
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as numpy
from urllib.error import HTTPError
from Bio import Seq, Entrez, SeqIO

#import dna_features_viewer as dfv
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)

#TODO: Additional error checking, ie cant plot sequence before it's set

class DNABindingMap:

	def __init__(self, sequence = None):
		#TODO: handle multiple optional arguments
		self.sequence = ""
		self.binding_data = None

	def setSequence(self, sequence):
		seq = None
		if DNABindingMap.__checkDNA(sequence):
			seq = sequence
		elif re.search("fasta", sequence) != None:
			seq = DNABindingMap.__parseFasta(sequence)
		elif re.search("[^0-9]+", sequence) == None:
			seq = DNABindingMap.__getSeqFromGenbank(sequence)
		else:
			print("The input does not match any format.")
			sys.exit(-1)
		self.sequence = seq.lower()

	def getSequence(self):
		return self.sequence

	# Check if the input DNA only contains ATGC
	def __checkDNA(dna):
		check = re.search("[^atgcATGC]+", dna)
		if check == None:
			return True
		else:
			return False
		
	def __parseFasta(file_name):
		fasta_file = open(os.path.abspath(file_name),'r')
		input_dna = ""
		for line in fasta_file.readlines():
			if line[0] != '>':
				input_dna += line.strip()
		return input_dna

	# Get the DNA sequence from NCBI using genbankid
	def __getSeqFromGenbank(genbankid):
		Entrez.email = 'A.N.Other@example.com'
		try:
			with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=genbankid) as handle:
				seq_record = SeqIO.read(handle, "fasta")
			return seq_record.seq
		except HTTPError:
			print("The genbank id does not exist!")

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
		return match

	def getBindingProteins(self, type = None):
		binding_proteins = []
		for key, value in self.binding_data.items():
			name = df.at[key, '#PDB_ID']
			classification = df.at[key, 'CLASSIFICATION']
			subtype = df.at[key, 'SUBTYPE']
			binding_proteins.append((name, classification, subtype, value))
		return binding_proteins


	def __randomColor():
		r = lambda: random.randint(0,255)
		return '#{:02x}{:02x}{:02x}'.format(r(), r(), r())

	def showResults(self, classification = None, subtype = None):
		global df

		if self.binding_data == None:
			print("ERROR: No binding data. Did you forget to call findBindingProteins()?")
			return None

		if classification is not None and classification not in df['CLASSIFICATION'].values:
			print("ERROR: Invalid classification")
			print("Allowed classifications: " + ", ".join(set(df['CLASSIFICATION'].values)))
			return

		if subtype is not None and subtype not in df['SUBTYPE'].values:
			print("ERROR: Invalid subtype")
			print("Allowed subtypes: " + ", ".join(set(df['SUBTYPE'].values)))
			return
		# TODO: make PDB ids clickable, go to rcsb page?

		feats=[]

		#consolidate proteins that bind to same indeces
		bars = {}
		for key, value in self.binding_data.items():
			new_key = str(value[0]) + "," + str(value[1]) 
			if classification != None:
				if df.at[key, 'CLASSIFICATION'] != classification:
					continue
			if subtype != None:
				if df.at[key, 'SUBTYPE'] != subtype:
					continue
			if new_key in bars.keys():
				bars[new_key].append(df.at[key, '#PDB_ID'])
			else:
				bars[new_key] = [df.at[key, '#PDB_ID']]

		if len(bars.keys()) == 0:
			print("Warning: no binding proteins to show")

		#generate bars
		for indices in bars.keys():
			start, end = indices.split(',')
			gf = GraphicFeature(start=int(start), end=int(end), strand=+1, color=DNABindingMap.__randomColor(),
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
	print(df['SUBTYPE'].head())
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

def test1():
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults()

def test2():
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.showResults()

def test3():
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(classification = 'Enzyme')

def test4():
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(subtype = 'Zinc finger')

def test5():
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(classification = 'fake classification')

if __name__ == "__main__":
	#input_dna = DNABindingMap._parge_args(sys.argv[1])
	test1()
	test2()
	test3()
	test4()
	test5()


