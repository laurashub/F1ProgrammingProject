import requests, io, re, os, sys, argparse, random
import pandas as pd 
import matplotlib.pyplot as plt
#import numpy as numpy
from urllib.error import HTTPError
from Bio import Seq, Entrez, SeqIO

#import dna_features_viewer as dfv
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)

#TODO: Additional error checking, ie cant plot sequence before it's set

class DNABindingMap:
	"""This class finds the proteins that will bind to a dna sequence.

		Attributes: 
			sequence(str): the dna sequence used for binding
			binding_data(dict): the binding map of the dna """

	def __init__(self, sequence = None):
		"""
		Args:
			sequence(str): a dna sequence
		"""
		#TODO: handle multiple optional arguments
		self.sequence = None
		if sequence is not None:
			self.setSequence(sequence)
		self.binding_data = None

	def setSequence(self, sequence):
		""" This function generate the input DNA sequence that will be used for protein binding and store it in the sequence attribute.
			It can read the user input in three formats: DNA sequence string, fasta file, and genbank id

			Parameters:
				sequence(string): the DNA sequence / the path of a fasta file / the genbank id


			Example:
				setSequence("agtctgctgactgacgtttg")
				setSequence("../filename.fasta")
				setSequence("568815597") """
		seq = None
		if DNABindingMap.__checkDNA(sequence):
			seq = sequence
		elif re.search("fasta", sequence) != None:
			seq = DNABindingMap.__parseFasta(sequence)
		elif re.search("[^0-9]+", sequence) == None:
			seq = DNABindingMap.__getSeqFromGenbank(sequence)
		else:
			print("ERROR: Invalid sequence")
			raise DNABindingMapError
		if seq != None:
			seq = seq.lower()
		self.sequence = seq

	def getSequence(self):
		"""
			This function returns the dna sequence used for binding."""
		try:
			self.__checkSeq()
			return self.sequence
		except DNABindingMapError: 
			pass

	def __checkDNA(dna):

		check = re.search("[^atgcATGC]+", dna)
		if check == None:
			return True
		else:
			return False

	def __checkSeq(self):

		if self.sequence is None:
			print("ERROR: Current DNABindingMap has no sequence. Did you forget to call setSequence()?")
			raise DNABindingMapError

	def __checkBinding(self):
		if self.binding_data is None:
			print("ERROR: No binding data. Did you forget to call findBindingProteins()?")
			raise DNABindingMapError
		
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
			return str(seq_record.seq)
		except HTTPError:
			print("ERROR: The genbank id does not exist!")

	def findBindingProteins(self, dna_len_threshold = None):
		"""
			This function finds the proteins that will bind to the dna sequence.

			Parameters:
				dna_len_threshold(int): the shortest length of a matched dna sequence

			Returns:
				match(dict): 
					keys: the matched protein indexes in the input dataframe
					values: the starting and end indexes of the dna that the protein binds
			Example:
				match = findBindingProtein(dna_len_threshold = 5) """
		try:
			self.__checkSeq()

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
		except:
			pass

	def getBindingProteins(self, classification = None, subtype = None):
		"""This function get the specified group of proteins for dna binding.

		Parameters:
			classification(str): The classfication of proteins
			subtype(str): The sybtype of proteins

		Returns:
			A panda dataframe with the group of proteins

		Examples:
			getBindingProteins(classification = "Enzyme")
			getBindingProteins(classification = "Transcription factor", subtype = "Endonuclease")"""

		try:
			self.__checkBinding()
			DNABindingMap.__checkClassification(classification)
			DNABindingMap.__checkSubtype(subtype)

			binding_proteins = []
			for key, value in self.binding_data.items():
				name = df.at[key, '#PDB_ID']
				cur_classification = df.at[key, 'CLASSIFICATION']
				cur_subtype = df.at[key, 'SUBTYPE']
				binding_proteins.append([name, cur_classification, cur_subtype, value[0], value[1]])

			#filter based on parameters
			if classification is not None:
				binding_proteins = [x for x in binding_proteins if x[1] == classification]
			if subtype is not None:
				binding_proteins = [x for x in binding_proteins if x[2] == subtype]

			return pd.DataFrame(binding_proteins, columns=['#PDB_ID', 'CLASSIFICATION', 'SUBTYPE', 'START', 'STOP']) 
		
		except DNABindingMapError:
			pass

	def __checkClassification(classification):
		if classification is not None and classification not in df['CLASSIFICATION'].values:
			print("ERROR: Invalid classification")
			print("Allowed classifications: " + ", ".join(sorted(list(set(df['CLASSIFICATION'].values)))))
			raise DNABindingMapError

	def __checkSubtype(subtype):
		if subtype is not None and subtype not in df['SUBTYPE'].values:
			print("ERROR: Invalid subtype")
			print("Allowed subtypes: " + ", ".join(sorted(list(set(df['SUBTYPE'].values)))))
			raise DNABindingMapError

	def __checkPDB(self, pdb):
		pdbs = list(self.getBindingProteins()['#PDB_ID'].values)
		if pdb != None and pdb not in pdbs:
			print("ERROR: PDB \"" + str(pdb) + "\" not found")
			print("Current binding PDB IDs: " + ", ".join(pdbs))
			raise DNABindingMapError

	def __randomColor():
		r = lambda: random.randint(0,255)
		return '#{:02x}{:02x}{:02x}'.format(r(), r(), r())

	def showResults(self, classification = None, subtype = None, pdb = None, circular = False, show_sequence = False, filename = None):
		"""
		This function gives a cartoon plots showing where the proteins bind the dna.

		Parameters:
			classification(str): The classfication of proteins
			subtype(str): The sybtype of proteins
			pdb(str): PDB ID
			circular(bool): Display the results in circular plot
			show_sequence(bool): Show sequence on graph, discouraged for long sequences
			filename(str): filename to save graph to, otherwise will only show

		Examples:
			showResults(classification = "Enzyme")
			showResults(classification = "Transcription factor", subtype = "Endonuclease")"""

		global df

		try:
			self.__checkSeq()
			self.__checkBinding()

			DNABindingMap.__checkClassification(classification)
			DNABindingMap.__checkSubtype(subtype)
			self.__checkPDB(pdb)

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
				if pdb != None:
					if df.at[key, '#PDB_ID'] != pdb:
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

			if circular and show_sequence:
				print("ERROR: Unable to show sequence on circular record.")
				raise DNABindingMapError
				
			if circular:
				record = CircularGraphicRecord(sequence_length=len(self.sequence), sequence=self.sequence,  features=feats)
			else:
				record = GraphicRecord(sequence=self.sequence,  features=feats)

			fig, ax = plt.subplots()
			if show_sequence:
				record.plot_sequence(ax)
			ax1, _ = record.plot(ax=ax)

			title_constraints = []
			if classification is not None:
				title_constraints.append("Classification={0} ".format(classification))
			if subtype is not None:
				title_constraints.append("Subtype={0} ".format(subtype))
			if pdb is not None:
				title_constraints.append("PDB={0} ".format(pdb))
			ax.set_title(", ".join(title_constraints))

			fig.tight_layout()

			if filename is not None:
				filename = filename.strip('.png')
				plt.savefig(filename + '.png')
			else:
				plt.show()

		except DNABindingMapError:
			print("Unable to show results")

	def reset(self):
		self.sequence = None
		self.binding_data = None

	def __str__(self):
		temp_seq = ""
		temp_binding = []
		if self.sequence != None:
			temp_seq = self.sequence
		if self.binding_data != None:
			temp_binding = self.getBindingProteins()
		return("Seq: {0} \n Binding Data: \n  {1}".format(temp_seq, temp_binding))

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

class DNABindingMapError(Exception):
	pass

#### test cases ###
def test1():
	print("---Test 1---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	print(tester_map.getBindingProteins())
	tester_map.showResults()

def test2():
	print("---Test 2---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.showResults()

def test3():
	print("---Test 3---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	print(tester_map.getBindingProteins(classification = 'Enzyme'))
	tester_map.showResults(classification = 'Enzyme')

def test4():
	print("---Test 4---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.getBindingProteins(subtype = 'Zif268 Zinc Finger')
	tester_map.showResults(subtype = 'Zif268 Zinc Finger')

def test5():
	print("---Test 5---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(classification = 'fake classification')

def test6():
	print("---Test 6---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(pdb= '3BAM')

def test7():
	print("---Test 7---")
	tester_map = DNABindingMap()
	print(tester_map.getSequence())

def test8():
	print("---Test 8---")
	tester_map = DNABindingMap()
	tester_map.setSequence(sys.argv[1])
	tester_map.findBindingProteins()
	tester_map.showResults(pdb= 'xxxxxxx')

def test9():
	print("---Test 9---")
	tester_map = DNABindingMap()
	tester_map.setSequence("adskfhas;dfkhjalskdjhfalksdjhfalksdfjh")

if __name__ == "__main__":
	#input_dna = DNABindingMap._parge_args(sys.argv[1])

	test1()
	#test2()
	#test2()
	#test4()
	#test5()
	#test6()
	#test9()


