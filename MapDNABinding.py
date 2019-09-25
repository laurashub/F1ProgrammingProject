import requests, io, re, os, sys, argparse, random
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as numpy
from urllib.error import HTTPError
from Bio import Seq, Entrez, SeqIO

#import dna_features_viewer as dfv
from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)

# Generate the DNA sequence(s) for matching
# Take dna sequence(string), fasta file, or genbackid as input
def get_seq(input):
	seq = None
	if _checkDNA(input):
		seq = input
	elif re.search("fasta", input) != None:
		seq = parse_fasta(input)
	elif re.search("[^0-9]+", input) == None:
		seq = get_seq_from_genbank(input)
	else:
		print("The input does not match any format.")
	return seq

# Get the DNA sequence from NCBI using genbankid
def get_seq_from_genbank(genbankid):
	Entrez.email = 'A.N.Other@example.com'
	try:
	    with Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=genbankid) as handle:
	        seq_record = SeqIO.read(handle, "fasta")
	    return seq_record.seq
	except HTTPError:
	    print("The genbank id does not exist!")

# Parse the input fasta file
def parse_fasta(file_name):
	fasta_file = open(os.path.abspath(file_name),'r')
	input_dna = ""
	for line in fasta_file.readlines():
		if line[0] != '>':
			input_dna += line.strip()
	return input_dna


def dbRead():
	url = ("http://melolab.org/pdidb/web/download/pdidb.txt")
	test = requests.get(url).content
	df = pd.read_csv(io.StringIO(test.decode('utf-8')), sep='\t', skiprows=5)
	df1 = df[['#PDB_ID', 'PROT_SEQS', 'DNA_SEQS']].copy()
	lower = lambda x: x.lower()
	df1['DNA_SEQS']=df1['DNA_SEQS'].apply(lower)
	return df1

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

def findBindingProteins(input_dna, df, dna_len_threshold = None):
    dna_seq = df["DNA_SEQS"].str.split(";").values
    match = {} 
    # Record the row number of the protein matched as keys
    # Record the starting and end positions of the matched dna as values
    for i in range(len(dna_seq)):
    	for dna in dna_seq[i]:
    		if (len(dna) >= dna_len_threshold  or dna_len_threshold == None) and dna in input_dna:
    			start = input_dna.find(dna)
    			end = start + len(dna) - 1
    			match[i] = [start,end]
    			break;
    return match

# Check if the input DNA only contains ATGC
def _checkDNA(dna):
      check = re.search("[^atgcATGC]+", dna)
      if check == None:
            return True
      else:
            return False

def random_color():
	r = lambda: random.randint(0,255)
	return '#{:02x}{:02x}{:02x}'.format(r(), r(), r())


def showResults(seq, df, binding_data):
	# TODO: make PDB ids clickable, go to rcsb page?

	feats=[]

	#consolidate proteins that bind to same indeces
	bars = {}
	for key, value in binding_data.items():
		new_key = str(value[0]) + "," + str(value[1]) 
		if new_key in bars.keys():
			bars[new_key].append(df.at[key, '#PDB_ID'])
		else:
			bars[new_key] = [df.at[key, '#PDB_ID']]

	#generate bars
	for indices in bars.keys():
		start, end = indices.split(',')
		gf = GraphicFeature(start=int(start), end=int(end), strand=+1, color=random_color(),
                   label=", ".join(bars[indices]))
		feats.append(gf)

	record = GraphicRecord(sequence=seq, features=feats)

	fig, ax = plt.subplots()
	ax1, _ = record.plot(ax=ax)
	record.plot_sequence(ax1)

	fig.tight_layout()
	plt.savefig("result.png")


df = dbRead()
input_dna = get_seq('6873')
print(input_dna)
#input_dna = "aaaaaaaaaagtcgcagcgtgggaccgtagctgaGTaattaCGgcagcgcac"
# if _checkDNA(input_dna):
# 	dna_lower = input_dna.lower()
# 	matched_proteins = findBindingProteins(dna_lower, df, dna_len_threshold = 5 )
# else:
# 	print("ERROR: unrecognized characters in input sequence.")
# 	sys.exit(-1)

# showResults(dna_lower, df, matched_proteins)

