import requests, io
import pandas as pd 
import matplotlib.pyplot as mp

def db_read():
	url = ("http://melolab.org/pdidb/web/download/pdidb.txt")
	test = requests.get(url).content
	df = pd.read_csv(io.StringIO(test.decode('utf-8')), sep='\t', skiprows=5)
	df1 = df[['#PDB_ID', 'PROT_SEQS', 'DNA_SEQS']].copy()
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

def find_binding_proteins(seq, database):
	return None

def show_results(seq, binding_data):
	return None

db = db_read()
print(db.head())

