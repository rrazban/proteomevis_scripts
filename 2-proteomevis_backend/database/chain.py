#!/usr/bin/python

help_msg = 'generate nodes data'

import sys, os
import numpy as np
import collections

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from properties import database
from write_sqlite3 import SQLite3
from output import print_next_step


def read_in_index():
	d = initialize_dict('dict')
	for organism in organism_list:
		pre_d = read_in('pdb', 'uniprot', filename='../0-identify_structure/0-identify_pdb/{0}/output.txt'.format(organism))
		pre_d = collections.OrderedDict(sorted(pre_d.items()))
		d[organism] = {i:pdb for i,pdb in enumerate(pre_d)}
	return d

def get_pdb_label(pdb):
	pdb1 = pdb
	pdb2 = pdb[:pdb.index('.')]
	pdb3 = pdb[pdb.index('.')+1:]
	return [pdb1, pdb2, pdb3]

def prepare_sql(d_org, d_index, d_val, protein_property_list, log_zero_list):
	line_list = []
	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p in range(total):
			line = [p, o]

			pdb = d_index[organism][p]
			pdb_label_list = get_pdb_label(pdb)				
			line.append(pdb)

			for pp, protein_property in enumerate(protein_property_list):
				try:
					val = float(d_val[organism][pp][pdb])
				except:
					line.append('')	#no value
					continue

				if protein_property=='dosage_tolerance':
					line.append(val)
				elif val!=0:
					line.append(np.log10(val))
				else:
					line.append(log_zero_list[pp])
			line_list.append(line)
	return line_list					


if __name__ == "__main__":
	help_message(help_msg, bool_org_dir = False)
	d_org = int2organism() 
	d_index = initialize_dict('dict')
	d_val = initialize_dict('list') 

	protein_property_list = ['length', 'abundance', 'evolutionary_rate', 'contact_density', 'PPI_degree', 'dosage_tolerance']
	log_zero_list = [-1, -1, -4, 1, -1]	#make into dict	#dont log dosage tolerance, already logged for yeast, ecoli is discrete
	for organism in organism_list:
		pre_d_i = read_in('pdb', 'uniprot', organism=organism)
		pre_d_i = collections.OrderedDict(sorted(pre_d_i.items()))
		d_index[organism] = {i:pdb for i,pdb in enumerate(pre_d_i)}

		d_ref = read_in('oln', 'pdb', organism=organism)
		for protein_property in protein_property_list:
			x_input = database(organism, protein_property)
			d = read_in(*x_input)

			d_subset = {pdb: d[oln] for oln,pdb in d_ref.iteritems() if oln in d}
			d_val[organism].append(d_subset)

	line_list = prepare_sql(d_org, d_index, d_val, protein_property_list, log_zero_list)

	columns = ['chain_id', 'species', 'pdb']
	columns.extend(protein_property_list)
	write_sqlite = SQLite3('proteomevis_chain', columns, line_list)
	write_sqlite.run()

	print_next_step('../')	
