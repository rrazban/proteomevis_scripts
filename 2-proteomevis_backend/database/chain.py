#!/usr/bin/python

help_msg = 'generate nodes data'

import sys, os
import numpy as np
import collections
import sqlite3

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from protein_property import database
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

def writeout_sql(d_org, d_index, d_val, zero_list):
	table_name = 'proteomevis_chain'
	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
	c.execute('CREATE TABLE {0}(chain_id,species,pdb,length,abundance,evorate,conden,ppi,dostox)'.format(table_name))

	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p in range(total):
			line_list = [p, o]

			pdb = d_index[organism][p]
			pdb_label_list = get_pdb_label(pdb)				
			line_list.append(pdb)

			for pp in range(len(d_val[organism])):
				try:
					val = float(d_val[organism][pp][pdb])
				except:
					line_list.append('')	#no value
					continue
				if pp==len(d_val[organism])-1:
					line_list.append(val)	#dosage tolerance case
				elif val!=0:
					line_list.append(np.log10(val))
				else:
					line_list.append(zero_list[pp])
					
			c.execute("INSERT INTO {0} VALUES {1}".format(table_name, tuple(line_list))) #i am unable to pass NULL through python list

	for column in ['abundance', 'dostox', 'evorate']:
		c.execute("UPDATE {0} SET {1}=null where {1}=''".format(table_name, column))	

	conn.commit()
	conn.close()


if __name__ == "__main__":
	help_message(help_msg, bool_org_dir = False)	#add verbose option
	d_org = int2organism() 
	d_index = initialize_dict('dict')
	d_val = initialize_dict('list') 

	for organism in organism_list:
		pre_d_i = read_in('pdb', 'uniprot', organism=organism)
		pre_d_i = collections.OrderedDict(sorted(pre_d_i.items()))
		d_index[organism] = {i:pdb for i,pdb in enumerate(pre_d_i)}

		protein_property_list = ['length_pdb', 'abundance', 'evorate', 'contact_density', 'PPI', 'dosage_tolerance']	#add length/contact density
		zero_list = [-1, -1, -4, 1, -1]	#dont log dosage tolerance 	#already logged for yeast, ecoli is discrete
		d_ref = read_in('oln', 'pdb', organism=organism)
		for protein_property in protein_property_list:
			x_input = database(organism, protein_property)
			d = read_in(*x_input)

			d_subset = {pdb: d[oln] for oln,pdb in d_ref.iteritems() if oln in d}
			d_val[organism].append(d_subset)

	writeout_sql(d_org, d_index, d_val, zero_list)
	print_next_step('../')	
