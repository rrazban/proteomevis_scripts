#!/usr/bin/python

help_msg = 'generate protein function data'

import sys, os
import collections

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from properties import database
from uniprot_api import UniProtAPI
from write_sqlite3 import SQLite3
from output import print_next_step


def get_info(organism):
	d_ref = read_in('uniprot', 'pdb', organism = organism)
	
	d = {}
	columns = ['id','genes','comment(function)','go(molecular function)','comment(SUBCELLULAR LOCATION)']
	uniprot_api = UniProtAPI(columns)
	label_list, response = uniprot_api.organism_info(organism = organism)

	for line in response:
		word_list = line.split('\t')
		word_list = [word.strip() for word in word_list]
		uniprot = word_list[label_list.index('Entry')]	#maybe iterate through a list
		if uniprot in d_ref:
			genes = word_list[label_list.index('Gene names')]
			function = word_list[label_list.index('Function [CC]')]
			function2 = word_list[label_list.index('Gene ontology (molecular function)')]
			location = word_list[label_list.index('Subcellular location [CC]')]
			if '"' in function:	#sqlite3 cant handle "" marks
				function = function.replace('"', '')
	
			d[d_ref[uniprot]] = [genes, location[len('SUBCELLULAR LOCATION: '):], function[len("FUNCTION: "):], function2]
	return d

def prepare_sql(d_org, d_translate, d_index, d_val, args):	#make this general for all to import
	line_list = []
	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p in range(total):
			pdb = d_index[organism][p]
			line = [p, pdb, d_translate[organism][pdb]]
			line.extend(d_info[organism][pdb])
			line.append(int(o))
			line_list.append(line)
	return line_list

if __name__ == "__main__":
	args = help_message(help_msg, bool_add_verbose=True, bool_org_dir = False)	#add verbose option
	d_org = int2organism() 
	d_translate = initialize_dict('dict')
	d_index = initialize_dict('dict')
	d_info = initialize_dict('dict')

	for organism in organism_list:
		pre_d_i = read_in('pdb', 'uniprot', organism=organism)
		pre_d_i = collections.OrderedDict(sorted(pre_d_i.items()))
		d_translate[organism] = pre_d_i
		d_index[organism] = {i:pdb for i,pdb in enumerate(pre_d_i)}
		d_info[organism] = get_info(organism)

	line_list = prepare_sql(d_org, d_translate, d_index, d_info, args)

	columns = ['chain_id','pdb','uniprot','genes','location','function1','function2','species']
	write_sqlite = SQLite3('proteomevis_inspect', columns, line_list)
	write_sqlite.run()

	print_next_step('../')
