#!/usr/bin/python

help_msg = 'generate protein function data'

import sys, os
import collections
import urllib, urllib2
import sqlite3

from chain import get_pdb_label

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from protein_property import database
from output import print_next_step


def get_info(organism):	#similar to 0-identify_structure/3-length_check
	url = 'http://www.uniprot.org/uniprot/'
	d_ref = read_in('uniprot', 'pdb', organism=organism)
	batch_size = 350		#491 is limit
	batch = len(d_ref)/batch_size
	
	uniprot_list = d_ref.keys()
	d=collections.defaultdict(list)
	for batch_i in range(batch+1):
		params = {'query':','.join(uniprot_list[(batch_i)*batch_size:(batch_i+1)*batch_size]), 'columns':'id,genes,comment(function),go(molecular function),comment(SUBCELLULAR LOCATION)','format':'tab'}	#doesnt work for single entry	#id querry doesnt work	#,comment(subcellular location),go'
		data = urllib.urlencode(params)
		request = urllib2.Request(url, data)
		response = urllib2.urlopen(request)
		label_list = next(response).split('\t')
		label_list = [label.rstrip() for label in label_list]
		for line in response:
			word_list = line.split('\t')
			word_list = [word.strip() for word in word_list]
			uniprot = word_list[label_list.index('Entry')]	#maybe iterate through a list
			genes = word_list[label_list.index('Gene names')]
			function = word_list[label_list.index('Function [CC]')]
			function2 = word_list[label_list.index('Gene ontology (molecular function)')]
			location = word_list[label_list.index('Subcellular location [CC]')]
			d[d_ref[uniprot]].extend([genes, function[len("FUNCTION: "):], function2, location[len('SUBCELLULAR LOCATION: '):]])
	return d

def writeout_sql(d_org, d_translate, d_index, d_val, args):
	table_name = 'proteomevis_inspect'
	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE IF EXISTS {0}'.format(table_name))	#have drop if present line
	c.execute('CREATE TABLE {0}(id,pdb,uniprot,genes,location,function1,function2,species)'.format(table_name)) 

	count=0
	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p in range(total):
			line_list = [count]
			count+=1	

			pdb = d_index[organism][p]
	#		pdb_label_list = get_pdb_label(pdb)				
	#		line_list.append(pdb_label_list[1])
			line_list.append(pdb)

			line_list.append(d_translate[organism][pdb])
			line_list.append(d_info[organism][pdb][0])
			line_list.append(d_info[organism][pdb][3])
			line_list.append(d_info[organism][pdb][1])
			line_list.append(d_info[organism][pdb][2])
			line_list.append(int(o))
			try:
				c.execute("INSERT INTO {0} VALUES {1}".format(table_name, tuple(line_list))) 
			except:
				if args.verbose:
					print 'function value cannot be read because has quotation marks'
					print line_list
					print 'function set to empty string'
				line_list[5]=''
				c.execute("INSERT INTO {0} VALUES {1}".format(table_name, tuple(line_list))) 
	conn.commit()
	conn.close()


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

	writeout_sql(d_org, d_translate, d_index, d_info, args)
	print_next_step('../')
