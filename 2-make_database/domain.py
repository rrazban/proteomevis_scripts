#!/usr/bin/python
import sys, os
import collections
import urllib, urllib2
import sqlite3

from chain import get_pdb_label

sys.path.append('../utlts/')
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from protein_property import database


def get_info(organism):	#similar to 0-identify_structure/3-length_check
	url = 'http://www.uniprot.org/uniprot/'
	d_ref = read_in('uniprot', 'pdb', organism=organism)
	batch_size = 350		#491 is limit
	batch = len(d_ref)/batch_size
	
	uniprot_list = d_ref.keys()
	d=collections.defaultdict(list)
	for batch_i in range(batch+1):
		params = {'query':','.join(uniprot_list[(batch_i)*batch_size:(batch_i+1)*batch_size]), 'columns':'id,genes,comment(function)','format':'tab'}	#doesnt work for single entry	#id querry doesnt work	#,comment(subcellular location),go'
		data = urllib.urlencode(params)
		request = urllib2.Request(url, data)
		response = urllib2.urlopen(request)
		label_list = next(response).split('\t')
		label_list = [label.rstrip() for label in label_list]
		for line in response:
			word_list = line.split('\t')
			word_list = [word.strip() for word in word_list]
			uniprot = word_list[label_list.index('Entry')]
			genes = word_list[label_list.index('Gene names')]
			function = word_list[label_list.index('Function [CC]')]
			d[d_ref[uniprot]].extend([genes, function[len("FUNCTION: "):]])
	return d

def writeout_sql(d_org, d_translate, d_index, d_val):
	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE proteomevis_domain')
	c.execute('CREATE TABLE proteomevis_domain(id,domain,uniprot,genes,details,function1,function2,invis,obsolete,species)') 

	count=0
	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p in range(total):
			line_list = [count]
			count+=1	

			pdb = d_index[organism][p]
			pdb_label_list = get_pdb_label(pdb)				
			line_list.append(pdb_label_list[1])

			line_list.append(d_translate[organism][pdb])
			line_list.append(d_info[organism][pdb][0])
			line_list.append("")
			line_list.append(d_info[organism][pdb][1])
			line_list.extend(["",1,""])
			line_list.append(int(o))
			try:
				c.execute("INSERT INTO proteomevis_domain VALUES {0}".format(tuple(line_list))) 
			except:
				print line_list
				line_list[5]=''
				print line_list
				c.execute("INSERT INTO proteomevis_domain VALUES {0}".format(tuple(line_list))) 
	conn.commit()
	conn.close()


if __name__ == "__main__":
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

	writeout_sql(d_org, d_translate, d_index, d_info)
