#!/usr/bin/python
import sys, os
import numpy as np
import collections
import sqlite3

sys.path.append('../utlts/')
from read_in_file import read_in
from parse_data import initialize_dict, organism_list, int2organism
from protein_property import database


def read_in_ppi_partners():
	d = initialize_dict("dict")
	for organism in organism_list:
		with open("../1-property_proteomevis/ppi_partner/{0}/ppi_partners.txt".format(organism), "r") as rfile:
			label_list = next(rfile).split('\t')
			label_list = [x.rstrip() for x in label_list]
			for line in rfile:
				word_list = line.split()
				protein = word_list[label_list.index('protein')]	
				ppi = word_list[label_list.index('protein partners'):]
				d[organism][protein] = ppi
	return d	

def writeout_sql(d_org, d_index, d_ppi, d_val):
	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE proteomevis_edge')
	c.execute('CREATE TABLE proteomevis_edge(id,sourceID,species,targetID,tm,sid,ppi,mutant)')

	count=0
	for o in range(len(d_org)):
		organism = d_org[o]
		total = len(d_index[organism])
		for p1 in range(total):
			pdb1 = d_index[organism][p1]
			for p2 in range(p1+1, total):
				pdb2 = d_index[organism][p2]
				
				ppi_bool = 0
				if pdb1 in d_ppi[organism]:
					if pdb2 in d_ppi[organism][pdb1]:	#d_ppi might not have pdb1	
						ppi_bool = 1

				pdb_pair = pdb1 +','+ pdb2
				if pdb_pair not in d_val[organism][0]:
					pdb_pair = pdb2 +','+ pdb1
				line_list = [count, p1, o, p2, float(d_val[organism][0][pdb_pair]), float(d_val[organism][1][pdb_pair]), ppi_bool, 0]
				c.execute("INSERT INTO proteomevis_edge VALUES {0}".format(tuple(line_list))) 
				count+=1
	conn.commit()
	conn.close()


if __name__ == "__main__":
	d_org = int2organism() 
	d_index = initialize_dict('dict')
	d_val = initialize_dict('list')

	for organism in organism_list:
		pre_d_i = read_in('pdb', 'uniprot', organism=organism)
		pre_d_i = collections.OrderedDict(sorted(pre_d_i.items()))
		d_index[organism] = {i:pdb for i,pdb in enumerate(pre_d_i)}

		for x in ['tm', 'sid']:
			d_val[organism].append(read_in(*database(organism, x)))
	d_ppi = read_in_ppi_partners()
	writeout_sql(d_org, d_index, d_ppi, d_val)
