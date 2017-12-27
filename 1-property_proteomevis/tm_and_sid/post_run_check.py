#!/usr/bin/python
import sys, os
from collections import defaultdict

sys.path.append('../../../utlts/')
from read_in_file import read_in
from parse_data import organism, organism_list
from protein_property import database
from output import today, database_update_needed


def rewrite(d_ref, organism, one_organism_bool):
	if one_organism_bool:	
		d_ref1 = d_ref2 = d_ref[organism]
	else:
		d_ref1 = d_ref['yeast']
		d_ref2 = d_ref['ecoli']

	with open('PDB_{0}.txt'.format(today), 'w') as wfile:
		with open('PDB.txt', 'r') as rfile:
			labels = next(rfile)
			label_list = labels.split()
			wfile.write(labels)
			for line in rfile:
				word_list = line.split()
				pdb1 = word_list[label_list.index('pdb1')]
				pdb2 = word_list[label_list.index('pdb2')]
				if pdb1 in d_ref1 and pdb2 in d_ref2:
					wfile.write(line)
def parse(d):
	new_d = defaultdict(list)
	for pdb1pdb2, tm in d.iteritems():
		pdb_list = pdb1pdb2.split(',')
		new_d[pdb_list[0]].append(pdb_list[1])
		new_d[pdb_list[1]].append(pdb_list[0])
	return new_d

def find_duplicates(d):
	new_d = {}
	for pdb, pdb_list in d.iteritems():
		if len(pdb_list)!=len(set(pdb_list)):
			print '{0} has duplicate pairs'.format(pdb)

def reference():
	d_ref={}
	for organism in organism_list:
		d_ref[organism] = read_in('pdb', 'uniprot', organism=organism)
	return d_ref

def check(d_ref, d):
	need_extra_pairs_bool = False
	if organism in organism_list:
		num = len(d_ref[organism])-1
		for pdb,pdb_list in d[organism].iteritems():
			if len(pdb_list)!=num:
				print "{0} has {1} pairs. Should be {2}".format(pdb, len(pdb_list), num)
				need_extra_pairs_bool = True 
	else:
		for pdb,pdb_list in d[organism].iteritems():
			if pdb in d_ref['yeast']:
				num = len(d_ref['ecoli'])	
			else:
				num = len(d_ref['yeast'])	
				
			if len(pdb_list)!=num:
				print "{0} has {1} pairs. Should be {2}".format(pdb, len(pdb_list), num)
				need_extra_pairs_bool = True 
	return need_extra_pairs_bool


if __name__ == "__main__":
	d_ref = reference() 
	rewrite(d_ref, organism, organism in organism_list)
	database_update_needed("PDB")	#not clear, cuz of message to move on

	raw_d = read_in(*database(organism, 'tm'))

	d = {}
	d[organism] = parse(raw_d)
	find_duplicates(d[organism])
	if check(d_ref, d):
		print "\nTo obtain missing pairs:"
		print '1) ../pre_run_extra.py to generate extra.txt' 
		print '2) ../run --extra'
		print '2a) if in yeast_ecoli directory, also ../run --EXTRA'
		print '3) manually copy output.txt contents and paste at the end of PDB.txt'
	else:
		print 'All pairs are present!'
