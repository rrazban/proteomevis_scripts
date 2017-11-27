#!/usr/bin/python

#selection criterion 2: remove all those identified in pre_seq2struc whose length is less than
#			80% of UniProt chain length 

import sys, os, collections
from Bio.PDB import PDBParser
from Bio import PDB                                                       
import urllib, urllib2

sys.path.append('../../../utlts/')
from read_in_file import read_in, set_up_read_in
from parse_data import organism
from output import writeout, database_update_needed


def get_length(d, d_len):
	d_output = {}	
	for name, uniprot in d.iteritems():
		pdb = PDBParser().get_structure(name, "../../../0-identify_structure/2-get_pdb_chain/{0}/{1}.pdb".format(organism, name))
		for c, chain in enumerate(pdb.get_chains()):	
			if c>0:	#should just be one chain per file
				print "WARNING multiple chains found: {0}".format(name)
			length = len([_ for _ in chain.get_residues() if PDB.is_aa(_)])	#omits missing residues
			if 0.8<length/float(d_len[uniprot]):#<1.2:	#allow for overlength in case pdb structure file 
									#includes pre-post-translationally modified residues 
									#-rare (0 cases last time i checked)
				d_output[name] = length
	return d_output

def get_info(d_ref):
	url = 'http://www.uniprot.org/uniprot/'
	batch_size = 350		#491 is limit
	batch = len(d_ref)/batch_size
	
	uniprot_list = d_ref.values()
	d={}
	for batch_i in range(batch+1):
		params = {'query':','.join(uniprot_list[(batch_i)*batch_size:(batch_i+1)*batch_size]), 'columns':'id,feature(CHAIN)','format':'tab'}	
		data = urllib.urlencode(params)
		request = urllib2.Request(url, data)
		response = urllib2.urlopen(request)
		label_list = next(response).split('\t')
		label_list = [label.rstrip() for label in label_list]
		for line in response:
			word_list = line.split()
			word_list = [word.strip() for word in word_list]
			uniprot = word_list[label_list.index('Entry')]
			if len(word_list)==1:	#does not capture UniProt peptide case
				print 'No chain found: {0}'.format(word_list)
				length = 1000000000
			elif word_list[label_list.index('Chain')+1]=='?':
				print 'No starting residue for chain: {0}'.format(word_list)
				length = int(word_list[label_list.index('Chain')+2])
			else:	
				length = int(word_list[label_list.index('Chain')+2])-int(word_list[label_list.index('Chain')+1])+1
			d[uniprot] = length 
	return d

def prepare_writeout(d_ref, d_output):
	d = {} 
	d_ref2 = read_in('pdb', 'oln', 'pre_seq2struc')
	for pdb, length in d_output.iteritems():
		d[pdb] = [d_ref[pdb], d_ref2[pdb], length]
	return d


if __name__ == "__main__":
	d_ref = read_in('pdb', 'uniprot', 'pre_seq2struc')
	d_len = get_info(d_ref) 
	d_output = get_length(d_ref, d_len)
	d = prepare_writeout(d_ref, d_output)
	filename = 'seq2struc'
	writeout(['pdb', 'uniprot', 'oln', 'length'], collections.OrderedDict(sorted(d.items())), filename=filename, date_bool=True)
	update_bool = database_update_needed(filename)
	if update_bool: print 'keep current seq2struc.txt as old_seq2struc.txt for efficient running of extra tm.sid jobs'
