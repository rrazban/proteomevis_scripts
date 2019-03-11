#!/usr/bin/python

help_msg = 'remove pre_seq2struc proteins whose PDB length is 0.8 times less than reported UniProt chain length' #(selection criterion 2) 

import sys, os
from Bio.PDB import PDBParser
from Bio import PDB                                                       
import urllib, urllib2

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import organism
from properties import database
from output import writeout, database_update_needed


def parse_pdb_length(name):
	pdb = PDBParser().get_structure(name, "../../../0-identify_structure/2-get_pdb_chain/{0}/{1}.pdb".format(organism, name))
	chain = list(pdb.get_chains())[0]	#only 1 chain present	
	return len([_ for _ in chain.get_residues() if PDB.is_aa(_)])	#omits missing residues


class UniProt2PDB():
	def __init__(self, verbose):
		self.verbose = verbose
		self.d_ref = read_in('oln', 'pdb', 'pre_seq2struc')
		self.d_ref2 = read_in('oln', 'uniprot', 'pre_seq2struc')

		self.d_input = read_in(*database(organism, 'length'))
		self.d_output = {}

	def get_pdb_length(self):
		for oln, pdb in self.d_ref.iteritems():
			if oln in self.d_input:	#some have no reported lengths in UniProt
				pdb = self.d_ref[oln]
				pdb_length = parse_pdb_length(pdb)
				frac = pdb_length/float(self.d_input[oln])
				if frac > 0.8 and frac < 1.2:	#need upperbound as well
					self.d_output[self.d_ref2[oln]] = [pdb, oln, pdb_length]

	def run(self):
		self.get_pdb_length()
		return self.d_output


if __name__ == "__main__":
	args = help_message(help_msg, bool_add_verbose=True)
	uniprot2pdb = UniProt2PDB(args.verbose)
	d_output = uniprot2pdb.run()

	filename = 'seq2struc'
	writeout(['uniprot', 'pdb', 'oln', 'pdb_length'], d_output, filename = 'new_{0}'.format(filename))
	update_bool = database_update_needed(filename)
	if update_bool: print 'keep old_seq2struc.txt for efficient running of extra tm.sid jobs'
