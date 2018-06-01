#!/usr/bin/python

help_msg = "PDB complex file to PDB chains' files"

import os, sys, subprocess, glob
from Bio.PDB import PDBParser, PDBIO
from collections import defaultdict

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import organism
from output import print_next_step


DIR = '../../1-get_pdb/{0}'.format(organism) 

def untar():
	filename_list = glob.glob('{0}/*tar'.format(DIR))
	for filename in filename_list:
		item = '-pdb-bundle.tar'
		mapping_file = filename[:-len(item)] + '-chain-id-mapping.txt' 
		if not os.path.exists(mapping_file):	#do not  decompress if already done
		        subprocess.call(["tar", "-xvf", filename, "-C", DIR])	#alternatively, could just remove .tar file

def read_in_mapping(pdb):
	d = {}
	with open('{0}/{1}-chain-id-mapping.txt'.format(DIR, pdb), 'r') as rfile:
		next(rfile)
		for line in rfile:
			word_list = line.split()
			if len(word_list)==2:
				d[word_list[1]] = word_list[0]
	return d

class PDBChain():
	def __init__(self, pdbs):
		self.pdbs = pdbs
		self.d_input = defaultdict(list) 
		self.io = PDBIO()

	def get_pdb_chain(self):
		for pdb_chain in self.pdbs:
			if not os.path.exists("{0}.pdb".format(pdb_chain)):	#check if already present from previous run
				pdb = pdb_chain[:4]
				pdb_chain_id = pdb_chain[5:]
				self.d_input[pdb].append(pdb_chain_id)	

	def save_pdb_chain_file(self, translate):
		pre_chain_list = list(self.pdb_file.get_chains())	#store this length
		chain_list = [chain.get_id() for chain in pre_chain_list]
		want_chain_list = self.chain_list[:] 
		for chain in want_chain_list:
			if translate!=None:
				if translate[chain] in chain_list:
					want_i = chain_list.index(translate[chain])
					self.chain_list.remove(chain)
			else:
				want_i = chain_list.index(chain)
			self.io.set_structure(pre_chain_list[want_i])
			self.io.save(self.pdb_file.get_id() + "." + chain + ".pdb")

	def run(self):
		self.get_pdb_chain()
		for pdb, self.chain_list in self.d_input.iteritems():
			if os.path.exists("{0}/{1}.pdb".format(DIR, pdb)):
				self.pdb_file = PDBParser().get_structure(pdb, "{0}/{1}.pdb".format(DIR, pdb))	
				self.save_pdb_chain_file(None)
			else:
				pdb_bundle = glob.glob("{0}/{1}-pdb-bundle*pdb".format(DIR, pdb))
				for sub_file in pdb_bundle:
					translate_chain = read_in_mapping(pdb)
					self.pdb_file = PDBParser().get_structure(pdb, sub_file)
					self.save_pdb_chain_file(translate_chain)

if __name__ == "__main__":
	help_message(help_msg)
	untar()	

	d_input = read_in('pdb', 'uniprot', filename='pre_seq2struc')
	pdbchain = PDBChain(d_input.keys())
	pdbchain.run()

	print_next_step()
