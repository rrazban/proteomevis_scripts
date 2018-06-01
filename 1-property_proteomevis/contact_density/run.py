#!/usr/bin/python

help_msg = 'calculate contact density from PDB structure file'

import os, sys, glob
from Bio.PDB import NeighborSearch, PDBParser, Atom, Residue, Polypeptide
from Bio import PDB
import numpy as np


CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message, false_or_true
from read_in_file import read_in
from parse_data import organism
from protein_property import database
from output import writeout, print_next_step


class Protein(object):
	def __init__(self, pdb, contact_defn):
		self.pdb = pdb 
		self.contact_defn = contact_defn

		pdbs_dir = '../../../0-identify_structure/2-get_pdb_chain'
		self.structure = PDBParser().get_structure('X', "{0}/{1}/{2}.pdb".format(pdbs_dir, organism, pdb))
		self.residues = []	#only consider actual residues
		self.atoms = []

		self.parse_structure()

	def parse_structure(self):
		for residue in self.structure.get_residues():
		#	if PDB.is_aa(residue, standard=True):	 # only the standard 20
			if PDB.is_aa(residue):
				res = residue.id[1]
				if res not in self.residues:	#dont doublecount mutated residues
					self.residues.append(res)
					self.atoms.extend(atoms_method(self.contact_defn, residue))

	def get_residues(self):
		return self.residues
	
	def get_max_residue(self):
		return max(self.residues)+1

def atoms_method(contact_defn, residue):
	if contact_defn=='Bloom':
		atoms = residue.get_atoms()
	elif contact_defn=='Shakh':
		if 'CB' in residue:
			atoms = [residue['CB']]
		elif Residue.Residue.get_resname(residue)=='GLY':	#glycine case
			if 'CA' in residue:
				atoms = [residue['CA']]
			else:
				atoms = []
		else:
			atoms = []
	else:
		print 'method index not supported'
		sys.exit()
	return atoms

class ProteinContact(Protein):
	def __init__(self, pdb, contact_defn):
		super(ProteinContact, self).__init__(pdb, contact_defn)

		self.residues_max = self.get_max_residue() 

	def contact_pairs(self, radius = ''):	#define radius here in case want to change but keep
		if not radius:
			d_radii = {"Bloom": 4.5, "Shakh": 7.5}	#in Angstroms
			radius = d_radii[self.contact_defn]
		ns = NeighborSearch(self.atoms)
		self.contact_pairs_list = ns.search_all(radius, level = 'R')

	def contact_matrix(self, radius = ''):
		self.contact_pairs(radius)
		M = np.zeros((self.residues_max, self.residues_max), dtype=int)	#in case radius changed

		for contact in self.contact_pairs_list:
			res1 = contact[0].id[1]
			res2 = contact[1].id[1]
			if not abs(res1 - res2) in [1, 0]:	#no nearest neighbors
				M[res1][res2] = 1
		return M + M.T
	

if __name__ == "__main__":
	help_message(help_msg)
	extra = ''
	method = false_or_true("Calculate contact density like Shakh2006 [default Zhou2008]?")
	if false_or_true("Relax selection criterion 2"):
		extra += 'pre_output'

	contact_defn = ['Bloom', 'Shakh'][method]
	d_input = read_in('pdb', 'oln', filename = extra)
	d_input1 = read_in('pdb', 'uniprot', filename = extra)
	d_output = {}
	for pdb, oln in d_input.iteritems():
		protein_contact = ProteinContact(pdb, contact_defn)
		residues = protein_contact.get_residues()
		contact_density = protein_contact.contact_matrix().sum() / float(len(residues))
		if organism=='protherm':
			d_output[d_input1[pdb]] = contact_density
			x_name = 'uniprot'
		else:
			d_output[oln] = contact_density 
			x_name = 'oln'

	filename = 'PDB'
	if method:
		filename+='_shakh'
	writeout([x_name, 'contact_density'], d_output, filename = '{0}{1}'.format(filename, extra))
	print_next_step()



