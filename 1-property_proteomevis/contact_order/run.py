#!/usr/bin/python

help_msg = 'calculate contact density from PDB structure file'

import os, sys, glob
import imp
from Bio.PDB import NeighborSearch, PDBParser, Atom, Residue, Polypeptide
from Bio import PDB
import numpy as np


CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message, false_or_true
from read_in_file import read_in
from parse_data import organism
from output import writeout, print_next_step


def contact_order(contact_matrix):
	CO = 0
	for res1, contact_res in enumerate(contact_matrix):
		for res2, contact in enumerate(contact_res):
			if contact:
				CO+= abs(res1-res2)
	return CO / (float(len(contact_matrix)*contact_matrix.sum()))		


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

	module = imp.load_source("run", "../../contact_density/run.py")	#normal import doesnt work
	for pdb, oln in d_input.iteritems():
		protein_contact = module.ProteinContact(pdb, contact_defn)
		CO = contact_order(protein_contact.contact_matrix())
		if organism=='protherm':
			d_output[d_input1[pdb]] = CO 
			x_name = 'uniprot'
		else:
			d_output[oln] = CO
			x_name = 'oln'

	filename = 'PDB'
	if method:
		filename+='_shakh'
	writeout([x_name, 'contact_order'], d_output, filename = '{0}{1}'.format(filename, extra))
	print_next_step()



