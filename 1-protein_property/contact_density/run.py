#!/usr/bin/python
import os, sys
from Bio.PDB import NeighborSearch, PDBParser, Atom, Residue
from Bio import PDB

sys.path.append('../../../utlts/')
from parse_user_input import false_or_true
from read_in_file import read_in
from parse_data import organism
from protein_property import database
from output import writeout 


def get_contacts(d_ref, method):
	d_output = {}
	for oln, name in d_ref.iteritems():
		structure = PDBParser().get_structure('X', "../../../0-identify_structure/2-get_pdb_chain/{0}/{1}.pdb".format(organism, name))

		atom_list = []
		for residue in structure.get_residues():
			if not PDB.is_aa(residue): continue
			
			if method==0:
				atom_list.extend(residue.get_atoms())
			elif method==1:	
				if 'CB' in residue:
					atom_list.append(residue['CB'])
				elif Residue.Residue.get_resname(residue)=='GLY':	#glycine case	#can use is_aa to check	
					if 'CA' in residue:
						atom_list.append(residue['CA'])
					else: pass
			else:
				print 'method index not supported'
				sys.exit()
		if len(atom_list)==0: continue
		ns = NeighborSearch(atom_list)

		radius_list = [4.5, 7.5]
		contact_list = ns.search_all(radius = radius_list[method], level = 'R')	#residues are returned
		contact_list_copy = contact_list[:]
		for contact in contact_list:
			if abs(contact[0].id[1]-contact[1].id[1]) in [1, 0]:
				contact_list_copy.remove(contact)
		d_output[oln] = len(contact_list_copy)
		if len(set(contact_list_copy))!=len(contact_list_copy): print "uh oh"	#seems to be nicely taken care of py NeighborSearch, doesnt output same contact pairs
	return d_output

def get_contact_density(d_contact, d_len):
	d = {}
	for oln, contact in d_contact.iteritems():
		d[oln] = contact/float(d_len[oln]) 
	return d


if __name__ == "__main__":
	extra = ''
	if false_or_true("Relax selection criterion 2"):
		extra += 'pre_output'

	d_ref = read_in('oln', 'pdb', filename = extra)
	d_len = read_in('oln', 'length')

	method = false_or_true("Calculate contact density like Shakh2006 [default Zhou2008]")
	d_contact = get_contacts(d_ref, method)

	d_cd = get_contact_density(d_contact, d_len)
	filename = 'PDB'
	if method:
		filename+='_shakh'
	writeout(['oln', 'contact_density'], d_cd, filename = '{0}{1}'.format(filename, extra))
