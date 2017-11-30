#!/usr/bin/python

#set up extra.txt files for run -e jobs

import sys, os

sys.path.append('../../../utlts/')
from read_in_file import set_up_read_in, read_in
from parse_data import organism
from output import writeout


if __name__ == "__main__":
	if organism=='yeast_ecoli':	#dependent on yeast/extra.txt being present
		d = read_in('uniprot', 'pdb', organism='ecoli')
		d_old = read_in('uniprot', 'pdb', filename='../ecoli/extra.txt')
	else:
		d = read_in('uniprot', 'pdb')
		d_old = read_in('uniprot', 'pdb', filename = '../../../0-identify_structure/3-length_check/{0}/{1}'.format(organism, 'old_seq2struc.txt'))
	pdb_list = set(d.items()) - set(d_old.items())
	d_output = dict(x for x in pdb_list)
	writeout(['uniprot', 'pdb'], d_output, filename='extra')


