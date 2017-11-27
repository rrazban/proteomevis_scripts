#!/usr/bin/python

#set up extra.txt files for run -e jobs

import sys, os

sys.path.append('../../../utlts/')
from read_in_file import set_up_read_in, read_in
from parse_data import organism
from output import writeout


if __name__ == "__main__":
	if organism=='yeast_ecoli':	#dependent on yeast/extra.txt being present
		d = read_in('pdb', 'uniprot', organism='ecoli')
		d_old = read_in('pdb', 'uniprot', filename='../ecoli/extra.txt')
	else:
		d = read_in('pdb', 'uniprot')
		d_old = read_in('pdb', 'uniprot', filename = '../../../0-identify_structure/3-length_check/{0}/output/{1}'.format(organism, 'old_seq2struc.txt'))
	pdb_list = set(d.items()) - set(d_old.items())
	d_output = dict(x for x in pdb_list)
	writeout(['pdb', 'uniprot'], d_output, filename='extra.txt')

