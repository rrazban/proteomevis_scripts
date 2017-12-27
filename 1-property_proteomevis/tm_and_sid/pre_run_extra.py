#!/usr/bin/python

help_msg = 'write extra.txt files for ./run --extra/--EXTRA'

import sys, os

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import set_up_read_in, read_in
from parse_data import organism
from output import writeout


if __name__ == "__main__":
	help_message(help_msg)	#need to adjust help message to allow yeast_ecoli case
	if organism=='yeast_ecoli':	#dependent on yeast/extra.txt being present
		d = read_in('uniprot', 'pdb', organism='ecoli')
		d_old = read_in('uniprot', 'pdb', filename='../ecoli/extra.txt')
	else:
		d = read_in('uniprot', 'pdb')
		d_old = read_in('uniprot', 'pdb', filename = '../../../0-identify_structure/3-length_check/{0}/{1}'.format(organism, 'old_seq2struc.txt'))
	pdb_list = set(d.items()) - set(d_old.items())
	d_output = dict(x for x in pdb_list)
	writeout(['uniprot', 'pdb'], d_output, filename='extra')
	print 'proceed to enter ../run --extra'

