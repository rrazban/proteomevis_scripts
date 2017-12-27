#!/usr/bin/python

help_msg = 'remove extra PDB images in pdb_image/'

import sys, os, glob

from get_image import get_path

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import organism_list, initialize_dict


def get_file():
	f = []
	for (dirpath, dirnames, filenames) in os.walk("pdb_image/"):
		f.extend(filenames)
	f.remove('.gitkeep')
	return f

def update_file_list(file_list, d):
	update = file_list[:]
	for organism, d_pdb in d.iteritems():
		for pdb in d_pdb:
			pdb_file = pdb+'.png'
			if pdb_file in file_list:
				update.remove(pdb_file)
	return update

def remove_image(update):
	for pdb_file in update:
		path = '{0}/{1}'.format(get_path(pdb_file), pdb_file)
		print path
		os.remove(path)

if __name__ == "__main__":
	help_message(help_msg, bool_org_dir=False)
	file_list = get_file()
	d = initialize_dict('dict')
	for organism in organism_list:
		d[organism] = read_in('pdb', 'uniprot', organism=organism)
	update = update_file_list(file_list, d)
	remove_image(update)
