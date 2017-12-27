#!/usr/bin/python

help_msg = 'generate images of PDB structures in pdb_image/'
custom_usage = 'PYMOLPATH/pymol/pymol.exe -qc get_image.py'

import sys, os, subprocess

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import organism_list, initialize_dict
from output import print_next_step


def pnghack(filepath, width=2000, height=2000):	
    """Workaround if cmd.png() doesn't work"""	#cmd.png() doesnt work with api
    cmd.set('ray_trace_frames', 1)  # Frames are raytraced before saving an image.
    cmd.viewport(width, height)  # Set resolution
    cmd.mpng(filepath, 1, 1)  # Use batch png mode with 1 frame only
    cmd.mplay()  # cmd.mpng needs the animation to 'run'

def get_image(pdb, path, organism):
	pdb_file = '../../0-identify_structure/2-get_pdb_chain/{0}/{1}.pdb'.format(organism, pdb)
	cmd.load(pdb_file, pdb)
	cmd.show_as('cartoon')
	cmd.disable("all")
	cmd.enable(pdb)
	cmd.set('ray_opaque_background', 0)
	cmd.spectrum('count', 'rainbow')
	pnghack("{0}.png".format(path))
	cmd.delete(pdb)

def make_intermediate_dir(dir1, dir2):
	if not os.path.isdir("pdb_image/{0}".format(dir1)):
		os.makedirs("pdb_image/{0}/{1}".format(dir1, dir2))
	elif not os.path.isdir("pdb_image/{0}/{1}".format(dir1, dir2)):
		os.mkdir("pdb_image/{0}/{1}".format(dir1, dir2))
	else:
		pass

def get_path(pdb):
	dir1 = pdb[0]
	dir2 = pdb[1]
	path = "pdb_image/{0}/{1}".format(dir1, dir2)
	if os.path.isdir(path):
		make_intermediate_dir(dir1, dir2)
	return path

def save_image(d):
	for organism, d_pdb in d.iteritems():
		for pdb in d_pdb:	#first check if image already exists
			path = '{0}/{1}.png'.format(get_path(pdb), pdb)
			if os.path.exists(path): continue 	#already exists. easier to be true than false
			get_image(pdb, path, organism)
			subprocess.call(['mv', "{0}0001.png".format(path), "{0}".format(path)])	#seems that appending frame number to output irreversible since it is a movie
												#extra .png because if chain is a number, that info is lost
	cmd.quit()

if __name__ == "pymol":	#not __main__ cuz running through pymol
	d = initialize_dict('dict')
	for organism in organism_list:
		d[organism] = read_in('pdb', 'uniprot', organism=organism)
	save_image(d)
	print_next_step()

if __name__== "__main__":
	args = help_message(help_msg, bool_org_dir=False, custom_usage = custom_usage)
