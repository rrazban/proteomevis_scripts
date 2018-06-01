import sys, os
import argparse

from parse_data import organism, organism_list


def run_message(bool_add_verbose, bool_org_dir, custom_usage):
	if bool_org_dir:
		usage = '../run.py [-h]'
	elif custom_usage:	#why not make into bool?
		usage = 'PYMOLPATH/pymol/pymol.exe -qc get_image.py'
	else:
		usage = '%(prog)s [-h]'

	if bool_add_verbose:
		usage += ' [-v]'
	return usage 

def help_message(message, bool_add_verbose=False, bool_org_dir=True, custom_usage = ''):
	parser = argparse.ArgumentParser(description = message, usage=run_message(bool_add_verbose, bool_org_dir, custom_usage))
	if bool_add_verbose:
		parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true")

	if bool_org_dir:
		if organism not in organism_list and organism!='yeast_ecoli':
			print 'Make sure to run in valid organism directory'
			parser.print_help()			
			sys.exit()
	if custom_usage:
		parser.print_help()			
		sys.exit()
	return parser.parse_args()
	
def false_or_true(message):	#False is default
	if (raw_input("{0} (False/True)? ".format(message))).strip() in ['1', 'True', 'T']:
		control = 1
	else:	
		control = 0
	return control

def which_organism(yeast_ecoli_bool=False):
	if yeast_ecoli_bool == True:
		message = "Which organism's protein pairs (yeast/ecoli/yeast_ecoli)? "
	else:
		message = "Which organism's proteome (yeast/ecoli)? "

	pre_organism = raw_input(message)
	pre_organism = pre_organism.strip().lower()
	if pre_organism in ['yeast', '0', 'y']:	
		organism = 'yeast'
	elif pre_organism in ['ecoli', '1', 'e']:
		organism = 'ecoli'
	elif pre_organism in ['protherm', '2', 'pro']:
	#	organism = 'yeast_ecoli'
		organism = 'protherm'
	else:
		print "sorry, {0} organism not supported".format(pre_organism)
		sys.exit()	
	return organism
