#!/usr/bin/python

help_msg = 'obtain ppi degree for proteome'

import os, sys
import imp
from bioservices import PSICQUIC

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import taxid, organism
from output import writeout


if __name__ == "__main__":
	help_message(help_msg)

	module = imp.load_source("run", "../../../1-property_proteomevis/ppi_partner/run.py") #normal import doesnt work, proly cuz name of module the same as this script 
	d_ppi, error_list, filename = module.get_physical_ppi(partner_bool=False)

	writeout(['oln','PPI_degree'], d_ppi, filename=filename)
	print "Error list: {0} ({1})".format(error_list, len(error_list))
