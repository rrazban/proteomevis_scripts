#!/usr/bin/python
from __future__ import absolute_import

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

	sys.path.append('../../../1-property_proteomevis/ppi_partner/')
	module = imp.load_source("run", "../../../1-property_proteomevis/ppi_partner/run.py")
	d_ppi, error_list, filename = module.get_physical_ppi(partner_bool=False)

	writeout(['oln','ppi'], d_ppi, filename=filename)	#maybe have module
	print "Error list: {0} ({1})".format(error_list, len(error_list))
