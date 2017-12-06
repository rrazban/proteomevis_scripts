#!/usr/bin/python

help_msg = "set organism parameters"

import sys, os
import sqlite3

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from parse_data import int2organism
from output import print_next_step


def printout_sql():
	d_org = int2organism()

	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE proteomevis_species')
	c.execute('CREATE TABLE proteomevis_species(id,name,has_localization,has_function,has_mutant_data)')

	line_list_list = [[0, d_org[0], 0, 0, 0], [1, d_org[1], 0, 0, 0]]
	for line_list in line_list_list:
		c.execute("INSERT INTO proteomevis_species VALUES {0}".format(tuple(line_list))) 
	conn.commit()
	conn.close()


if __name__ == "__main__":
	help_message(help_msg, bool_org_dir = False)	#add verbose option
	printout_sql()
	print_next_step('../')	
