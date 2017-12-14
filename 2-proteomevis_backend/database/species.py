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

	table_name = 'proteomevis_species'
	conn = sqlite3.connect('db.sqlite3')
	c = conn.cursor()
	c.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
	c.execute('CREATE TABLE {0}(id,name,has_localization)'.format(table_name))

	line_list_list = [[0, d_org[0], 1], [1, d_org[1], 1]]
	for line_list in line_list_list:
		c.execute("INSERT INTO {0} VALUES {1}".format(table_name, tuple(line_list))) 
	conn.commit()
	conn.close()


if __name__ == "__main__":
	help_message(help_msg, bool_org_dir = False)	#add verbose option
	printout_sql()
	print_next_step('../')	
