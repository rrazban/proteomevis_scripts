#!/usr/bin/python

help_msg = "'db.sqlite3' to csv files in csv/"

import sys, os
import sqlite3
import csv

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from output import print_next_step


def sqlite3tocsv():
	column_list = ['proteomevis_species', 'proteomevis_chain', 'proteomevis_inspect', 'proteomevis_edge']

	with sqlite3.connect("db.sqlite3") as connection:
		for column in column_list:
			csvWriter = csv.writer(open("csv/{0}.csv".format(column), "w"))
			c = connection.cursor()

			c.execute("SELECT * FROM {0}".format(column))
			rows = c.fetchall()
			csvWriter.writerow([d[0] for d in c.description])
			csvWriter.writerows(rows)

if __name__ == "__main__":
	help_message(help_msg, bool_org_dir = False)
	sqlite3tocsv()
	print_next_step('../')
