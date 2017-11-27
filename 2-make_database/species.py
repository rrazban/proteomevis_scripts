#!/usr/bin/python
import sys, os
import sqlite3

sys.path.append('../utlts/')
from parse_data import int2organism


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
	printout_sql()
