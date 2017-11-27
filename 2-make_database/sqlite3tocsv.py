#!/usr/bin/python

#'db.sqlite3' -> csv files found in https://github.com/rrazban/proteomevis/make_database

import sys, os
import sqlite3
import csv


column_list = ['proteomevis_species', 'proteomevis_chain', 'proteomevis_domain', 'proteomevis_edge']

with sqlite3.connect("db.sqlite3") as connection:
	for column in column_list:
		csvWriter = csv.writer(open("csv/{0}.csv".format(column), "w"))
		c = connection.cursor()

		c.execute("SELECT * FROM {0}".format(column))
		rows = c.fetchall()
		csvWriter.writerow([d[0] for d in c.description])
		csvWriter.writerows(rows)

