#!/usr/bin/python
import os,sys

from parse_user_input import false_or_true


def database(organism, option, bool_more=False):
	DIR = os.path.dirname(os.path.realpath(__file__)) 
	DIR_pv = DIR + '/../1-property_proteomevis'
	DIR_proteome = DIR + '/../property_proteome'
	d_db = {}

	if option in ['tm', 'sid']:
		d_db[option] = [['pdb1', 'pdb2'], option.upper(), "{0}/tm_and_sid/{1}/PDB.txt".format(DIR_pv, organism)]

	elif option=='length':
		d_db[option] = ['oln', 'length', '', '', organism] 

	elif option in ["contact_density", "con_den"]:
		d_db[option] = ['oln', 'contact_density', "{0}/{1}/{2}/PDB.txt".format(DIR_pv, 'contact_density', organism)]

	elif option=="abundance":
		if organism=='ecoli': filename = 'Arike_2012'
		elif organism=='yeast': filename = 'Sina_2003' 
		else: pass
		d_db[option] = ["string_external_id", "abundance", '{0}/abundance/{1}/{2}.txt'.format(DIR_proteome, organism, filename), ('l', '.')]

	elif "ppi" in option.lower():
		filename = 'intact'
		d_db[option] = ['oln', 'ppi', '{0}/PPI_degree/{1}/{2}.txt'.format(DIR_proteome, organism, filename)]	

	elif "dos" in option.lower():
		if organism=='ecoli':
			d_db[option] = ["ensembl", "dt", '{0}/dosage_tolerance/ecoli/Kitagawa_2005.txt'.format(DIR_proteome)]
		elif organism=='yeast':
			d_db[option] = ['Probe ID', 'Log2 Ratio (20Gen)', '{0}/dosage_tolerance/yeast/Douglas_2012.txt'.format(DIR_proteome), ('r', ":")]
		else: pass

	elif "rate" in option.lower():
		d_db[option] = ["oln", "evorate", "{0}/evolutionary_rate/{1}/Zhang_2015.txt".format(DIR_proteome, organism)]

	else:
		"{0},{1} option not found".format(organism, option)

	return d_db[option]


