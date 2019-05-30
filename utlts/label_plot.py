#!/usr/bin/python
import sys	


def get_label(dir_input, option=''):
	if dir_input=='proteome':
		return "length (UniProt)"
	elif dir_input=='':
		return
	else: pass
	dir_input_list = dir_input.split('/')
	print dir_input
	quantity = dir_input_list[-3]
	source = dir_input_list[-1][:-4]
	if quantity == "tm.sid":
		quantity = option.upper()
	elif quantity=="dosage_tolerance" and False:	#also need bin==True
		if option==True: quantity = "ratio of dosage sensitive genes"	#have not yet implemented
	elif quantity=="ppi":
		quantity = "PPI_degree"
	else: pass
	return "{0} ({1})".format(quantity, source)

def get_title(organism, proteome_subset_bool):
	if proteome_subset_bool:
		data_amount = 'proteome'
	else:
		data_amount = 'ProteomeVis'
	d = organism_latin()
	return "{0} ({1})".format(d[organism], data_amount)

def organism_latin():
	d={}
	d['yeast'] = r'$\mathit{S.cerevisiae}$'
	d['ecoli'] = r'$\mathit{E.coli}$'
	d['human'] = r'$\mathit{H.sapiens}$'
	d['protherm'] = "ProTherm"

	d['yeast_literature'] = r'$\mathit{S. cerevisiae}$ literature'
	d['ecoli_literature'] = r'$\mathit{E. coli}$ literature'
	return d

def organism_color():	
	color = {}
	color['ecoli'] = 'blue'
	color['yeast'] = 'orange'
	color['yeast_ecoli'] = 'green'
	color['human'] = 'purple'
#	color['protherm'] = 'purple'

	color['ecoli_literature'] = 'purple'
	color['yeast_literature'] = 'red'
	return color

def organism_hatch():
	d = {}
	d['ecoli'] = '+'	#literature -
	d['yeast'] = 'x'	#literature \

	d['ecoli_literature'] = '-'
	d['yeast_literature'] = '\\'
	return d

def TM_label():	
	d = organism_latin()
	label = {}
	label['ecoli'] = '({0}, {0})'.format(d['ecoli'])
	label['yeast'] = '({0}, {0})'.format(d['yeast'])
	label['yeast_ecoli'] = '({0}, {1})'.format(d['yeast'], d['ecoli'])
	return label 

def full_organism_name():	#taken directly from UniProt
	d = {}
	d['yeast'] = 'Saccharomyces cerevisiae (strain ATCC 204508 / S288c)' 
	d['ecoli'] = 'Escherichia coli (strain K12)'
	return d
