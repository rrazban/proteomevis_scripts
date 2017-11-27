#!/usr/bin/python
import sys, os


def introduction(line, delineator):
	if line[0]!='#':
		label_list = line.split(delineator)
		label_list = [x.rstrip() for x in label_list]
		return label_list

def set_up_read_in(filename, organism=''):
	delineator = '\t'
	if '.txt' not in filename and '.csv' not in filename and '.tsv' not in filename:
		DIR = os.path.dirname(os.path.realpath(__file__)) 
		DIR += '/../0-identify_structure/'
		if not organism:
			from parse_data import organism

		if filename=='':
			filename = "{0}/3-length_check/{1}/seq2struc.txt".format(DIR, organism)
		else:
			filename = "{0}/0-identify_pdb/{1}/{2}.txt".format(DIR, organism, filename)
	else:
		filename = filename
		if '.csv' in filename:
			delineator = ','
	return filename, delineator

def duplicate_check(filename, value, d):
	continue_bool = False
	if value in d:
		print '{0} - duplicate found: {1}'.format(os.path.basename(filename), value)
		continue_bool = True
	return continue_bool

def read_in(want_x, want_y, filename='', strip='', organism=''):	#maybe just pass organism here for default cases
	filename, delineator = set_up_read_in(filename, organism)

	d = {}
	label_list = []
	with open(filename, 'r') as rfile:
		for line in rfile:
			if not label_list:
				label_list = introduction(line, delineator)
				continue
			word_list = line.split(delineator)
			value_list = []
			for w,want in enumerate([want_x, want_y]):
				word_list = [word.strip() for word in word_list] #extra strip to get rid of \n since defining delineator
				if type(want)==list:
					want_list = []
					for want_item in want:
						want_list.append(word_list[label_list.index(want_item)])
					value_list.append(",".join(want_list))
				elif ";" in word_list[label_list.index(want)] and w==0:	#OLN has to be x axis to be fully accounted for seperately. if y axis, left as list 
					genes = word_list[label_list.index(want)]	#dont want to go here for hist formation cuz overcount, although not by much
					gene_list = genes.split('; ')
					value_list.append(gene_list)
				else:
					value = word_list[label_list.index(want)]
					if w==0 and strip: 
						if strip[0]=='l': value = value[value.index(strip[1])+1:]	#strip overstrips sometimes
						else: value = value[:value.index(strip[1])]	
					value_list.append(value)

			if type(value_list[0]) == list:
				for gene in value_list[0]:
					if duplicate_check(filename, gene, d): continue
					d[gene] = value_list[1]
			else:
				if duplicate_check(filename, value_list[0], d): continue
				d[value_list[0]] = value_list[1]
	return d
