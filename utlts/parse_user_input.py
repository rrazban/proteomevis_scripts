#!/usr/bin/python
import sys, os


def false_or_true(message):	#False is default
	if (raw_input("{0} (False/True)? ".format(message))).strip() in ['1', 'True', 'T']:
		control = 1
	else:	
		control = 0
	return control

def which_organism(yeast_ecoli_bool=False):
	if yeast_ecoli_bool == True:
		message = "Which organism's protein pairs (yeast/ecoli/yeast_ecoli)? "
	else:
		message = "Which organism's proteome (yeast/ecoli)? "

	pre_organism = raw_input(message)
	pre_organism = pre_organism.strip().lower()
	if pre_organism in ['yeast', '0', 'y']:	
		organism = 'yeast'
	elif pre_organism in ['ecoli', '1', 'e']:
		organism = 'ecoli'
	elif pre_organism in ['yeast_ecoli', '2', 'ye']:
		organism = 'yeast_ecoli'
	else:
		print "sorry, {0} organism not supported".format(pre_organism)
		sys.exit()	
	return organism
