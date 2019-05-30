#!/usr/bin/python

help_msg = 'make 2D scatterplots'

import sys, os
import numpy as np
from scipy.stats import spearmanr, linregress 
import operator
import matplotlib.pyplot as plt

sys.path.append('../utlts/')
from parse_user_input import which_organism, false_or_true
from read_in_file import read_in
from properties import database
from label_plot import get_title, organism_color



def get_label(user_input, input_dir):
	return "{0} ({1})".format(user_input, os.path.basename(input_dir)[:-4])

def get_operation(user_input, organism):
	if user_input[0] == '-':	#only first character can by -
		sign = operator.neg
		user_input = user_input[1:]
	elif user_input[0] == '|':
		sign = operator.abs
		user_input = user_input[1:-1]
	else:
		sign = operator.pos

	if '/' in user_input:
		divide_i = user_input.index('/')
		z_inputs = (user_input[:divide_i], user_input[divide_i+1:])
		operation = operator.div
	elif '*' in user_input:
		multiply_i = user_input.index('*')
		z_inputs = (user_input[:multiply_i], user_input[multiply_i+1:])
		operation = operator.mul
	else:
		z_inputs = [user_input]
		operation = ''		
	return [database(organism, z_input)  for z_input in z_inputs], sign, operation

def reconcile(d1, d2, sign, operation=''):
	d_input = {}
	for key, value in d1.iteritems():
		if key in d2:
			try:
				if operation:
					d_input[key] = sign(operation(float(value), float(d2[key])))
				else:
					d_input[key] = sign(float(value))
			except:
				pass
	return d_input

def parse_input(user_input, organism):
	inputs, sign, operation = get_operation(user_input, organism)
	if len(inputs)>1:
		d = reconcile(read_in(*inputs[0]), read_in(*inputs[1]), sign, operation)
	else:
		d = reconcile(read_in(*inputs[0]), read_in(*inputs[0]), sign)
	return d, get_label(user_input, inputs[0][2]), 'PDB' not in inputs[0][2]

def set_condition(d_ref, d_condition, min_, max_):
	total = 0
	d = {}
	for oln, condition in d_condition.iteritems():
		if oln in d_ref:
			total+=1
			if min_ < condition < max_:
				d[oln] = d_ref[oln]
	print '% of original dataset present: {0}'.format(len(d)/float(total))
	return d

def parse_condition(d_ref, organism):
	while True:
		condition = str(raw_input("Set min/max condition: "))
		if condition:
			d_condition, _, _ = parse_input(condition, organism)
			min_ = float(raw_input("Min: "))
			max_ = float(raw_input("Max: "))
			d_ref = set_condition(d_ref, d_condition, min_, max_)
		else:
			break
	return d_ref

def reference(proteome_subset_bool, organism):
	if organism=='protherm':
		d_ref = read_in('uniprot', 'pdb', organism=organism)
	elif proteome_subset_bool:
		d_ref = read_in('Gene names  (ordered locus )', 'Entry', 'proteome', organism=organism)	#based on oln cuz proteome data usually have oln only
	else:
		d_ref = read_in('oln', 'pdb', organism=organism)
	return  parse_condition(d_ref, organism)

def merge(d1, d2, proteome_subset_bool, organism):
	d_ref = reference(proteome_subset_bool, organism)

	data = []
	labels = []
	for key in d1.iterkeys():
		if key not in d_ref:
			continue	
		if key in d2:
			data.append(tuple(float(d[key]) for d in [d1, d2]))
			labels.append(key)
	return zip(*data), labels


def plotout(organism, data, num_list, title, label_list, labels):
	color = organism_color()

	r, pvalue = spearmanr(data[0], data[1])
	plt.plot(data[0], data[1], 'bo', c = color[organism], label="$\\rho=${0:.2f} ({1:.2E})\n$n=${2} ({3:.0f}%)".format(r, pvalue, num_list[0], 100*num_list[0]/num_list[1]))

	plt.title(title)
	plt.xlabel(label_list[0]), plt.ylabel(label_list[1]), plt.legend()
	plt.show()


if __name__ == "__main__":
	organism = which_organism() 
	d_x, xlabel, proteome_subset_bool_x = parse_input(str(raw_input("Property x: ")), organism)
	d_y, ylabel, proteome_subset_bool_y = parse_input(str(raw_input("Property y: ")), organism)

	proteome_subset_bool = proteome_subset_bool_x and proteome_subset_bool_y
	if proteome_subset_bool:
		proteome_subset_bool = false_or_true("Include proteins not in ProteomeVis")

	data, labels = merge(d_x, d_y, proteome_subset_bool, organism)
	num = len(data[0])
	total = len(read_in('Entry','Entry', 'proteome', organism=organism))

	plotout(organism, data, [num, total], get_title(organism, proteome_subset_bool), [xlabel, ylabel], labels)
