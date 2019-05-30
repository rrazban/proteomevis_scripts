#!/usr/bin/python

help_msg = 'make 1D histogram'

import sys, os
import numpy as np
from scipy.stats import spearmanr, hmean
import matplotlib.pyplot as plt

from versus import parse_input, merge

sys.path.append('../utlts/')
from read_in_file import read_in
from parse_user_input import false_or_true, which_organism
from properties import database
from label_plot import get_label, get_title, organism_color 


def plotout(data, tot_num, title, xlabel):
	color = organism_color()
	num = len(data)
	data = np.array([float(datum) for datum in data])

	plt.hist(data, color = color[organism], density = True, label="n={0} ({1:.0f}%)".format(num, 100*num/tot_num))
	plt.title(title), plt.xlabel(xlabel), plt.ylabel("probability"), plt.legend()
	plt.show()



if __name__ == "__main__":
	organism = which_organism()
	d_x, xlabel, proteome_subset_bool = parse_input(str(raw_input("Property x: ")), organism)
	if proteome_subset_bool:
		proteome_subset_bool = false_or_true("Include proteins not in ProteomeVis")

	data, _ = merge(d_x, d_x, proteome_subset_bool, organism)

	tot_num =  len(read_in('Entry', 'Entry', 'proteome', organism=organism))
	plotout(data[0], tot_num, get_title(organism, proteome_subset_bool), xlabel)	#TM or SID case cant be determined
