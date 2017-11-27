#!/usr/bin/python

#make figure in Edge Filtering panel for all organisms

import sys, os
import matplotlib.pyplot as plt

sys.path.append('../../utlts/')
from read_in_file import read_in
from parse_data import organism_list 
from protein_property import database


def merge(d1, d2, data_ppi):
	data = []
	color = []
	for key in d1.iterkeys():
		data.append(tuple(d[key] for d in [d1, d2]))

		if key in data_ppi:
			color.append('g')
		else:
			color.append('k')
	return data, color


def read_in_ppi_partners():	
	data = []
	with open("../../1-protein_property/ppi_partner/{0}/ppi_partners.txt".format(organism), "r") as rfile:
		label_list = next(rfile).split('\t')
		label_list = [x.rstrip() for x in label_list]
		for line in rfile:
			word_list = line.split()
			protein = word_list[label_list.index('protein')]	
			ppi = word_list[label_list.index('protein partners'):]
			for protein_partner in ppi:
				name = protein+','+protein_partner
				data.append(name)
	return data	

def plotout(data_list, color, d_label):
	fig = plt.figure()
	ax=fig.add_axes([0,0,1,1])
	ax.set_axis_off()
	ax.scatter(data_list[1], data_list[0], c=color, s=1)
	ax.set_xlim([0,1])
	ax.set_ylim([0,1])
	ax.set_aspect('auto')
	plt.savefig('species.{0}.png'.format(int(d_label[organism])), transparent="True")
#	plt.show()


if __name__ == "__main__":
	d_label = read_in('name', 'id', filename='../../2-make_database/csv/proteomevis_species.csv')
	protein_property_list = ['sid', 'tm']
	for organism in organism_list:
		predata_list = []
		for protein_property in protein_property_list:
			x_input = database(organism, protein_property)
			d = read_in(*x_input)
			predata_list.append(d)
		data_ppi = read_in_ppi_partners()	
		data_tup, color = merge(predata_list[0], predata_list[1], data_ppi)
		data_list = zip(*data_tup)
		plotout(data_list, color, d_label)
