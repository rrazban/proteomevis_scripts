#!/usr/bin/python

help_msg = 'get uniprot length of entire proteome'

import os, sys

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in 
from parse_data import organism
from uniprot_api import UniProtAPI
from output import writeout


def parse_chain_length(words, i, verbose):	#put this in class
	if len(words)==1:	#does not capture UniProt peptide case
		if verbose:
			print 'No chain found: {0}. Structure is discarded'.format(words)
		length = '' 
	elif '>' in words[i+1]:
		length = '' 
	elif '?' in words[i+1]:
		length = ''
	elif '?' in words[i] or '<' in words[i]:
		if verbose:
			print 'No starting residue for chain: {0}'.format(words)
		length = int(words[i+1])
	else:	
		length = int(words[i+1]) - int(words[i]) + 1
	return length

class UniProtLength():
	def __init__(self, verbose, d_ref):
		self.verbose = verbose
		self.d_ref = d_ref

		uniprotapi = UniProtAPI(['id', 'feature(CHAIN)'])
		if organism=='new_protherm':
			print len(d_ref)
			self.labels, self.raw_data = uniprotapi.uniprot_info(d_ref.keys())
		else:
			self.labels, self.raw_data = uniprotapi.organism_info()

		self.d_output = {}

	def run(self):
		for line in self.raw_data:
			words = line.split()
			uniprot = words[self.labels.index('Entry')]
			if uniprot in self.d_ref:
				chain_length_i = self.labels.index('Chain')+1
				chain_length = parse_chain_length(words, chain_length_i, self.verbose)
				if chain_length:
					self.d_output[uniprot] = chain_length
		return self.d_output


if __name__ == "__main__":
	args = help_message(help_msg, bool_add_verbose = True)
        d_ref = read_in('Entry', 'Gene names  (ordered locus )', filename = 'proteome')

	uniprot_length = UniProtLength(args.verbose, d_ref)
	d_output = uniprot_length.run()

	if organism!='protherm':
		d_output = {d_ref[uniprot]: res for uniprot, res in d_output.iteritems()}
		xlabel = 'oln'
	else:	#not supported for ProTherm
		xlabel = 'uniprot'
	writeout([xlabel, 'length'], d_output, filename = 'UniProt')
