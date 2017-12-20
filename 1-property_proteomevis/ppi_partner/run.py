#!/usr/bin/python

help_msg = 'find which protein pairs in seq2struc form PPIs'

import os, sys
import collections
from bioservices import PSICQUIC

CWD = os.getcwd()
UTLTS_DIR = CWD[:CWD.index('proteomevis_scripts')]+'/proteomevis_scripts/utlts'
sys.path.append(UTLTS_DIR)
from parse_user_input import help_message
from read_in_file import read_in
from parse_data import organism, taxid
from output import writeout, print_next_step 


def get_score(ppis, score_crit):
	new_ppis = ppis
	for ppi in ppis:
		score = ppi[-1][len('intact-miscore:'):]	#by construction, always is last element
		try: 			
			score = float(score)
		except:
			score = 0
		if score<score_crit:
			new_ppis.remove(ppi)			
	return new_ppis

def get_ppi_partner(uniprot, pdb, ppis, d_ref):
	d = collections.defaultdict(list)	
	for ppi in ppis:			#make sparse list, only have those with ppi interaction. everyone else is false
		want = u'uniprotkb'
		for ppi_info in ppi:
			if want in ppi_info:
				uniprot_partner = ppi_info[len(want)+1:].encode('utf-8')
				if uniprot_partner != uniprot:
					if uniprot_partner in d_ref: 
						d[pdb].append(d_ref[uniprot_partner])
				break	#usually second item in ppi_info
	return d

def get_ppi_degree(ppis):
	if ppis == [[u'']]:
		num_ppi = 0
	else:
		num_ppi = len(ppis)
	return num_ppi

def get_physical_ppi(partner_bool=True):
	if partner_bool:
		d_ref = read_in('uniprot', 'oln')
		d_ref2 = read_in('uniprot', 'pdb')
	else:
		d_ref = read_in('Entry', 'Gene names  (ordered locus )', 'proteome')
	taxonomy = taxid()[organism]
	db = 'intact'		#oln not supported for ecoli for mint and biogrid
	score_crit = None

	d = {}
	error_list = []
	s = PSICQUIC(verbose=False)
	for uniprot, oln in d_ref.iteritems():
		try:
			ppis = s.query(db, "{0} AND taxid:{1} AND affinity".format(uniprot, taxonomy))
		except:
			ppis = []
			error_list.append(uniprot)
			print "error! can't find ppis for {0}".format(uniprot)
			continue

		if score_crit != None:
			ppis = get_score(ppis, score_crit)	
		
		if partner_bool:
			d.update(get_ppi_partner(uniprot, d_ref2[uniprot], ppis, d_ref2))
		else:
			d[oln] = get_ppi_degree(ppis)

	if score_crit:
		db+= '_{0}'.format(score_crit)	
	return d, error_list, db


if __name__ == "__main__":
	help_message(help_msg)
	d_ppi, error_list, filename = get_physical_ppi()	
	writeout(['protein', 'protein partners'], d_ppi, filename='ppi_partners')
	print "Error list: {0} ({1})".format(error_list, len(error_list))
	print_next_step()
