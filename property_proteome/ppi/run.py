#!/usr/bin/python
import sys
from bioservices import PSICQUIC

sys.path.append('../../../utlts/')
from read_in_file import read_in
from parse_data import taxid, organism
from output import writeout


def get_score(ppis, score_crit):
	new_ppis = ppis
	for ppi in ppis:
		score = ppi[-1][len('intact-miscore:'):]	#by construction, always is last element
		try: 			
			score = float(score)
		except:
			score = 0
		if score<score_crit:	#threshold 0.2 chosen because where probability diminishes
			new_ppis.remove(ppi)			
	return score, new_ppis

def get_physical_ppi(taxid, ensembl_list, db, score_crit):
	d = {}
	data = []	#make into default dict? append #no just use a unique identifier
	error_list = []
	s = PSICQUIC(verbose=False)

	for ensembl in ensembl_list:
		try:
			ppis = s.query(db, "{0} AND taxid:{1} AND affinity".format(ensembl, taxid))
			if score_crit != None:
				score, ppis = get_score(ppis, score_crit)	#put a bool here
		except:
			ppis = []
			error_list.append(ensembl)
			print "error! can't find ppis for {0}".format(ensembl)
			continue
		if ppis == [[u'']]:
			num_ppi = 0
		else:
			num_ppi = len(ppis)
		d[ensembl] = num_ppi
	return d, error_list


if __name__ == "__main__":
	d_ref = read_in('Gene names  (ordered locus )', 'Entry', filename = 'proteome')
	ensembl_list = d_ref.keys() 

	d_taxid = taxid()
	db = 'intact'		#ensembl not supported for ecoli for mint and biogrid
	score_crit = None
	d_ppi, error_list = get_physical_ppi(d_taxid[organism], ensembl_list, db, score_crit)	#get ppi for every gene

	if score_crit != None:
		db+='_{0}'.format(score_crit)
	writeout(['oln','ppi'], d_ppi, filename=db)
	print "Error list: {0} ({1})".format(error_list, len(error_list))
