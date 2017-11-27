#!/usr/bin/python
import sys
from collections import defaultdict
from bioservices import PSICQUIC

sys.path.append('../../../utlts/')

from read_in_file import read_in
from parse_data import organism, taxid
from output import writeout


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

def get_physical_ppi(taxonomy, d_ref, db, score_crit):
	d = defaultdict(list)
	error_list = []
	s = PSICQUIC(verbose=False)

	for uniprot, pdb in d_ref.iteritems():
		try:
			ppis = s.query(db, "{0} AND taxid:{1} AND affinity".format(uniprot, taxonomy))
		except:
			ppis = []
			error_list.append(uniprot)
			print "error! can't find ppis for {0}".format(uniprot)
			continue

		if score_crit != None:
			ppis = get_score(ppis, score_crit)	

		for ppi in ppis:			#make sparse list, only have those with ppi interaction. everyone else is false
			want = u'uniprotkb'
			for ppi_info in ppi:
				if want in ppi_info:
					uniprot_partner = ppi_info[len(want)+1:].encode('utf-8')
					if uniprot_partner != uniprot:
						if uniprot_partner in d_ref: 
							d[pdb].append(d_ref[uniprot_partner])
						break	#usually it is the second iteration
	return d, error_list


if __name__ == "__main__":
	d_ref = read_in('uniprot', 'pdb')

	d_taxid = taxid()
	db = 'intact'		#ensembl not supported for ecoli for mint and biogrid
	score_crit = None
	d_ppi, error_list = get_physical_ppi(d_taxid[organism], d_ref, db, score_crit)	#get ppi for every gene
	writeout(['protein', 'protein partners'], d_ppi, filename='ppi_partners')
	print "Error list: {0} ({1})".format(error_list, len(error_list))
