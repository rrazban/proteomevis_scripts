import os, glob


DIR = os.path.dirname(os.path.realpath(__file__)) 
proteomevis_DIR = DIR + '/../1-property_proteomevis'
proteome_DIR = DIR + '/../property_proteome'

def abundance(organism):
	if organism=='ecoli': filename = 'Arike_2012'
	elif organism=='yeast': filename = 'Sina_2003' 
	else: pass
	return ["string_external_id", "abundance", '{0}/abundance/{1}/{2}.txt'.format(proteome_DIR, organism, filename), ('l', '.')]	#identifier second in list

def dosage_tolerance(organism):
	if organism=='ecoli':
		return ["ensembl", "dt", '{0}/dosage_tolerance/ecoli/Kitagawa_2005.txt'.format(proteome_DIR)]
	elif organism=='yeast':
		return ['Probe ID', 'Log2 Ratio (20Gen)', '{0}/dosage_tolerance/yeast/Douglas_2012.txt'.format(proteome_DIR), ('r', ":")]
	
def dNdS(option, organism):
	if organism=='yeast':
		x_term = 'ORF'
		filename = 'Wall_2005.txt'
	elif organism=='ecoli':
		x_term = 'ecoli'
		filename = 'Dasmeh_2017.txt'
	return [x_term, option, "{0}/evolutionary_rate/{1}/{2}".format(proteome_DIR, organism, filename)]

def TM_align(option, organism):
	return [['pdb1', 'pdb2'], option, "{0}/tm_and_sid/{1}/PDB.txt".format(proteomevis_DIR, organism)]

def database(organism, option):					#make into class
	special_cases = ['abundance', 'dosage_tolerance']	#or make a dict such that strip value in place
	#extra special cases = ['dN', 'dS', 'dN/dS']
	if option in special_cases:
		return globals()[option](organism)
	elif 'dN' in option or 'dS' in option:
		return dNdS(option, organism)		
	elif option in ['TM', 'SID', 'nal', 'align1', 'align2']:
		return TM_align(option, organism)
	
	for sub_DIR in [proteomevis_DIR, proteome_DIR]:
		path = "{0}/{1}".format(sub_DIR, option)
		if os.path.exists(path):
			break
	filenames = glob.glob("{0}/{1}/*.txt".format(path, organism))
	if len(filenames)>1:	#try filenames to see if option is there
		files = [os.path.basename(filename) for filename in filenames]
		file_i = int(raw_input("Which file {0}, by index? ".format(files)))
	else:
		file_i = 0
	print filenames[file_i]	
	return ['', option, filenames[file_i]]
