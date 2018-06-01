import sys, os


def setup(filename, organism=''):	#only tab seperated files accepted
	if '.txt' not in filename:
		DIR = os.path.dirname(os.path.realpath(__file__)) 
		DIR += '/../0-identify_structure/'
		if not organism:
			from parse_data import organism

		if filename=='':
			filename = "{0}/3-length_check/{1}/seq2struc.txt".format(DIR, organism)
		else:	#proteome or pre_seq2struc
			filename = "{0}/0-identify_pdb/{1}/{2}.txt".format(DIR, organism, filename)
	return filename

def introduction(line):
	if line[0]!='#':
		labels = line.split('\t')
		labels = [x.rstrip() for x in labels]
		return labels

def duplicate_check(filename, value, d):
	continue_bool = False
	if value in d:
		print '{0} - duplicate found: {1}'.format(os.path.basename(filename), value)
		continue_bool = True
	return continue_bool

class ReadIn():
	def __init__(self, filename = '', organism = ''):
		self.labels = []
		self.d_output = {}
		self.filename =	setup(filename, organism) 		

	
	def run(self, xlabel, ylabel, strip = ''):
		with open(self.filename, 'r') as rfile:
			for line in rfile:
				if not self.labels:
					self.labels = introduction(line)
				else:
					if not xlabel:
						xlabel = self.labels[0]

					word_list = line.split('\t')
					value_list = []
					for w,want in enumerate([xlabel, ylabel]):
						word_list = [word.strip() for word in word_list] #extra strip to get rid of \n since defining delineator
						if w==0 and type(want)==list:	#for tm/sid read ins
							want_list = []
							for want_item in want:
								want_list.append(word_list[self.labels.index(want_item)])
							value = ','.join(want_list)
						else:
							value = word_list[self.labels.index(want)]
							if w==0 and strip: 
								if strip[0]=='l': value = value[value.index(strip[1])+1:]	#strip() overstrips sometimes
								elif strip[0]=='r': value = value[:value.index(strip[1])]	

						value_list.append(value)

					if duplicate_check(self.filename, value_list[0], self.d_output): continue
					self.d_output[value_list[0]] = value_list[1]
		return self.d_output

def read_in(xlabel, ylabel, filename='', strip='', organism=''):
	readin = ReadIn(filename, organism)
	return readin.run(xlabel, ylabel, strip)
