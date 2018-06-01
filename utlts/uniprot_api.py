import urllib, urllib2

from parse_data import taxid

#UniProt column names are found at
#https://www.uniprot.org/help/uniprotkb_column_names

class UniProtAPI():
	def __init__(self, columns):
		self.columns = columns
		self.url = 'https://www.uniprot.org/uniprot/'	
		self.batch_size = 350	#491 is limit

		self.raw_data = []

	def info(self):
		data = urllib.urlencode(self.params)
		request = urllib2.Request(self.url, data)
		response = urllib2.urlopen(request)
		labels = next(response).split('\t')
		self.raw_data.extend(response)
		self.labels = [label.rstrip() for label in labels]

	def uniprot_info(self, uniprots):
		for batch_i in range(len(uniprots) / self.batch_size + 1):
			self.params = {'query':','.join(uniprots[(batch_i)*self.batch_size:(batch_i+1)*self.batch_size]), 'columns':','.join(self.columns), 'format':'tab'}	#as of 5/29/18 cannot get multiple uniprots at once
			self.info()
		return self.labels, self.raw_data
		

	def organism_info(self, organism = ''):
		if not organism:
			from parse_data import organism
		self.params = {'query':'organism:{0} AND reviewed:yes'.format(taxid()[organism]), 'columns':','.join(self.columns), 'format':'tab'}
		self.info()
		return self.labels, self.raw_data
