
import os

# Represents genes which are mutually conserved between two organisms.
# Two input files comprise this gene-set: a query and target. The query
# is such that the first column is the query organism gene, and its second
# column is its BLASTP homology. Similarly, the target sequence is the
# reciprocal. For a gene to be mutually conserved, its homology must be
# the same given the query and target.  
class MutuallyConservedGeneSet():
	def __init__(self, query, target):
		self.fname_query = query # query filename
		self.fname_target = target # target filename
		self._is_query_valid_file() # test query file is valid
		self._is_target_valid_file() # test target file is valid
	
	# Determine if the query filename is valid or not
	def _is_query_valid_file(self):
		if os.path.exists(self.fname_query):
			return True
		else:
			raise OSError(self.fname_query +' is not a valid file [ERROR]')
	
	# Determine if the target filename is valid or not
	def _is_target_valid_file(self):
		if os.path.exists(self.fname_target):
			return True
		else:
			raise OSError(self.fname_target +' is not a valid file [ERROR]')

# Get the number of species. Each has a query and target file to it.
def get_num_species():
	num_species = raw_input('How many species to analyze? (incl. control): ')
	if num_species.isdigit(): # only integers are accepted
		return int(num_species)
	else:
		raise ValueError('#/species must be an integer [ERROR]')

if __name__ == '__main__':
	print('--- Identify conserved genes across multi-species BLASTP hits ---')
	try:
		n = get_num_species() # get number of analyses
		for i in range(n): # per analysis number, create a conserved gene-set
			print('\n--- Analysis #' + str(i+1) + ' ---')
			query = raw_input('Query BLASTP filename: ')
			target = raw_input('Target BLASTP filename: ')
			geneset = MutuallyConservedGeneSet(query, target) # wrap both files
	except ValueError as e:
		print(e)
	except OSError as e:
		print(e)	