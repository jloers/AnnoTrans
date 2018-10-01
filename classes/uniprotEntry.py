#!/usr/bin/env python3.6

__author__  = "Jens Loers"
__email__   = "jens.loers@uni-bielefeld.de"

class UniprotEntry:
	"""
	Datastructure to store all information of interest from a Uniprot
	Entry. Also allows basic manipulation of the data like asigning
	additional identifier.
	"""
	
	def __init__(self, ID, aa):
		"""
		inititaliza and store all information of interest in a uniprot 
		entry. Might be extendend over developemental progress
		"""
		self.ID						= ID	# store primary ID
		self.secondaryIDs			= []	# secondary IDs e.g. use for mapping
		self.amino_acids			= aa	# contains core amnio acid range of uniprot annotation
		self.functions				= []	# stores all functions with regard to amino acid sites

	def createIsoformAnotation(self):
		"""
		adjust the nomenclature for an annotated isoform
		"""
		print('TODO')
		
