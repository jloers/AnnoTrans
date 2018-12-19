### Jens Loers ###
### v0.1 ###
### jens.loers@uni-bielefeld.de ###

from sys import argv
import subprocess, os
import classes.uniprotEntry as uE

def load_multiple_fasta_file(multiple_fasta_file, mode):
	""" loads multiple fasta file into dict """
	
	content = {}
	
	with open( multiple_fasta_file ) as f:
		if mode == 'annotation':
			header = f.readline().strip()[1:].split(' ')[0]
		else:
			header = f.readline().strip()[1:].split('|')[1]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				if mode == 'annotation':
					header = line[1:].split(' ')[0]
				else:
					header = line[1:].split('|')[1]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	
	return content

def loadUniprotGFF(path_to_file):
	"""
	Function to load uniprot information in GFF Format
	"""
	allEntrys = {}
	with open(path_to_file, "r") as DataFile:
		entrys = DataFile.read().split('##sequence-region ')
		for entry in entrys[1:]:
				
			# evaluate first line of entry
			lines = entry.split('\n')
			uniID = lines[0].split(' ')[0]
			start = lines[0].split(' ')[1]
			end   = lines[0].split(' ')[2]

			## create a new uniprot object
			newEntry = uE.UniprotEntry(uniID, [int(start), int(end)])
			
			for line in lines[1:]:
				columns 	= line.split('\t')
				function	= []
				if len(columns) > 2:
					try:
							function.append(columns[2])
							function.append(int(columns[3]))
							function.append(int(columns[4]))
							function.append((columns[8]))
					except:
						#should add a logfile to store those conditions somewere (TODO)
						continue
				if len(function) != 0:						
					newEntry.functions.append(function)
			allEntrys[uniID] = newEntry
	
	return allEntrys

def loadAlignment(path):
	allAlignemnts = {}
	infile  = open(path, 'r')
	blocks  = infile.read().split('>')
	for block in blocks:
		if len(block) > 0:
			allAlignemnts[block.split('\n')[1]] = block.split('\n')[1:]
	return allAlignemnts

def generateRBHs(fasta1, fasta2, tmpFolder):
	'''
	Start the script to obtain RBHs
	'''
	# create tmp outputFolder 
	if not os.path.exists(tmpFolder):
		os.makedirs(tmpFolder)
	dir_path = os.path.dirname(os.path.realpath(__file__))
	print(dir_path)
	p = subprocess.Popen('python '+dir_path+'/classes/RBH_BHH_identification.py --prefix '+tmpFolder+' --input1 '+fasta1+' --input2 '+fasta2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	for stdout_line in iter(p.stdout.readline, ""):
		print(stdout_line)
	
	# wait for process to finish
	retval = p.wait()

def generateAlignments(path_to_fasta, tmpFolder, nr_entries, processes):
	"""
	Generate alignments using blast, outformat 0
	TODO: parallelisation in efficient way
	"""
	#generateTmpFolder
	if not os.path.exists(tmpFolder):
		os.makedirs(tmpFolder)
	
	process = processes
	for i in range(1, nr_entries+1):
		#print('blastp  -out '+tmpFolder+str(i)+' -query ' + path_to_fasta + '/' +str(i)+'a -subject ' +path_to_fasta+ '/'+str(i)+'b -outfmt 0')
		p = subprocess.Popen('blastp  -out '+tmpFolder+str(i)+' -query ' + path_to_fasta + '/' +str(i)+'a -subject ' +path_to_fasta+ '/'+str(i)+'b -outfmt 0', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		retval = p.wait()

def getRBHs(path_to_file):
	"""
	Function to map ids. This might be necessary to link
	the information from pdb database to the used gene idenfier
	"""
	with open(path_to_file, "r") as DataFile:
		mapping = {}
		entrys = DataFile.read().split('\n')
		for entry in entrys[1:]:
			tmp 	= entry.split('\t')
			if len(tmp) > 1:
				tmp2	= tmp[1].split('|')
				if tmp[0].split('.')[0] in mapping:
					if mapping[tmp[0].split('.')[0]][2] < tmp[3]:
						mapping[tmp[0].split('.')[0]]= [tmp2[1], tmp[2], tmp[3], tmp[0]]
				else:
					mapping[tmp[0].split('.')[0]]= [tmp2[1], tmp[2], tmp[3], tmp[0]]
		return mapping

def transitiveMapping(mapPOI_TRD, mapPOI_IMD):
	'''
	generates two outputs, produces the relationship between TRD and IMD
	'''
	mapIMD_TRD, mapTRD_IMD = {}, {}
	for i in mapPOI_TRD:
		if i in mapPOI_IMD:
			mapTRD_IMD[mapPOI_TRD[i][0]] = mapPOI_IMD[i]
			mapIMD_TRD[mapPOI_IMD[i][0]] = mapPOI_TRD[i]
	
	for i in mapPOI_IMD:
		if i in mapPOI_TRD:
			mapTRD_IMD[mapPOI_TRD[i][0]] = mapPOI_IMD[i]
			mapIMD_TRD[mapPOI_IMD[i][0]] = mapPOI_TRD[i]
					
	return mapIMD_TRD, mapTRD_IMD

def reverseMapping(mapping):
	mapping2 = {}
	for i in mapping:
		mapping2[mapping[i][0]] = [i, mapping[i][1], mapping[i][2], mapping[i][3]] 
	return mapping2	

def splitFastaFilesPairwise(fasta1, fasta2, mapping1, tmpFolder):
	'''
	Create a folder with files for pairwise blast alignment
	'''

	# create tmp outputFolder 
	if not os.path.exists(tmpFolder):
		os.makedirs(tmpFolder)	
	
	# create actuall split files	
	counter	  = 1
	for entry in fasta1:
		#print(entry)
		if entry.split('.')[0] in mapping1:
			outfile1 = open(tmpFolder+str(counter)+'a', 'w')
			outfile1.write('>'+entry.split('.')[0]+'\n')
			outfile1.write(fasta1[entry])
			outfile1.close()
			outfile2 = open(tmpFolder+str(counter)+'b', 'w')
			outfile2.write('>'+mapping1[entry.split('.')[0]][0]+'\n')
			outfile2.write(fasta2[mapping1[entry.split('.')[0]][0]])			
			outfile2.close()
			counter += 1

	return counter
	
def mergeSplit(inpath, outpath, nr_entries):
	"""
	merge into nice and parseable format (necessary controll output)
	TODO: check if to large proteins eventually collide with the format
	"""
	# define outpath
	outfile			= open(outpath, 'w')
	allAlignments 	= {}
	
	# iterate throug entires and load/parse them
	for i in range(1,nr_entries):
		infile		= open(inpath+str(i), 'r')
		lines 		= infile.read().split('\n')
		infile.close()
		queryID 	= lines[3].split('= ')[1]
		subjectID 	= lines[7].split('= ')[1]
		info		= lines[13]
		query		= ''
		meta		= ''
		subject		= ''

		nextMeta 	= False
		lengthMeta	= 0
		for line in lines[13:]:
			if len(line.split(' ')) == 0:
				nextMeta = False
				continue
			if nextMeta == True:
				meta		+= line[12:lengthMeta+12]
				nextMeta = False
			if line[0:5] == 'Query':
				query 		 += line[12:].split('  ')[0]
				nextMeta 	 = True
				lengthMeta = len(line[12:].split('  ')[0])
			if line[0:5] == 'Sbjct':
				subject		 += line[12:].split('  ')[0]
				nextMeta 	 = False
		outfile.write('>\n')
		outfile.write(queryID+'\n')
		outfile.write(subjectID+'\n')
		outfile.write(info+'\n')
		outfile.write(query+'\n')
		outfile.write(meta+'\n')
		outfile.write(subject+'\n')	
		allAlignments[queryID] = [queryID, subjectID, info, query, meta, subject]

	return 	allAlignments		

def checkAnnotationQualitiy(seq1, seq2, meta, start, end, function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv):
	'''
	Function to obtain alignment quality
	Important TODO: Deal with Disulphide-Bonds
	'''
	start_new = start
	end_new   = end
		
	# adapt start position of annotation
	gaps = seq2[:start].count('-')
	start_new += gaps
	for i in range(start, start+gaps):
		try:
			if seq2[i] == '-':
				start_new += 1
		except:
			continue
	
	# adapt end position of annotation
	gaps = seq2[:end].count('-')
	end_new += gaps
	for i in range(end, end+gaps):
		try:
			if seq2[i] == '-':
				end_new += 1
		except:
			continue

	region = meta[start_new-1:end_new]
		
	if len(region) > 0:
		positives 		= (len(region)-region.count(' '))/len(region)
		idents	  		= (len(region)-region.count(' ')-region.count('+'))/len(region)
		positives_all 	= (len(meta)-meta.count(' '))/len(meta)
		idents_all		= (len(meta)-meta.count(' ')-meta.count('+'))/len(meta)
		quality = str(round(idents,2))+'_'+str(round(idents_all,2))+'_'+str(round(positives, 2))+'_'+str(round(positives_all, 2))
		#print(str(round(idents,2)),str(round(idents_all,2)),str(round(positives, 2)),str(round(positives_all, 2)))

		# get position of reference:
		gaps = seq1[:start_new].count('-')
		start_ref = start_new-gaps
		gaps = seq1[:end_new].count('-')
		end_ref = end_new-gaps
		
		if end_ref-start_ref == 0:
			add = sizeOfEnvironment
			region2						 = meta[(start_new-add) : (end_new+add+1)]
			#print(region2)
			try:	
				positives_environment 		 = (len(region2)-region2.count(' '))/len(region2)
				idents_environment 			 = (len(region2)-region2.count(' ')-region2.count('+'))/len(region2)
			except:
				positives_environment = 0
				idents_environment	  = 0
				
			add2 = sizeOfEnvironment
			region3						 = meta[(start_new-add2) : (end_new+add2+1)]
			#if idents == 1 and positives >= amountOfPositivesAll:
			#	print(idents_environment, amountOfPositivesEnv)
			try:	
				positives_environment2 		 = (len(region3)-region3.count(' '))/len(region3)
				idents_environment2 		 = (len(region3)-region3.count(' ')-region3.count('+'))/len(region3)
			except:
				positives_environment2 = 0
				idents_environment2	  = 0	
						
		if positives >= amountOfPositivesAll:
			if end_ref-start_ref == 0:
				if idents == 1 and idents_environment >= amountOfPositivesEnv:
					#if function[0] == 'Active site':
					#	print(region2, region2.count(' '))
					return(start_ref, end_ref, quality)
			else:
				#if function[0] == 'Domain':
				#	print(positives)
				return(start_ref, end_ref, quality)

def checkAnnotationQualitiy3(seq1, seq2, meta, start, end, function):
	'''
	Function to obtain alignment quality
	Important TODO: Deal with Disulphide-Bonds
	'''
	start_new = start
	end_new   = end
		
	# adapt start position of annotation
	gaps = seq2[:start].count('-')
	start_new += gaps
	for i in range(start, start+gaps):
		try:
			if seq2[i] == '-':
				start_new += 1
		except:
			continue
	
	# adapt end position of annotation
	gaps = seq2[:end].count('-')
	end_new += gaps
	for i in range(end, end+gaps):
		try:
			if seq2[i] == '-':
				end_new += 1
		except:
			continue

	region = meta[start_new-1:end_new]
	if len(region) > 0:
		positives 		= (len(region)-region.count(' '))/len(region)
		idents	  		= (len(region)-region.count(' ')-region.count('+'))/len(region)
		positives_all 	= (len(meta)-meta.count(' '))/len(meta)
		idents_all		= (len(meta)-meta.count(' ')-meta.count('+'))/len(meta)
		quality = str(round(idents,2))+'_'+str(round(idents_all,2))+'_'+str(round(positives, 2))+'_'+str(round(positives_all, 2))
		#print(str(round(idents,2)),str(round(idents_all,2)),str(round(positives, 2)),str(round(positives_all, 2)))
		
		# get position of reference:
		gaps = seq1[:start_new].count('-')
		start_ref = start_new-gaps
		gaps = seq1[:end_new].count('-')
		end_ref = end_new-gaps
				
		if positives >= 0.8:
			if end_ref-start_ref == 0:
				if idents == 1:
					return(start_ref, end_ref, quality)
			else:
				return(start_ref, end_ref, quality)

def createDatabaseImprovement(outPath, trdGFF_entries, imdGFF_entries, alignments_IMD_TRD, alignments_POI_TRD, mapIMD_TRD, mapTRD_POI, mapTRD_IMD, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv):
	'''
	Assumptions: If uniport started to annotate, we assume that the correspoinding tag was sucessuflly anotated. 
	'''
	counter 		= 0
	outfile 		= open(outPath+'IMD.gff', 'w') 
	sharedEntries 	= {}
	
	# statistic and count variables:
	gene_nr 	= {'TRD':0, 'IMD_0':0, 'IMD_1':0, 'IMD_2':0}
	tag_count	= {'TRD_0':{}, 'TRD_1':{}, 'IMD_0':{}, 'IMD_1':{}, 'IMD_2':{}}
	
	# starting annotation process
	outfile.write('##gff-version 3'+'\n')	
	for gff in imdGFF_entries:
		gene_nr['IMD_0'] += 1
		outfile.write('##sequence-region '+gff+' '+str(imdGFF_entries[gff].amino_acids[0])+' '+str(imdGFF_entries[gff].amino_acids[1])+'\n')
		if gff in alignments_IMD_TRD:
			gene_nr['IMD_2'] += 1
			counter += 1
			allTags = {}
			improvement = False
			sharedEntries[mapIMD_TRD[gff][0]] = True
			# get and output all entries, which are already annotated
			for function in imdGFF_entries[gff].functions:
				if function[0] not in allTags:
					allTags[function[0]] = True
				outfile.write(gff+'\tUniProtKB'+'\t'+str(function[0])+'\t'+str(function[1])+'\t'+str(function[2])+'\t.\t.\t.\t'+function[3]+'\n')
				# count tags
				if function[0] in tag_count['IMD_2']:
					tag_count['IMD_2'][function[0]] += 1
				else:
					tag_count['IMD_2'][function[0]]  = 1
				if function[0] in tag_count['IMD_0']:
					tag_count['IMD_0'][function[0]] += 1
				else:
					tag_count['IMD_0'][function[0]]  = 1					
			# check whether something from TRD can be Transfered to IMD
			for function in trdGFF_entries[mapIMD_TRD[gff][0]].functions:
				if function[0] not in allTags:
					if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant']:
						if function[0] not in ['Disulfide bond']:
							values = checkAnnotationQualitiy(alignments_IMD_TRD[gff][3], alignments_IMD_TRD[gff][5], alignments_IMD_TRD[gff][4], function[1], function[2], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
						else:
							values1 = checkAnnotationQualitiy(alignments_IMD_TRD[gff][3], alignments_IMD_TRD[gff][5], alignments_IMD_TRD[gff][4], function[1], function[1], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
							values2 = checkAnnotationQualitiy(alignments_IMD_TRD[gff][3], alignments_IMD_TRD[gff][5], alignments_IMD_TRD[gff][4], function[2], function[2], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
							if values1 != None and values2 != None:
								values = (values1[0], values2[0], values1[2]+'_'+values2[2])
							else:
								values = None							
						if values != None:
							outfile.write(gff+'\tUniProtKB'+'\t'+str(function[0])+'\t'+str(values[0])+'\t'+str(values[1])+'\t.\t.\t.\t'+function[3]+'|annoTrans:'+alignments_IMD_TRD[gff][1]+'_'+str(values[2])+'\n')
							# count tags
							if function[0] in tag_count['IMD_1']:
								tag_count['IMD_1'][function[0]] += 1
							else:
								tag_count['IMD_1'][function[0]]  = 1
							if function[0] in tag_count['IMD_2']:
								tag_count['IMD_2'][function[0]] += 1
							else:
								tag_count['IMD_2'][function[0]]  = 1								
							improvement = True
			if 	improvement == True:
				gene_nr['IMD_1'] += 1
		else:
			for function in imdGFF_entries[gff].functions:
				i=0
				outfile.write(gff+'\tUniProtKB'+'\t'+str(function[0])+'\t'+str(function[1])+'\t'+str(function[2])+'\t.\t.\t.\t'+function[3]+'\n')
				# add tag counts
				if function[0] in tag_count['IMD_2']:
					tag_count['IMD_2'][function[0]] += 1
				else:
					tag_count['IMD_2'][function[0]]  = 1
				if function[0] in tag_count['IMD_0']:
					tag_count['IMD_0'][function[0]] += 1
				else:
					tag_count['IMD_0'][function[0]]  = 1					

	# identify new entries which were mapped bu do not exists in the dababase which shall be improved
	alignments_POI_TRD_nonImprovend = {}
	for i in alignments_POI_TRD:
		gene_nr['TRD'] += 1
		if alignments_POI_TRD[i][1] not in sharedEntries:
			alignments_POI_TRD_nonImprovend[i] = alignments_POI_TRD[i]
	
	# create new database for missing entries
	NWD_nr, NRW_tags = createNewDatabase(outPath, trdGFF_entries, alignments_POI_TRD_nonImprovend, mapTRD_POI, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
	
	# get TRD_TAGS
	for gff in trdGFF_entries:
		for function in trdGFF_entries[gff].functions:
			if function[0] in tag_count['TRD_0']:
				tag_count['TRD_0'][function[0]] += 1
			else:
				tag_count['TRD_0'][function[0]]  = 1
			if gff in mapTRD_IMD:
				if function[0] in tag_count['TRD_1']:
					tag_count['TRD_1'][function[0]] += 1
				else:
					tag_count['TRD_1'][function[0]]  = 1				
						
	# process count data for return
	tag_count['NWD_1']	= NRW_tags['NWD_1']
	tag_count['NWD_0']	= NRW_tags['NWD_0']
	#tag_count['TRD']	= NRW_tags['TRD']
	gene_nr['NWD'] 		= NWD_nr['NWD']

	#return statistics
	return gene_nr, tag_count
	
def createNewDatabase(outPath, trdGFF_entries, alignments_POI_TRD, mapTRD_POI, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv):
	"""
	Create a new database by evaluating all entries in the database which shall be transferred
	"""
	gene_nr   = {'NWD':0}
	tag_count = {'TRD':{}, 'NWD_0':{}, 'NWD_1':{}}
	
	outfile 		= open(outPath+'NWD.gff', 'w') 
	outfile.write('##gff-version 3'+'\n')	
	for gff in trdGFF_entries:
		existingEntry = False
		try:
			if mapTRD_POI[gff][0] in alignments_POI_TRD:
				existingEntry = True
		except:
			continue
		if existingEntry == True:
			strings = []
			strings.append('##sequence-region '+mapTRD_POI[gff][0]+' '+str(trdGFF_entries[gff].amino_acids[0])+' '+str(trdGFF_entries[gff].amino_acids[1])+'\n')
			improvement = False
			
			for function in trdGFF_entries[gff].functions:
				# count tags
				if function[0] in tag_count['TRD']:
					tag_count['TRD'][function[0]] += 1
				else:
					tag_count['TRD'][function[0]]  = 1					
				if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant']:
					if function[0] not in ['Disulfide bond']:
						values = checkAnnotationQualitiy(alignments_POI_TRD[mapTRD_POI[gff][0]][3], alignments_POI_TRD[mapTRD_POI[gff][0]][5], alignments_POI_TRD[mapTRD_POI[gff][0]][4], function[1], function[2], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
					else:
						values1 = checkAnnotationQualitiy(alignments_POI_TRD[mapTRD_POI[gff][0]][3], alignments_POI_TRD[mapTRD_POI[gff][0]][5], alignments_POI_TRD[mapTRD_POI[gff][0]][4], function[1], function[1], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
						values2 = checkAnnotationQualitiy(alignments_POI_TRD[mapTRD_POI[gff][0]][3], alignments_POI_TRD[mapTRD_POI[gff][0]][5], alignments_POI_TRD[mapTRD_POI[gff][0]][4], function[2], function[2], function, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
						if values1 != None and values2 != None:
							values = (values1[0], values2[0], values1[2]+'_'+values2[2])
						else:
							values = None	
					if function[0] in tag_count['NWD_0']:
						tag_count['NWD_0'][function[0]] += 1
					else:
						tag_count['NWD_0'][function[0]]  = 1
					if values != None:
						improvement = True
						strings.append(gff+'\tUniProtKB'+'\t'+str(function[0])+'\t'+str(values[0])+'\t'+str(values[1])+'\t.\t.\t.\t'+function[3]+'|annoTrans:'+alignments_POI_TRD[mapTRD_POI[gff][0]][1]+'_'+str(values[2])+'\n')
						# count tags
						if function[0] in tag_count['NWD_1']:
							tag_count['NWD_1'][function[0]] += 1
						else:
							tag_count['NWD_1'][function[0]]  = 1	
			if improvement == True:
				gene_nr['NWD'] += 1
				for i in strings:
					outfile.write(i)
	return gene_nr, tag_count

def statistics(outPath, counts):
	outfile 		= open(outPath+'statistics.csv', 'w') 
	outfile.write('Abbreviations:\n')
	outfile.write('TRD_0: All entires from database you want to Transfer\n')
	outfile.write('TRD_1: All entires from database you want to Transfer which match with a IMD entry\n')
	outfile.write('IMD_0: All database entries from the database you want to improve\n')
	outfile.write('IMD_1: Intersection of IMD_0 and TRD_1\n')
	outfile.write('IMD_2: Entries from TRD_1 which passed quality controll\n')
	outfile.write('NWD_0: All tags which can potentially be annotated in a new database\n')
	outfile.write('NWD_1: All tags from NWD_= which passed quality criteria\n')
	outfile.write('NWD_Non: Entries which did not pass quality filter for NWD\n\n')

	outfile.write('subset\tgene_count\n')
	#for i in counts[0]:
	#	outfile.write(i+'\t'+str(counts[0][i])+'\n' )
	#outfile.write('NWD_Non\t'+ str(counts[0]['TRD']-counts[0]['IMD_1']-counts[0]['NWD'])+'\n\n')

	outfile.write('Tags per subset:\n')
	allTags = {}
	# get all possible Tags
	
	for i in counts[1]:
		for tag in counts[1][i]: 
			if tag not in allTags and tag not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant', 'Sequence uncertainty', 'Non-standard residue', 'Non-terminal residue']:
				allTags[tag] = 0
	outfile.write('Tag\t')
	for i in counts[1]:
		outfile.write(i+'\t')
	outfile.write('\n')
	
	for tag in allTags:
		string = tag+'\t'
		for i in counts[1]:
			if tag in counts[1][i]:
				string+= str(counts[1][i][tag])+'\t'
			else:
				string+='0\t'
		outfile.write(string+'\n')

def getMaxOverlap(function, targets):
	maxOverlapSites_function = 0
	length_of_target   = 0
	for target in targets:
		maxOverlap = 0
		if function[1]  <= target[1] and function[2] >= target[1] and function[2] <= target[2]:
			maxOverlap = function[2]-target[1]
			#print('1', target, maxOverlap)
		if function[1]  >= target[1] and function[1] <= target[2] and function[2] >= target[2]:
			maxOverlap = target[2]-function[1]
			#print('2', target, maxOverlap)			
		if function[1]  <= target[1] and function[2] >= target[2]:
			maxOverlap = target[2] - target[1] 
			#print('3', target, maxOverlap)			
		if function[1]  >= target[1] and function[2] <= target[2]:
			maxOverlap = function[2] - function[1] 
			#print('4', target, maxOverlap)

		if maxOverlap > maxOverlapSites_function:
			maxOverlapSites_function = maxOverlap
			maxOverlapSites_target = maxOverlap
			length_of_target = target[2]-target[1]
	
	#print(function, targets, maxOverlapSites_function, float((function[2]-function[1])), maxOverlapSites_function/float((function[2]-function[1])), maxOverlapSites_function/float(length_of_target))
	
	try:
		return 	maxOverlapSites_function/float((function[2]-function[1])), maxOverlapSites_function/float(length_of_target)#, maxOverlapSites_target
	except:
		return 	maxOverlapSites_function/float(1+(function[2]-function[1])), maxOverlapSites_function/float(1+length_of_target)#, maxOverlapSites_target
		
def compareToGoldStandard(outpath, gff_ref, gff_new, mapPOI_TRD, mapTRD_POI):
	'''
	Function to compare a new dabase with a gold standard. Be careful with interpretation, gold standards are rare
	'''
	
	outfile 	= open(outPath+'comparence_to_gold_standard.csv', 'w')
	TP_sites	= 0		
	FP_sites	= 0
	FN_sites	= 0
	
	Overlap_reg	= 0
	All_regions = 0
	
	c1 = 0
	c2 = 0

	for tag_sample in [['Domain'], ['Motif'],  ['Transit peptide','Peptide','Signal peptide','Propeptide'], ['Region'], ['Transmembrane', 'Intramembrane'],['DNA binding', 'Zinc finger'],['Active site', 'Binding site','Site', 'Metal binding', 'Calcium binding', 'Modified residue','Lipidation','Glycosylation'],['Turn', 'Helix', 'Beta strand'] , ]:
	
		#print(tag_sample)
		#for functions:		
		if tag_sample in [['Domain'], ['Motif'],['Transit peptide','Peptide','Signal peptide','Propeptide'],['Region'], ['Transmembrane', 'Intramembrane'],['DNA binding', 'Zinc finger'], ['Turn', 'Helix', 'Beta strand']]:

			overlaps_new   = []
			overlaps_ref   = []
			functions_count = 0
			for i in gff_new:
				functions_new 	= []
				functions_ref 	= []				
				# get all positions in comparable format
				tmp = i.split('|')[1]
				if tmp in gff_ref:
					for function in gff_ref[tmp].functions:
						#if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant', 'Sequence uncertainty', 'Non-standard residue', 'Non-terminal residue', 'Repeat', 'Cross-link', 'Initiator methionine', 'Mutagenesis']:
						if function[0] in tag_sample:
							functions_ref.append([function[0], function[1], function[2]])
					for function in gff_new[i].functions:
						#if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant', 'Sequence uncertainty', 'Non-standard residue', 'Non-terminal residue', 'Repeat', 'Cross-link', 'Initiator methionine', 'Mutagenesis']:
						if function[0] in tag_sample:
							functions_new.append([function[0], function[1], function[2]])		
					
					for function in functions_new:
						functions_count += 1
						max_overlap = getMaxOverlap(function, functions_ref)
						overlaps_new.append(max_overlap[0])
						overlaps_ref.append(max_overlap[1])

			#print(sum(overlaps_new)/len(overlaps_new), functions_count)
			#print(sum(overlaps_ref)/len(overlaps_ref))
		else:
			functions_new = []
			functions_ref = []
			
			#print(len(gff_new))
			
			# get all positions in comparable format
			for i in gff_new:	
				functions_new = []
				functions_ref = []
				
				# get all positions in comparable format
				tmp = i.split('|')[1]
				if tmp in gff_ref:
					#print('e')
					for function in gff_ref[tmp].functions:
						#if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant', 'Sequence uncertainty', 'Non-standard residue', 'Non-terminal residue', 'Repeat', 'Cross-link', 'Initiator methionine', 'Mutagenesis']:
						if function[0] in tag_sample:
							functions_ref.append([function[0], function[1], function[2]])
					for function in gff_new[i].functions:
						#if function[0] not in ['Alternative sequence', 'Chain', 'Sequence conflict', 'Topological domain','Compositional bias', 'Natural variant', 'Sequence uncertainty', 'Non-standard residue', 'Non-terminal residue', 'Repeat', 'Cross-link', 'Initiator methionine', 'Mutagenesis']:
						if function[0] in tag_sample:
							functions_new.append([function[0], function[1], function[2]])		
				tp_local = 0
				for function in functions_new:
					if function in functions_ref:
						TP_sites += 1
						tp_local += 1
					else:
						FP_sites += 1
					FN_sites += (len(functions_ref)-tp_local)		
			#print(TP_sites, FP_sites, FN_sites, TP_sites/(TP_sites+FP_sites), TP_sites/(TP_sites+FN_sites))
			try:
				print(TP_sites, TP_sites+FP_sites, TP_sites/(TP_sites+FP_sites))
			except:
				print(TP_sites, TP_sites+FP_sites, 0)
							
if __name__ == '__main__':
	'''
	AnnoTrans is a tool to improve an uniport entry unsing orthologue 
	and paralogue gene annototation. The procedure is the following:
	
	Input:
		uniprot:
			- uniprot protein database of annotations you want to transfer (gff format) (TRD)
			- fasta files of the uni proteins (fasta format)
			 optional:
			- uniprot protein database of annotations you want to improve (gff format)  (IMD)
			- fasta files of the uni proteins you want to improve
		structural annoation:
			- fasta file of the transkripts from the structural annotation
	Output:
		- uniprot protein database of annotations for the protein of interest
		 optional:
		- improved uniprot protein database with additional annotations for the proteins of interest (POIs)
			
	1.	 a. identify RBHs 				(orthologue genes)
		 b. identify best blast hits	(paralogue genes)
	2.	 a. pairwise alignment of POI and TRD
		 (optional)
		 b. pairwise alignment of IMD and TRD
		 c. pairwise alignment of POI and IMD
	3. 	 a. evaluate annoation and transfer if quality criteria are fullfilled
	4.	 a. Output new gff annotation file, outpuz mapping file POIs -> TRD and IMD, output improvement statisics
	'''
	# Settings
	outPath 				= None
	trdGFF  				= None
	imdGFF 			 		= None
	trdFAS					= None
	imdFAS					= None
	sAnFAS					= None
	sizeOfEnvironment 		= 5
	amountOfPositivesAll	= 0.8
	amountOfPositivesEnv	= 0.8
	
	# get settings and path
	for i in range(1,len(argv)):
		if argv[i] == '-outPath':
			outPath = argv[i+1]
		if argv[i] == '-trdGFF':
			trdGFF = argv[i+1]		
		if argv[i] == '-imdGFF':
			imdGFF = argv[i+1]				
		if argv[i] == '-trdFAS':
			trdFAS = argv[i+1]
		if argv[i] == '-imdFAS':
			imdFAS = argv[i+1]
		if argv[i] == '-sAnFAS':
			sAnFAS = argv[i+1]
		if argv[i] == '-env':
			sizeOfEnvironment = int(argv[i+1])
		if argv[i] == '-posAll':
			amountOfPositivesAll = float(argv[i+1])
		if argv[i] == '-posEnv':
			amountOfPositivesEnv = float(argv[i+1])	
		
	### 1. identify RBHs	###
	# POI vs TRD
	#generateRBHs(sAnFAS, trdFAS, outPath+'/POIvsTRD')
	#print('POI vs TRD done!')
	# POI vs IMD
	#generateRBHs(sAnFAS, trdFAS, outPath+'/POIvsIMD')
	#print('POI vs TRD done!')
	
	# create (transitive) mapping
	mapPOI_TRD 				= getRBHs(outPath+'/POIvsTRD/RBH_file.txt')
	#mapPOI_IMD 				= getRBHs(outPath+'/POIvsIMD/RBH_file.txt')
	#mapIMD_TRD, mapTRD_IMD	= transitiveMapping(mapPOI_TRD, mapPOI_IMD)
	mapTRD_POI				= reverseMapping(mapPOI_TRD)

	### 2. do pairwise alignments to obtain blast meta information in a parsable way ###
	
	## split up fasta files
	
	# IMD vs TRD
	#nr_entries_IMD_TRD = splitFastaFilesPairwise(load_multiple_fasta_file(imdFAS, 'uniprot'), load_multiple_fasta_file(trdFAS, 'uniprot'), mapIMD_TRD, outPath+'/IMD_TRD_SPLIT/')
	# POI vs TRD
	#nr_entries_POI_TRD = splitFastaFilesPairwise(load_multiple_fasta_file(sAnFAS, 'annotation'), load_multiple_fasta_file(trdFAS, 'uniprot'), mapPOI_TRD, outPath+'/POI_TRD_SPLIT/')

	#print(nr_entries_IMD_TRD, nr_entries_POI_TRD)
	#print(nr_entries_POI_TRD)
	#nr_entries_IMD_TRD = 20267	#
	nr_entries_POI_TRD = 8727#3789#40340
	# start BLAST for IMD vs TRD
	#generateAlignments(outPath+'/IMD_TRD_SPLIT/', outPath+'/IMD_TRD_SPLIT/Alignments/', nr_entries_IMD_TRD, 50)
	#generateAlignments(outPath+'/POI_TRD_SPLIT/', outPath+'/POI_TRD_SPLIT/Alignments/', nr_entries_POI_TRD, 50)
	
	# merge blast alignment results
	#alignments_IMD_TRD = mergeSplit(outPath+'/IMD_TRD_SPLIT/Alignments/', outPath+'alignments_IMD_TRD.aln' , nr_entries_IMD_TRD)
	#alignments_POI_TRD = mergeSplit(outPath+'/POI_TRD_SPLIT/Alignments/', outPath+'alignments_POI_TRD.aln' , nr_entries_POI_TRD)
	
	#alignments_IMD_TRD = loadAlignment(outPath+'alignments_IMD_TRD.aln')
	alignments_POI_TRD = loadAlignment(outPath+'alignments_POI_TRD.aln')
	
	### 3. Transfer annotation and output new databases in gff format
	trdGFF_entries	=	loadUniprotGFF(trdGFF)
	#imdGFF_entries	=	loadUniprotGFF(imdGFF)

	#counts = createDatabaseImprovement(outPath, trdGFF_entries, imdGFF_entries, alignments_IMD_TRD, alignments_POI_TRD, mapIMD_TRD, mapTRD_POI, mapTRD_IMD, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
	counts = createNewDatabase(outPath, trdGFF_entries, alignments_POI_TRD, mapTRD_POI, sizeOfEnvironment, amountOfPositivesAll, amountOfPositivesEnv)
	### 4. create file with counts and statistics about the annotation process
	#statistics(outPath, counts)
	#print(amountOfPositivesEnv)
	### 5. evaluation with gold standard
	#compareToGoldStandard(outPath, loadUniprotGFF('/vol/cluster-data/jloers/annoTrans/evalutation_1/data/rice.gff'), loadUniprotGFF('/vol/cluster-data/jloers/annoTrans/evalutation_1/NWD.gff'), mapPOI_TRD, mapTRD_POI)
	compareToGoldStandard(outPath, loadUniprotGFF('/vol/cluster-data/jloers/annoTrans/evalutation_2/data/arabidopsis-uniprot.gff'), loadUniprotGFF('/vol/cluster-data/jloers/annoTrans/evalutation_2/NWD.gff'), mapPOI_TRD, mapTRD_POI)
