#!/usr/bin/python

import xml.etree.ElementTree as ET
import re
import numpy as np
import difflib 
import time
startTime = time.time()
import os
import sys
import chimera
from chimera import runCommand as rc
import pickle 

class SSM_search_results():
	"""Class representing the outcome of SSM search - the xml file - the results 
	of a structural comparison between a target against all the hits SSM yields."""

	def __init__(self):

		self.xml_file = ''

		# info on target
		self.root = []
		self.target_pdb = []
		self.target_chain = []
		self.target_nRes = []
		self.target_matching_res = ''  	# the residues matched on target stored as structural array
										# 1st col==hit pdb, 2nd col==target's starting res of match, 3rd col==target's ending res of match
		# for storing info on hits
		self.match = []
		self.hit_nRes_percSSEs = {}
		self.hit_names = {}
		self.matched_target_percSSEs = {}
		self.RTMatrices = {}
		self.Q_scores = {}
		self.hits_matching_res = {}  # the residues matched on hits
		
		# results of after their correspodning filter being applied
		self.first_filter_pdbs = {}
		self.second_filter_pdbs = {}
		self.third_filter_pdbs = {}

		self.model0 = []
		self.all_hits_clusters = {}

	def parse_xml(self):
		"""To parse the xml file - the SSM output file.
		The SSM is confusingly calling our target as query, and our hits as targets. Be aware!"""

		tree = ET.parse(self.xml_file)
		self.root = tree.getroot()

		# reading the input parameters
		ssmInput = self.root.find("SSMInput")
		self.target_perc_param = ssmInput.find("percent1")
		self.hit_perc_param = ssmInput.find("percent2")
		target = ssmInput.find('query') # the target (SSM calls it query, SSM calls our hits as targets, a bit confusing)
		self.target_pdb = target.find('pdbcode')
		selection = target.find('selection')
		self.target_chain = selection.find('chains')

		# to find number of residues in our target
		target_for_nRes = self.root.find('Query')
		chain_for_nRes = target_for_nRes.find('chain')
		self.target_nRes = chain_for_nRes.find('Nres')

		# finding out the total number of SSE alignments (for use in self.target_matching_res)
		no_of_align = 0
		for SSEpair in self.root.iter("SSEpair"):
			no_of_align += 1

		# getting the info on matches (between our target and hits found by SSM)
		self.no_of_matches = self.root.find("NofMatches")
		self.match = self.root.findall("Match")
		self.target_matching_res = np.zeros((no_of_align,),dtype=('a4,i4,i4')) # to store all SSE alignments separately
		SSEpairs_fill_pos = 0 # for filling the structural array self.target_matching_res 
		for m in self.match:
			hit = m.find("Target")
			hit_pdb = hit.find("name")

			# for the 1st step in removing redundancy - nRes+%SSE
			hit_nRes = hit.find("Nres")
			hit_percSSE = hit.find("matched-SSE-percent")
			self.hit_nRes_percSSEs[hit_pdb.text] = (int(hit_nRes.text),int(hit_percSSE.text))


			# for the 2nd step in removing redundancy - names of the targets
			hit_name = hit.find("title")
			hit_name = hit_name.text
			hit_name = hit_name.replace('\n','') # getting rid of ALL newline characters
			hit_name = re.sub('\s{2,}', ' ', hit_name) # using regexes remove two or more spaces
			if hit_name[-1] == ' ': # if there is still a single empty space in the end, remove it
				hit_name = hit_name[:-1]
			# how to break this stupid line below into many!?
			poss_words = 'IN COMPLEX WITH EXPLORING BACKBONE PATTERN RESOLUTION HUMAN HOMO SAPIENS THE AN REDUCED BOUND TO A CRYSTAL STRUCTURE OF MUTANT NMR THREE-DIMENSIONAL STRUCTURES PROPERTIES SOLUTION CATALYTIC CYCLE FROM CATALYTIC CORE COMPONENT SUBSTRATE AND CATALYSIS DOMAIN EDITING INACTIVATED BY COMPLEXED BINDING STATE PUTATIVE E.COLI C-TERMINAL OF DROSOPHILA MELANOGASTER MODIFIED BY CLOSTRIDIUM ACETOBUTYLICUM NORTHEAST STRUCTURAL GENOMICS CONSORTIUM RESPONSE REGULATOR TARGET DESIGN MOUSE PROTEIN NUCLEOTIDE MUTANT INSIGHTS INTO INHIBITOR DESIGN'
			word_list = poss_words.split()
			for word in word_list:
				hit_name = hit_name.replace(word+' ','') # removing words from the names if present
			self.hit_names[hit_pdb.text] = hit_name # cool, now here having lovely strings

			# for the 3rd step in removing redundancy - cap on percentage of target SSE being matched
			matching_target = m.find("Query")
			matched_target_percSSE = matching_target.find("matched-SSE-percent")
			self.matched_target_percSSEs[hit_pdb.text] = int(matched_target_percSSE.text)

			# getting the RTMatrices for superposition of our target wits hits
			RTMatrix = m.find('RTMatrix')
			self.RTMatrices[hit_pdb.text] = ([[RTMatrix.find('Rxx').text, RTMatrix.find('Rxy').text, RTMatrix.find('Rxz').text, RTMatrix.find('Tx').text]])
			self.RTMatrices[hit_pdb.text].append([RTMatrix.find('Ryx').text, RTMatrix.find('Ryy').text, RTMatrix.find('Ryz').text, RTMatrix.find('Ty').text])
			self.RTMatrices[hit_pdb.text].append([RTMatrix.find('Rzx').text, RTMatrix.find('Rzy').text, RTMatrix.find('Rzz').text, RTMatrix.find('Tz').text])

			# the Q scores of matches (the geometrical measure of similarity between them)
			Q_score = m.find('Q-score')
			self.Q_scores[hit_pdb.text] = float(Q_score.text)  # the most important, the geometrical measure of structural similarity (encompasses rmsd + number of aligned residues)
	
			# gettting the individual SSE pairs of matches (by looking at the start and end residue of that particular SSE on the target being matched)
			SSEMatch = m.find("SSEMatch")
			SSEpairSeparately = SSEMatch.findall("SSEpair")	    
			for pair in SSEpairSeparately:
				t = pair.find("Query")
				iRes = t.find("initRes")
				iNum = iRes.find("seqNum")
				eRes = t.find("endRes")
				eNum = eRes.find("seqNum")
				self.target_matching_res[SSEpairs_fill_pos] = ( hit_pdb.text, int(iNum.text), int(eNum.text) )
				SSEpairs_fill_pos += 1

				# getting the start and end of matching SSEs for hits as well
				h = pair.find("Target")
				hiRes = h.find("initRes")
				hiNum = hiRes.find("seqNum")
				heRes = h.find("endRes")
				heNum = heRes.find("seqNum")
				
				if hit_pdb.text not in self.hits_matching_res.keys():
					self.hits_matching_res[hit_pdb.text] = range(int(hiNum.text),int(heNum.text)+1)
				else:
					l = self.hits_matching_res[hit_pdb.text]
					for el in range(int(hiNum.text),int(heNum.text)+1):
						l.append(el)

	def remove_redundancy(self, q_score_cutoff, name_thresh, tSSE_cap):
		"""Removing redundancy: 
		1st - keeping only unique hit's nRes and the SSE matching percentage, with Q-score (q_score_cutoff) being larger than threshold
		2nd - remove hits that have ratio of their similarity larger than threshold (name_thresh)
		3rd - remove all remaining hits, whose percentage of target's matching SSE larger than threshold (tSSE_cap)
		Eventually saving hits that remained after filtering as well as number residue of the target to be used in 1st clustering."""

		# 1st filter
		lis = self.hit_nRes_percSSEs.values()
		s = set()
		first_filter_pass = {}
		for el in lis:
			for key,val in self.hit_nRes_percSSEs.iteritems():
				if float(self.Q_scores[key]) > q_score_cutoff:
					if val not in s:
						s.add(val)
						first_filter_pass[key] = val
		self.first_filter_pdbs = first_filter_pass.keys()
		print "\nUp until the end of 1st filtering time elapsed", time.time() - startTime, "seconds."
		print "Number of remaining hits", len(self.first_filter_pdbs)

		# 2nd filter
		crow_to_crow = {} # saying 'like attracts like'
		for el in self.first_filter_pdbs:
			crow_to_crow[el] = [] 
			for el2 in self.first_filter_pdbs:
				if el != el2:
					ratio = difflib.SequenceMatcher(None, self.hit_names[el], self.hit_names[el2]).ratio()
					if ratio > name_thresh:
						crow_to_crow[el].append(el2)
		self.second_filter_pdbs = self.first_filter_pdbs[:]
		for k,v in crow_to_crow.iteritems():
			if k in self.second_filter_pdbs:
				for el in v:
					if el in self.second_filter_pdbs:
						self.second_filter_pdbs.remove(el)
		print "\nUp until the end of 2nd filtering time elapsed", time.time() - startTime, "seconds."
		print "Number of remaining hits", len(self.second_filter_pdbs)

		# 3rd filter
		self.third_filter_pdbs = self.second_filter_pdbs[:]
		for i in self.second_filter_pdbs:
			if self.matched_target_percSSEs[i] >= tSSE_cap: # want target %SSE below threshold 
				self.third_filter_pdbs.remove(i)
		print "\nUp until the end of 3rd filtering time elapsed", time.time() - startTime, "seconds."
		print "Number of remaining hits", len(self.third_filter_pdbs)

		# saving the hits that passed all 3 filters for 1st clustering 
		outfile = open('third_filter_pdbs_'+ self.target_pdb.text + "_" + self.target_chain.text +'.pckl', 'w')
		pickle.dump(self.third_filter_pdbs, outfile)
		outfile.close()

		# saving number of target residues to be used in 1st clustering
		outfile = open('nRes_'+ self.target_pdb.text + "_" + self.target_chain.text +'.pckl', 'w')
		pickle.dump(self.target_nRes.text, outfile)
		outfile.close()

	def write_RTMatrices(self):
		"""This is writing text files that will be used by Chimera
		to superpose the tatget with a hit - saved in the hit's folder"""

		for hit in self.third_filter_pdbs:
			my_path = os.path.abspath(os.curdir)+ "/hit_" + hit
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)
  
			matrix_file = 'RTMatrix_'+self.target_pdb.text+'_'+self.target_chain.text+'_against_'+hit+'.txt'			
			f = open(matrix_file,'wb')

			# opening the model so know the numbering of models
			model1 = chimera.openModels.open(hit,type="PDB")

			if len(model1) == 1:
				f.write('Model 1.0')
			else:
				f.write('Model 1.1') # if there are more than 1 model (as in NMR structures, take the first one)
			rc('~open #1') 	

			for values in self.RTMatrices[hit]:
				f.write('\n\t')
				for val in values:         
					f.write(str(val)+' ')
			f.close()

			# moving back up the directory 
			os.chdir("..")

	def open_target_in_Chimera(self):
		"""Opens the target chain in chimera."""
		
		rc("~open all")
		self.model0 = chimera.openModels.open(self.target_pdb.text,type="PDB")
		rc("sel #0:." + self.target_chain.text.lower())
		rc("select invert sel")
		rc("delete sel")

	def write_all_SSE_centroids(self):
		"""Infile is the file with starting and ending residues of all SSEs (helices and beta strands).
		The function will calculate the centroids of the SSEs and output them into text file for 1st clustering."""
		
		try:
			infile_name = 'all_SSE_start-end_residues_'+ self.target_pdb.text + "_" + self.target_chain.text +'.txt'
			infile = open(infile_name,'r')
		except IOError:
			print "The file", infile_name, "has not been created."
			print "You need to run the Arun's module first."
			sys.exit()

		outfile = open('all_SSE_centroids_'+ self.target_pdb.text + "_" + self.target_chain.text +'.txt','wb')

		prot = self.model0[0]
		for line in infile:
			line = line.strip() 
			columns = line.split() 
			rng = range(int(columns[0]), int(columns[1])+1)
			the_SSE_coord = []
			for r in rng:
				res = prot.findResidue(r)
				caAtom = res.findAtom('CA')
				the_SSE_coord.append(caAtom.coord())

			the_SSE_coord_arr = np.array(the_SSE_coord)
			for el in np.average(the_SSE_coord_arr,axis=0):
				outfile.write(str(el)+' ')
			outfile.write('\n')
		   
		infile.close()
		outfile.close()

	def superpose_target_against_hits(self):
		""" To superpose the target against a hit of choice to see the quality of SSM superpositon.
		Will report the Q-score of the superposition in the name of the image that will also be saved
		in the corresponding hit's folder."""

		clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
		'#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC' ,
		'#D2691E' , '#DC143C' , '#556B2F' , '#FF1493' , '#98FB98']

		self.open_target_in_Chimera()

		for hit in self.third_filter_pdbs:
		
			model1 = chimera.openModels.open(hit,type="PDB") 
			if len(model1) > 1:    # if more models than one (as in NMR structures) keep only model 1.1      
				for i in range(2,len(model1)+1):
					rc("sel #1."+str(i))
					rc("delete sel") # all others delete
		
			# coloring the matching regions:
			target_residues = []
			a = self.target_matching_res[self.target_matching_res['f0']==hit]
			for i in a:
				rng = range(i['f1'],i['f2']+1)
				for res in rng:
					target_residues.append(res)
					rc("color " + clr[0] + " #0:" + str(res)) 
		  
			hit_residues = self.hits_matching_res[hit]
			for res in hit_residues:
				rc("color " + clr[1] + " #1:" + str(res)) 

			# superimposing them taking the RTMatrix from the file saved above
			my_path = os.path.abspath(os.curdir)+ "/hit_" + hit
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)
			matrix_file = 'RTMatrix_' + self.target_pdb.text + '_' + self.target_chain.text + '_against_' + hit + '.txt'
			
			rc('matrixset '+ matrix_file) 

			rc ("focus")
			png_name = hit + "_superposed_Q-"+ str(self.Q_scores[hit]) + ".png"
			rc("copy file " + png_name + " supersample 3")

			# moving back up the directory 
			os.chdir("..")

			# saving the image again for convenience of viewing them all at once 
			my_path = os.path.abspath(os.curdir)+ "/superpositions" 
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)
			rc("copy file " + png_name + " supersample 3")
			os.chdir("..")
			
			try:
				rc('~open #1')
			except:
				rc('~open #1.1')
			rc('~color') 

	def write_matching_SSE_centroids(self):
		"""The routine will write and save text file with matching SSEs for every hit in
		their corresponding folder. The text file will be used in 1st Clustering."""

		prot = self.model0[0]

		for hit in self.third_filter_pdbs:
			
			my_path = os.path.abspath(os.curdir)+ "/hit_" + hit
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)

			outfile = open('matching_SSE_centroids_'+hit+'.txt','wb')

			for match in self.target_matching_res:
				if match['f0'] == hit:
					temp_list = []
					ran = range(match['f1'],match['f2']+1)
					for r in ran:
						res = prot.findResidue(r)
						caAtom = res.findAtom('CA')
						temp_list.append(caAtom.coord())

					matching_SSE_centroid = np.average(np.array(temp_list),axis=0)
					for cent in matching_SSE_centroid:
						outfile.write(str(cent)+' ')
					outfile.write('\n')

			outfile.close()
			os.chdir("..")

	def get_results_clust_SSE(self):
		"""It will visualize the clustered SSEs for every hit separately and 
		save the images in the corresponding hit's folders.
		It also saves all clusters with more than 2 SSEs in a text file 
		for use as input for 2nd clustering, saving it in target's folder."""

		self.open_target_in_Chimera() 

		clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
		'#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC' ,
		'#D2691E' , '#DC143C' , '#556B2F' , '#FF1493' , '#98FB98']
		
		raw_input("Adjust the view, how you want to be saving your clusters and press any key to continue.")

		all_cluster_count = 1 # will count number of acceptable clusters (> SSE large for all hits)
		for hit in self.third_filter_pdbs:

			my_path = os.path.abspath(os.curdir)+ "/hit_" + hit
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)

			# getting the ids of clustered SSEs
			infile = open('SSEs_clustered_'+ hit +'.txt','r')
			clus = []
			for line in infile:
				line = line.strip()
				columns = line.split()
				clus.append(int(line))
			infile.close()

			clustered_SSE = np.array(clus)

			
			# getting the SSEs residues for hit
			SSEs_res = []
			for match in self.target_matching_res:
				if match['f0'] == hit:
					temp_list = []
					ran = range(match['f1'],match['f2']+1)
					SSEs_res.append(ran)
			SSEs_residues = np.array(SSEs_res)

			# coloring in the clustered SSEs
			rc('~color')
			for i in range(1,np.amax(clustered_SSE)+1): # going through clusters
				for sse in SSEs_residues[clustered_SSE==i]: # clustered_SSE contains the clustering info
					for res in sse:
						rc("color " + clr[i] + " #0:" + str(res))  
			
			# getting all clusters larger than 2 SSEs
			for i in range(1,np.amax(clustered_SSE)+1): # going through clusters
				if len(SSEs_residues[clustered_SSE==i]) > 2: # if the cluster contains more than 1 SSE
					self.all_hits_clusters[all_cluster_count] = []
					for sse in SSEs_residues[clustered_SSE==i]: # clustered_SSE contains the clustering info
						for res in sse:
							self.all_hits_clusters[all_cluster_count].append(res)
					all_cluster_count += 1

			
			png_name = hit + "_first_clust.png"
			rc("copy file " + png_name + " supersample 3")
			os.chdir("..")

			# saving the image again for convenience of viewing them all at once
			my_path = os.path.abspath(os.curdir)+ "/first_clust_results" 
			if not os.path.isdir(my_path):
				os.makedirs(my_path)
			os.chdir(my_path)
			rc("copy file " + png_name + " supersample 3")
			os.chdir("..")

		f = open('all_hits_clusters.pckl', 'w')
		pickle.dump(self.all_hits_clusters, f)
		f.close()


if __name__ == "__main__":
	# choose the xml file
	xml_name = '1f6m_E.xml'

	target = xml_name[:6]
	ssm = SSM_search_results()	
	ssm.xml_file = xml_name
	ssm.parse_xml()
	# creating the target's folder where all the hits' subfolder will be saved
	my_path = os.path.abspath(os.curdir)+ "/target_" + target
	if not os.path.isdir(my_path):
		os.makedirs(my_path)
	os.chdir(my_path)

	print "This is about target " + ssm.target_pdb.text + ", chain " + ssm.target_chain.text
	print "... having", ssm.target_nRes.text, "residues"
	print "... SSM found", ssm.no_of_matches.text, "matching proteins"
	print "Input parameters for lowest acceptable matches:"
	print "- at least " + ssm.target_perc_param.text + "% of the target's SSEs matching at least " + ssm.hit_perc_param.text +"% of the hits' SSEs."

	ssm.remove_redundancy(q_score_cutoff=0.06,name_thresh=0.6,tSSE_cap=50)
	ssm.open_target_in_Chimera()
	ssm.write_RTMatrices()	
	ssm.write_all_SSE_centroids()
	ssm.write_matching_SSE_centroids()
	ssm.superpose_target_against_hits() # can be skipped if no visual inspection of superposition required
										# as it is quite expensive
	
	print "Here run the 'first_clustering.py' from normal python terminal."
	raw_input("Then press any key to see the results of the 1st clustering.")
	ssm.get_results_clust_SSE()
		
	print "The script has finished. Next do the 2nd clustering."
	print "Elapsed time", time.time() - startTime, "seconds."




