#!/usr/bin/python

import os
import random
import chimera
from chimera import runCommand as rc 
import pickle
import numpy as np

# The script clusters 1st clusters on the basis of their overlap
# to generate final clusters - extracted evolutionary information
# if they overlap by more than 20% they become members of the same 
# final cluster

rc("~open all") 

# choose the target:
target = '1f6m_E'  

# the optimal value of 20% of the overlap between the 1st clusters 
perc_cutoff = 0.2 

# go to the target's folder
my_path = os.path.abspath(os.curdir)+ "/target_" + target
if not os.path.isdir(my_path):
	os.makedirs(my_path)
os.chdir(my_path)

# access the hits' 1st clusters that have more than 2 SSEs
f = open('all_hits_clusters.pckl')
hits_clusters = pickle.load(f)
f.close()

# iterative hierarchical clustering
pool = hits_clusters.keys() 
max_overlap_dict = {}
while pool:
	try: 
		num = random.choice(range(0,len(pool)-1))
	except IndexError:
		num = 0
	pick = pool[num]
	max_overlap_dict[pick] = [pick]  
	pool.remove(pick)  
	temp_pool = []
	copy_pool = pool[:]  
	for i in copy_pool:
		cutoff = (len(hits_clusters[pick]) + len(hits_clusters[i])) * perc_cutoff
		inters = set(hits_clusters[pick]).intersection(hits_clusters[i]) 
		if len(inters) >= cutoff:
			max_overlap_dict[pick].append(i)
			pool.remove(i)
			temp_pool.append(i)

	copy2_pool = pool[:]
	while temp_pool:
		copy_temp_pool = temp_pool[:] 
		for j in copy_temp_pool:
			copy2_pool = pool[:]  
			for k in copy2_pool:
				cutoff = (len(hits_clusters[j]) + len(hits_clusters[k])) * perc_cutoff
				inters = set(hits_clusters[j]).intersection(hits_clusters[k]) 
				if len(inters) >= cutoff:
					max_overlap_dict[pick].append(k)
					pool.remove(k)
					temp_pool.append(k)
			temp_pool.remove(j)

			copy_temp_pool = temp_pool[:] 	    
	pool = list(set(pool).intersection(copy2_pool)) 

print "\n", target
print "There are ", len(max_overlap_dict), " FINAL clusters:"
print "Each being covered by:"
for k,v in max_overlap_dict.iteritems():
	print len(v), "first clusters"
print

# seeing the different clusters in chimera:
model0 = chimera.openModels.open(target[:4],type="PDB")
rc("sel #0:." + target[-1].lower())
rc("select invert sel")
rc("delete sel")

clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
    '#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC' ,
    '#D2691E' , '#DC143C' , '#556B2F' , '#FF1493' , '#98FB98']

count = 0
cluster_residues = {}
for overlap in max_overlap_dict.values():
	cluster_residues[count] = []
	for hit in overlap:
		for res in hits_clusters[hit]:
			rc("color " + clr[count] + " #0:" + str(res))
			if res not in cluster_residues[count]:
				cluster_residues[count].append(res)
	count += 1
	if count < len(max_overlap_dict):
		raw_input('Want to see another cluster on top of it? Click any key to continue.')
	
raw_input("Adjust the view to save the image in the target's folder.")
png_name = "FINAL_overlap_clusters_"+ target +".png"
rc("copy file " + png_name + " supersample 3")

# getting the centroids of the final clusters and writing them in text file  
# to be used in comparison against the closest RIBFIND's RBs
prot = model0[0]
f = open('final_clusters_centroids_'+target+'.txt','wb')
for k,v in cluster_residues.iteritems():
	ca_coord = []
	for r in v:
	    res = prot.findResidue(r)
	    caAtom = res.findAtom('CA')
	    ca_coord.append(caAtom.coord())
	ca_coord_arr = np.array(ca_coord)
	centroid = np.average(ca_coord_arr,axis=0)

	for i in centroid:
		f.write(str(i) + ' ')
	f.write('\n')
f.close()

os.chdir("..")

