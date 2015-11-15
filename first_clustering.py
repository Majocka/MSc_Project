#!/usr/bin/python

import scipy.spatial.distance as dist
import numpy, scipy
import numpy as np
import sys
import os
import random
import pickle

# The script will obtain clusters of matching SSEs per every hit separately
# - to obtain localized structural matches of a hit
# The clusters are partitioned with the distance cutoff - scaled 60% of 
# the average distance between all target's SSEs.
# The scaled percentage was obtained by our linear model, which describes the 
# relationship between the number of target's residues and the 'ideal' scale -
# the scale of 60% of the average distances between all SSEs as manually determined  
# to lead to 'ideal' partitioning of clusters at the domain level

def cluster_SSEs(dist_mat, dist_cutoff):
	"""Will cluster the centroids of matching SSEs per hit.
	Inputs are distance matrix of the matching SSEs centroids and the distance cutoff
	as a scaled percentage of the average distance of all the target's SSEs.
	It will return a dictionary whose values will represent the ids of SSEs 
	that were clustered. Len of the dictionary will say how many clusters are obtained.
	"""

	pool = range(0,len(dist_mat)) # pool of my SSEs to be clustered (from 0 to max - as indexes in Python)
	clust_dict = {} # values of the dictionary will represent the indexes of SSEs being clustered

	while pool:
		num = random.choice(pool)
		clust_dict[num] = [num]  
		pool.remove(num)  
		temp_pool = []
		copy_pool1 = pool[:]  
		for i in copy_pool1:
			my_dist = dist_mat[num][i]
			if my_dist < dist_cutoff:
				print my_dist
				clust_dict[num].append(i)
				pool.remove(i)
				temp_pool.append(i)
			else:
				print 'my distance off', my_dist

		
		copy_pool2 = pool[:] 
		while temp_pool:
			copy_temp_pool = temp_pool[:] 
			for j in copy_temp_pool:
				copy_pool2 = pool[:]  
				for k in copy_pool2:
					my_dist = dist_mat[j][k]
					if my_dist < dist_cutoff:
						print my_dist
						clust_dict[num].append(k)
						pool.remove(k)
						temp_pool.append(k)
					else:
						print 'my distance off', my_dist
				temp_pool.remove(j)
				copy_temp_pool = temp_pool[:] 
	
		pool = list(set(pool).intersection(copy_pool2)) 

	print "Results - values in this dictionaries representing clustered SSEs:"
	c = 0
	for k,v in clust_dict.iteritems():
		c = c + len(v)
		print k,v
	print "Number of clustered SSEs: ", c

	return clust_dict

def get_estimated_scale(x):
	"""To calculate the equation of our linear model
	to predict the scale of 60 percent for calculation of distance cutoff."""
	p1 = -0.00237
	p2 = 1.902
	return p1*x + p2	

########################################## main
# choose the target:
target = '1f6m_E' 

# go to the target's folder
my_path = os.path.abspath(os.curdir)+ "/target_" + target
if not os.path.isdir(my_path):
	os.makedirs(my_path)
os.chdir(my_path)

# getting the number of residues in the target
f = open('nRes_'+ target +'.pckl')
nRes = int(pickle.load(f))
f.close()

# calculating the distance_cutoff for the target
scale = get_estimated_scale(nRes)
perc_cutoff = 0.6 * scale
print target
print scale
print perc_cutoff * 100

# getting the filtered hits for that target
f = open('third_filter_pdbs_'+ target +'.pckl')
hits_list = pickle.load(f)
f.close()

# getting the centroids of all SSEs in the target to calculate
# the distance cutoff
all_SSE_centroids = []
f = open('all_SSE_centroids_'+target+'.txt','r')
for line in f:
	temp = []
	line = line.strip()
	columns = line.split()
	for j in columns:
		temp.append(float(j))
	all_SSE_centroids.append(temp)
f.close()
all_SSE_centroids = numpy.array(all_SSE_centroids)
all_distances = dist.pdist(all_SSE_centroids, 'euclidean')
aver_dist = numpy.average(all_distances,axis=0)
distance_cutoff = perc_cutoff * aver_dist


# looping over the hits to cluster their SSEs independently from hits
for hit in hits_list:
	# getting the matching SSE centroids
	my_path = os.path.abspath(os.curdir)+ "/hit_" + hit
	if not os.path.isdir(my_path):
		print "The hit's folder has not been created by chimera_main.py."
		sys.exit()
	os.chdir(my_path)
	matching_SSE_centroids = []
	f = open('matching_SSE_centroids_'+hit+'.txt','r')
	for line in f:
		temp = []
		line = line.strip()
		columns = line.split()
		for j in columns:
			temp.append(float(j))
		matching_SSE_centroids.append(temp)
	f.close()

	data_mat = numpy.array(matching_SSE_centroids)
	distances = dist.pdist(data_mat, 'euclidean')
	dist_mat = dist.squareform(distances)

	print "\nStarting iterative hierarchical SSE clustering for:"
	print target, "against", hit
	print "With average distance:", aver_dist
	print "... and distance cutoff after scaling:", distance_cutoff

	# doing the clustering
	clusters = cluster_SSEs(dist_mat, distance_cutoff)

	# saving the SSEs clustered as indices in txt file to be loaded in 'main_chimera.py'
	outfile = open('SSEs_clustered_'+ hit +'.txt','wb')
	clust = range(0,len(dist_mat))
	clust_num = 1
	for k,v in clusters.iteritems():
		for el in v:
			clust[el] = clust_num
		clust_num += 1

	for i in clust:
		outfile.write(str(i)+'\n')
	outfile.close()

	print "Created clusters:", max(clust)
	os.chdir("..")

os.chdir("..")

