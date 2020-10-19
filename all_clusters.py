#######################################################################
# finding clusters for atomistic empty DPC micelle systems            #
#                                                                     #
# JJH 1.23.20                                                         #
#######################################################################

#####################
# loading libraries # 
#####################

import MDAnalysis as mda 
import MDAnalysis.analysis.distances 
import numpy as np
import scipy.cluster
from scipy.cluster.hierarchy import single, fcluster 
import matplotlib.pyplot as plt
import MDAnalysis.core.topologyattrs
import MDAnalysis.analysis.rdf 
import math 

###################
#  input options  #
###################

u = mda.Universe('input.gro', 'input.xtc') # specify input trajectory
cutoff = 8.5 # angstroms, cutoff distance for clustering 
cS = u.select_atoms('name C312') # last carbon in tail of surfactant
save = 'save_directory/clusters.npy' # specify where to save output

##############
# clustering #
##############

box = u.dimensions
r = [] # get positions 
for ts in u.trajectory:
    r.append((cS.positions))
r = np.array(r)
dist = [] # get all pair distances, use box definition to apply minimum image distance
for i in range(len(r)):
    dist.append(MDAnalysis.analysis.distances.self_distance_array(r[i], box=box))
dist = np.array(dist)
print('clustering.')
z = [] # perform hierarchical single-linkage clustering 
for i in range(len(dist)):
    z.append(single(dist[i]))
z = np.array(z)
print('clustering..')
hierarchy = [] # get clusters using cutoff (in angstroms)
for i in range(len(z)):
    hierarchy.append(fcluster(z[i], cutoff, criterion='distance'))
hierarchy = np.array(hierarchy) 
print('clustering...')
uni = [] # get counts of each unique cluster 
for i in range(len(hierarchy)):
    uni.append(np.unique(hierarchy[i], return_counts=True))
uni = np.array(uni)
clusters = [] # slice out one column so we only have counts
for i in range(len(uni)):
    clusters.append(uni[i, 1])
clusters = np.array(clusters)
print('saving...')
np.save(save, clusters)
print('clusters saved as '+save+'!')


#######################################################################################

