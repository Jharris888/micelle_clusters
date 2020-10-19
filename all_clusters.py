#######################################################################
# finding clusters for 5 reps of atomistic empty DPC micelle systems  #
# analysis for SANS_FSE                                               #
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

##############
# clustering #
##############

def_clust = 10 # the size of a cluster 

for i in range(1,6):
    print('loading trajectory...')
    u = mda.Universe('rep%i/tric.gro'%i, 'rep%i/tric.xtc'%i)
    save = 'Analysis/clusters/rep%i_cluster.npy'%i 
    box = u.dimensions
    cutoff = 8.5 # angstroms, based on second peak of rdf graph 
    cS = u.select_atoms('name C312') # select the atoms
    r = [] # get positions 
    for ts in u.trajectory:
        r.append((cS.positions))
        #print('loading positions frame %i of %i'%(u.trajectory.frame+1, len(u.trajectory)))
    r = np.array(r)
    dist = [] # get all pair distances formatted as flat upper triangles
    for i in range(len(r)):
        dist.append(MDAnalysis.analysis.distances.self_distance_array(r[i], box=box))
        #print('loading distances frame %i of %i'%(i+1, len(u.trajectory)))
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

# load the cluster data and calculate mean and std of cluster sizes 

import MDAnalysis as mda 
import MDAnalysis.analysis.distances 
import numpy as np
import scipy.cluster
from scipy.cluster.hierarchy import single, fcluster 
import matplotlib.pyplot as plt
import MDAnalysis.core.topologyattrs
import MDAnalysis.analysis.rdf 
import math 

##############
# clustering #
##############

def_clust = 10 # the size of a cluster 
last_frames=100

all_clusts = []
real_clusts = []
numclusts = []
for k in range(1,6):
    all_clust = []
    real_clust = []
    numclust = []
    clusters = np.load('Analysis/clusters/rep%i_cluster.npy'%k,allow_pickle=True)  
    for i in range(len(clusters)):
        numclust1 = []
        all_clust.append(clusters[i])
        for j in range(len(clusters[i])):
            if clusters[i][j] > def_clust:
                numclust1.append(clusters[i][j])
        numclust.append(len(numclust1))
        real_clust.append(np.array(numclust1))
    numclusts.append(np.array(numclust))
    real_clusts.append(real_clust)
    all_clusts.append(np.array(all_clust))

# appending only the last last_frames frames of each rep and concatenating 

end_clusts = []
histdatas = []
for j in range(0,5):
    histdata = []    
    end_clust = []
    for i in range(len(real_clusts[j])-last_frames,len(real_clusts[j])):
        histdata.append(all_clusts[j][i])
        end_clust.append(real_clusts[j][i])
    end_clusts.append(np.concatenate(end_clust))
    histdatas.append(np.concatenate(histdata))

# mean of clusters bigger than def_clust 

means = []
for i in range(0,5):
    means.append(np.mean(end_clusts[i]))


###############################
# graphing and saving figures #
###############################

# convergence graphs for each replicate


plt.figure()
plt.plot(numclusts[0], 'r')
plt.plot(numclusts[1], 'gold')
plt.plot(numclusts[2], 'green')
plt.plot(numclusts[3], 'blue')
plt.plot(numclusts[4], 'purple')
plt.xlabel("Time (ns)")
plt.ylabel("Number of Clusters")
plt.minorticks_on()
#plt.savefig('figures/conv/conv_reps.pdf')
plt.show()

# histograms for each replicate

maximum = np.amax(np.concatenate(histdatas))
binedges = np.linspace(-0.5,maximum+0.5,maximum+2)
bincenters = binedges[1:] - ((binedges[1] - binedges[0])/2.)

for i in range(0,5):
    hist, e = np.histogram(histdatas[i], bins=binedges)
    normhist = hist / np.sum(hist)
    plt.figure()
    plt.bar(bincenters, normhist, width=1.0, color='k')
    plt.ylabel('p(Cluster Size)')
    plt.xlabel('Cluster Size')
    plt.xlim(0,maximum+10)
    plt.title('cluster size (greater than 15) mean = %.2f'%means[i])
    plt.minorticks_on()
    plt.savefig('figures/hists/hist_rep%s.pdf'%i)
    plt.show()

# average conv. graph

convs = []
for i in range(0,5):
    conv1 = []
    for j in range(0,2000):
        conv1.append(numclusts[i][j])
    convs.append(conv1)

av_conv = np.mean(convs, axis=0)
std_conv = np.std(convs, axis=0)

plt.figure()
plt.plot(av_conv, 'b')
plt.errorbar(np.arange(len(av_conv)), av_conv, yerr=std_conv, fmt='none', elinewidth=1, ecolor='c')
plt.xlabel("Time (ns)")
plt.ylabel("Number of Clusters")
plt.minorticks_on()
plt.savefig('figures/av_conv.pdf')
plt.show()


# average histogram 

av_hist = np.concatenate(histdatas)
hist, e = np.histogram(av_hist, bins=binedges)
normhist = hist / np.sum(hist)
plt.figure()
plt.bar(bincenters, normhist, width=1.0, color='b')
plt.ylabel('p(Cluster Size)')
plt.xlabel('Cluster Size')
plt.xlim(0,maximum+10)
plt.title(r'cluster size (greater than 15) mean = %.2f $\pm %.2f$'%(np.mean(np.concatenate(end_clusts)), np.std(np.concatenate(end_clusts))))
plt.minorticks_on()
#plt.savefig('figures/av_hist.pdf')
plt.show()