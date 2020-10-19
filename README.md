# micelle_clusters
 Python code for converting trajectory information of a system of micelles into lists of the number of monomers in each cluster in each time step. The code uses MDAnalysis to load trajectories, select atoms, and calculate distances between atoms. It uses Scipy to perform hierarchical clustering on the distances with a defined cutoff. The output is return as a list of lists of each frame in the trajectory, containing the aggregation number of each cluster present in each frame.
 
 Input options:
 
 1.  coordinate and trajectory file in the line 'u = mda.Universe('input.gro', 'input.xtc')' 
 
 2. cutoff, a distance in angstroms which corresponds to the maximum distance for a surfactant to be considered part of the hierarchical cluster. 
 
 3. the name of the atom which will be used for clustering, this is the point on the surfactant which must be in range with the same point on other surfactants in order for the two surfactants to be clustered together, 'cS = u.select_atoms('name C312')', in this example, 'name C312 selects the last carbon tail bead in  
 
 4. save, defines where to write the final output file of clusters. 