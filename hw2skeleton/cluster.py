from .utils import Atom, Residue, ActiveSite
# Some utility classes to represent a PDB structure

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name

import glob
import os



def read_active_sites(dir):
    """
    Read in all of the active sites from the given directory.

    Input: directory
    Output: list of ActiveSite instances
    """
    files = glob.glob(dir + '/*.pdb')

    active_sites = []
    # iterate over each .pdb file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.pdb")):

        active_sites.append(read_active_site(filepath))

    print("Read in %d active sites"%len(active_sites))

    return active_sites


def read_active_site(filepath):
    """
    Read in a single active site given a PDB file

    Input: PDB file path
    Output: ActiveSite instance
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".pdb":
        raise IOError("%s is not a PDB file"%filepath)

    active_site = ActiveSite(name[0])

    r_num = 0

    # open pdb file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in f:
            if line[0:3] != 'TER':
                # read in an atom
                atom_type = line[13:17].strip()
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                atom = Atom(atom_type)
                atom.coords = (x_coord, y_coord, z_coord)

                residue_type = line[17:20]
                residue_number = int(line[23:26])

                # make a new residue if needed
                if residue_number != r_num:
                    residue = Residue(residue_type, residue_number)
                    r_num = residue_number

                # add the atom to the residue
                residue.atoms.append(atom)

            else:  # I've reached a TER card
                active_site.residues.append(residue)

    return active_site


def write_clustering(filename, clusters):
    """
    Write the clustered ActiveSite instances out to a file.

    Input: a filename and a clustering of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusters)):
        out.write("\nCluster %d\n--------------\n" % i)
        for j in range(len(clusters[i])):
            out.write("%s\n" % clusters[i][j])

    out.close()


def write_mult_clusterings(filename, clusterings):
    """
    Write a series of clusterings of ActiveSite instances out to a file.

    Input: a filename and a list of clusterings of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusterings)):
        clusters = clusterings[i]

        for j in range(len(clusters)):
            out.write("\nCluster %d\n------------\n" % j)
            for k in range(len(clusters[j])):
                out.write("%s\n" % clusters[j][k])

    out.close()

example_path="/Users/student/Documents/hw2-skeleton/data"
# test is all 136 active sites 
test=read_active_sites(example_path)

first=test[0]
# avg_coordinates returns a list of xyz cordinates for the backbone of each pdb file.
def avg_coordinates(PDB):
    residues_list=[]
    atoms_list=[]
    for residue in PDB.residues:
        if residue == ['CA','N','C','CB']:
            residues_list.append(residue)
        for atoms in residue.atoms:
            atoms_list.append(atoms.coords)
            
    return atoms_list
# Convert the backbone coordinate list into a dataframe: 
import pandas as pd
import numpy as np
def convert_dataframe(List):
    backbone_coordinates = pd.DataFrame(np.array(List).reshape(len(List),3), columns = list("xyz"))
    return(backbone_coordinates)

def all_active_sites(test):
    for i in test:
        example_list=avg_coordinates(i)
        convert = convert_dataframe(example_list)
        one_active_site=pd.DataFrame(convert.mean(),columns=[i.name]).T
        yield one_active_site
# a dataframe of the meansof all the active sites        
mean_of_all_activesites = pd.concat(list(all_active_sites(test)))

new_array=np.array(mean_of_all_activesites)
# Below I used euclidean distance by pairwise to obtain a similarity meteric

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    similarity = 0.0

    # Fill in your code here!
#Caculating euclidean distance by pairwise: using scipy.spatial distance
from scipy.spatial import distance

#Y = distance.pdist(new_array, 'euclidean')
Y = distance.squareform(distance.pdist(new_array, 'euclidean'))
    return Y


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    
#Randomly intitalize two data points called the center centeroids
#Use Euclidean distacne to find what data point is closest to the centeroids
#Based on the distance from c1 and c2 centeroids, the data point will be grouped into clusters
#Compute the datapoints of the centeroid inside cluster 1 
#Repostion the centeroid of the cluster 1 to the new centeroid
#Compute the centeroid of datapoints inside cluster 2
#Reposition the centeroid of cluster 2 to the new centeroid 
#Repeat the calculation of centeroids and repostioning until none of the cluster assignments change

    return []


def cluster_hierarchically(active_sites):
    
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    #pseudocode step by steo for agglomerative hierarchical clustering: 
        #Compute the proximity matrix 
        #Merge the two closest clusters
        #Update the proximty matrix between new and original clusters 
        #Repeat until one cluster remains
    import numpy as np
    from scipy.spatial import distance_matrix
    from matplotlib import pyplot as plt
    def k_means(new_array,K):
        nrow = new_array.shape[0]
        ncol = new_array.shape[1]
        
        # pick K random data points as intital centroids
        initial_centriods = np.random.choice(nrow, K, replace=False)
        centeroids = new_array[initial_centroids]
        centroids_old = np.zeros[K, ncol]
        cluster_assignments = np.zeros(nrows)
        while (centroids_old != centeroids).any():
            centeroids_old.append(centeroids)
            #compute the distances between data points and the centeriods 
            dist_matrix = distance_matrix(new_array, centeroids, p=2)
            #find closest centeroid
            for i in np.arrange(nrow): 
                d=dist_matrix[i]
                closest_centeroid = (np.where(d == np.min(d)))[0][0]
            #associate data points with closest centeriod
            cluster_assignments[i] = closest_centeroid
            #recompute centeriods
            for k in np.arange[K]:
                new_array_k = new_array[cluster_assignments==k]
                centeroids[k] = np.apply_along_axis(np.mean, axis=0, arr=new_array_k)
        return(centeroids, cluster_assignments)
                

#The similiarity metric I used is called Euclidean distance: I used it because to assign data points to a centeroid I needed a proximity measurement. Euclidean was the best option for me biologically because I can see if the backbone coordinates align, are the same close in distance. 
#In the end when thinking about using the backbone coordinates as a way to measure similarity in active site. This was a decision biologically. biologically it would be hard to overlapped coordinates because all the active sites are not the same length in atoms. For future directions it would be better if I chose either amino acid count or a count postive, negative , hydrophobic phenotype. These measurements would be more informative biologically. 
#I chose K-means clustering
#I chose agglomerative hierarchical clustering 


