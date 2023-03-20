#this script requires functions graph_density.py and Betti_k.py from TDA tutorial (Centeno et al., 2022)
    #modified version of graph_density can be found below (line 45)
#see Centeno, E.G.Z., Moreni, G., Vriend, C., Douw, L., & Santos, F.A.N. (2022). A hands-on tutorial on network and topological neuroscience. Brain Structure and Function, 227(3), 741-762.
import numpy as np # version 1.18.5 146
import networkx as nx # version 2.4 147
import community # version 0.13 (python-louvain) 148
import gudhi # version 3.3.0 149
import scipy.io # version 1.4.1 150
import os
import math
import graph_density as gd
import Betti_k as bet
import itertools

os.chdir('input_directory') #indicate directory where data.mat file is stored
path = os.getcwd()
mat = scipy.io.loadmat('data.mat')  #.mat file where matrix variable is stored

mat = mat['data'] #same as filename if matrix variable name was saved with the same name in matlab

matlength = mat.shape #gives shape of matrix variable (34716 x 163 x 12) 34716 edges, 163 subjects, 12 tasks
densval = np.arange(0, 1, 0.01) #set the filtration value (0-1 range, with 0.01 stepwise increments)
betti = 0 #set ith homology (0 for 0 dimension, 1 for 1 dimension)

output_mat = np.zeros((264, 163, 3), dtype=float, order='F') #output of TDA

for z in range(matlength[2]):
    for s in range(matlength[1]):
        temp = basemat[:, s, z]
        tempr = []
        # convert Z back to r-coefficients
        for i in range(len(temp)):
            tempr.append((math.exp(2 * temp[i]) - 1)/(1 + math.exp(2 * temp[i])))
        tempr = scipy.spatial.distance.squareform(tempr, checks=True)
        #tempr = abs(tempr) #You can choose to remove the absolute value and keep negative values in

        for dens in range(len(densval)):
            G_tda = gd.graph_density(densval[dens], tempr)
            B = bet.Betti_k(G_tda, betti, verbose=False)
            output_base[dens, s, z] = B
            if B < 2:
                break
            scipy.io.savemat('output_mat.mat', {'output_mat': output_mat}, appendmat=True)

#######################################################################################################################
#graph_density adapted from Centeno et al. (2022)
import numpy as np
import networkx as nx
def graph_density(d, matrix):
    np.fill_diagonal(matrix,0)
    # Flatten and rank the correlation values
    temp = sorted(matrix.ravel(), reverse=True)
    size = len(matrix)
    cutoff = np.ceil(d*(size*(size-1)))
    value = temp[int(cutoff)]
    #value = 1-d #this can be used instead of "value" above to filter by correlation threshold (1-Pearson's r)
    G0 = nx.from_numpy_matrix(matrix)
    G0.remove_edges_from(list(nx.selfloop_edges(G0)))
    G1 = nx.from_numpy_matrix(matrix)

    for u,v,a in G0.edges(data=True):
        if (a.get('weight')) <= value:
            G1.remove_edge(u, v)

    finaldensity = nx.density(G1)

    return G1