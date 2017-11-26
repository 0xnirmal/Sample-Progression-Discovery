import numpy as np
import pandas as pd
import sys
import argparse
import seaborn as sns
import matplotlib.pyplot as plt 
plt.style.use("ggplot")
import math
import scipy
import scipy.io as sio
from collections import defaultdict
import matplotlib as mpl
# import networkx as nx
# import visJS2jupyter.visJS_module
import copy


# Argument and data parsing (processing mat files in python proved to be more difficult than originally expected)
# 
# I subsetted the original dataframe with the genes found using the feature selection from the original source (I worked backwards--first I tried to reconstruct the progression using the subsetting genes, then I tried to implement the gene subsetting)


parser = argparse.ArgumentParser(description='Process display arguments')
parser.add_argument("-f", "--jupyter-json")
parser.add_argument("-mat-file", "--mat-file", default="../data/full_data.mat")
parser.add_argument("-gene-subset-file", "--gene-subset-file", default="../data/gene_names.txt")
args = parser.parse_args()

gene_subset = pd.read_csv(args.gene_subset_file, header=None, squeeze=True)
mat = sio.loadmat(args.mat_file)

all_genes = []
for x in mat["probe_names"]:
    if not isinstance(x[0][0], np.unicode_):
         all_genes.append(str(x[0][0][0]))
    else:
        all_genes.append(str(x[0][0]))

full_df = pd.DataFrame(mat["data"], index=all_genes, columns=[x[0] for x in mat["exp_names"][0]])
df = full_df.loc[gene_subset]


# Given the full data frame (genes x time_steps), construct an adjacency matrix where edges are the euclidean distance of one time step with another


def get_adj_mat_from_df(df):
    adj_mat = pd.DataFrame([], index=df.columns, columns=df.columns)
    for i in df:
        for j in df:
            adj_mat.loc[i, j] = float(scipy.spatial.distance.euclidean(df[i], df[j]))
    return adj_mat

adj_mat = get_adj_mat_from_df(df)


# Boruvka's algorithm from the paper to find Minimum Spanning- Tree of a given connected, undirected and weighted graph
# 
# Taken from GeeksforGeeks (http://www.geeksforgeeks.org/greedy-algorithms-set-9-boruvkas-algorithm/) w/ minor modifications 


class Graph:

   def __init__(self,vertices):
       self.V= vertices #No. of vertices
       self.graph = [] # default dictionary to store graph
       self.output_edges = [] #hold the outputs
        
   # function to add an edge to graph
   def addEdge(self,u,v,w):
       self.graph.append([u,v,w])

   # A utility function to find set of an element i
   # (uses path compression technique)
   def find(self, parent, i):
       if parent[i] == i:
           return i
       return self.find(parent, parent[i])

   # A function that does union of two sets of x and y
   # (uses union by rank)
   def union(self, parent, rank, x, y):
       xroot = self.find(parent, x)
       yroot = self.find(parent, y)

       # Attach smaller rank tree under root of high rank tree
       # (Union by Rank)
       if rank[xroot] < rank[yroot]:
           parent[xroot] = yroot
       elif rank[xroot] > rank[yroot]:
           parent[yroot] = xroot
       #If ranks are same, then make one as root and increment
       # its rank by one
       else :
           parent[yroot] = xroot
           rank[xroot] += 1

   # The main function to construct MST using Kruskal's algorithm
   def boruvkaMST(self):
       parent = []; rank = []; 

       # An array to store index of the cheapest edge of
       # subset. It store [u,v,w] for each component
       cheapest =[]

       # Initially there are V different trees.
       # Finally there will be one tree that will be MST
       numTrees = self.V
       MSTweight = 0

       # Create V subsets with single elements
       for node in range(self.V):
           parent.append(node)
           rank.append(0)
           cheapest =[-1] * self.V
    
       # Keep combining components (or sets) until all
       # compnentes are not combined into single MST

       while numTrees > 1:

           # Traverse through all edges and update
           # cheapest of every component
           for i in range(len(self.graph)):

               # Find components (or sets) of two corners
               # of current edge
               u,v,w =  self.graph[i]
               set1 = self.find(parent, u)
               set2 = self.find(parent ,v)

               # If two corners of current edge belong to
               # same set, ignore current edge. Else check if 
               # current edge is closer to previous
               # cheapest edges of set1 and set2
               if set1 != set2:    
                    
                   if cheapest[set1] == -1 or cheapest[set1][2] > w :
                       cheapest[set1] = [u,v,w] 

                   if cheapest[set2] == -1 or cheapest[set2][2] > w :
                       cheapest[set2] = [u,v,w]

           # Consider the above picked cheapest edges and add them
           # to MST
           for node in range(self.V):

               #Check if cheapest for current set exists
               if cheapest[node] != -1:
                   u,v,w = cheapest[node]
                   set1 = self.find(parent, u)
                   set2 = self.find(parent ,v)

                   if set1 != set2 :
                       MSTweight += w
                       self.union(parent, rank, set1, set2)
                       print ("Edge %d-%d with weight %f included in MST" % (u,v,w))
                       self.output_edges.append((u, v))
                       numTrees = numTrees - 1
            
           cheapest =[-1] * self.V


g = Graph(adj_mat.shape[0])
for i in range(adj_mat.shape[0]):
   for j in range(adj_mat.shape[1]):
       if j > i:
           g.addEdge(i, j, adj_mat.iloc[i, j])

g.boruvkaMST()


# G = nx.complete_graph(17)
# nodes = G.nodes()
# edges = [(0, 1),
#          (1, 2),
#          (1, 3),
#          (3, 4),
#          (4, 5),
#          (5, 6),
#          (6, 7),
#          (7, 8),
#          (8, 9),
#          (9, 10),
#          (10, 11),
#          (11, 12),
#          (12, 13),
#          (13, 14),
#          (14, 16),
#          (16, 15)
#         ]

# pos = nx.circular_layout(G)
# nodes_dict = [{"id":n, "x":pos[n][0]*200, "y":pos[n][1]*200} for n in nodes]

# node_map = dict(zip(nodes,range(len(nodes))))  
# edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]], 
#               "title":'test'} for i in range(len(edges))]

# visJS2jupyter.visJS_module.visjs_network(nodes_dict, edges_dict, graph_title="Paper Results", graph_id=0)


# G = nx.complete_graph(17)
# nodes = G.nodes()
# edges = G.edges()
# edges = g.output_edges

# pos = nx.circular_layout(G)
# nodes_dict = [{"id":n, "x":pos[n][0]*200, "y":pos[n][1]*200} for n in nodes]

# node_map = dict(zip(nodes,range(len(nodes))))  
# edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]], 
#               "title":'test'} for i in range(len(edges))]

# visJS2jupyter.visJS_module.visjs_network(nodes_dict,edges_dict, graph_title="My Results", graph_id=1)


# My results replicate the celluar progression process well--there is full continuity from time step 0 to 14 with a jump to 16 after 14 and then a reversion to 15. This is nearly identical to the results in the provided code--this is a slight disparity in which the provided code generates edges from 1 to 2 and then from 1 to 3 instead of directly from 1 to 2 and then from 2 to 3. 
# 
# I tried to understand why the MSTs generated were different, but couldn't seem to identify the difference between the two methods--both my code and the provided source code use Boruvska's algorithm using the selected genes with Euclidean distance as the edge weight to generate the MST. 

# Now I try to recreate the feature selection--there were 54 gene modules created using the source code provided:

# In[626]:

import sklearn.cluster
L = 10
c1 = 0.7
c2 = 0.9


def is_coherent(df):
    if df.shape[0] != 17:
        raise TypeError("Invalid Shape!")
    c = df.corr(method="pearson").mean().mean()
    return c > c1

def split_module(module):
    init_clustering_df = pd.DataFrame([], index=module.index, columns=range(L))
    for i in range(L):
        kmeans = sklearn.cluster.KMeans(n_clusters=2).fit(module)
        init_clustering_df[i] = kmeans.labels_
    
    kmeans = sklearn.cluster.KMeans(n_clusters=2).fit(init_clustering_df)
    mask = (kmeans.labels_ == 1)
    return (module[mask], module[~mask])

module_list = []
init_modules = split_module(full_df)
module_list.append(init_modules[0])
module_list.append(init_modules[1])

stop_condition = False
while not stop_condition:
    new_module_list = []
    stop_condition = True
    for module in module_list:
        if is_coherent(module.transpose()):
            new_module_list.append(module)
        else:
            stop_condition = False
            splits = split_module(module)
            new_module_list.append(splits[0])
            new_module_list.append(splits[1])
    if not stop_condition:
        module_list = new_module_list



def find_first_match(module_list):
    for i in range(len(module_list)):
        for j in range(i+1, len(module_list)):
            c = pd.DataFrame([module_list[i].mean(), module_list[j].mean()]).transpose().corr().loc[0, 1]
            if c > c2:
                return (i, j)
    return (-1, -1)

def find_first_match_revised(module_list):
    df = pd.DataFrame([], columns=range(len(module_list)), index=module_list[0].columns)
    for i in range(len(module_list)):
        df[i] = module_list[i].mean()
    corr_gram = df.corr()
    for i in range(len(module_list)):
        for j in range(i+1, len(module_list)):
            if corr_gram.iloc[i, j] > c2:
                return (i, j)
    return (-1, -1)

comb_module_list = copy.deepcopy(module_list)
stop_condition = False
while not stop_condition:
    stop_condition = True
    match = find_first_match_revised(comb_module_list)
    if match[0] != -1:
        stop_condition = False
        new_comb_module_list = []
        for i in range(len(comb_module_list)):
            if i != match[0] and i != match[1]:
                new_comb_module_list.append(comb_module_list[i])
        concat = pd.concat([module_list[match[0]], module_list[match[1]]])
        new_comb_module_list.append(concat)
        comb_module_list = new_comb_module_list


print("Number of Modules:")
print(len(comb_module_list))

# Using my code, there are 627 "coherent" modules generated from the data.
# This is extremely different from the 54 gene modules created using the provided source code, so there must be a bug in my code when I try to recombine the modules. 
# I think this may be due to the ambiguity with how the paper worded the recombination process of modules: 
# "If the Pearson correlation of two modules' centers is higher than a pre-specified thershold c2 these two modules are merged"
# I assumed this meant to take the means of modules at each time step and then find the correlation among the means of two modules, but this is slightly ambiguous. It is also possible I have a bug in my code in this step; however, my four hours were up after spending 20 minutes attempting to debug this.
