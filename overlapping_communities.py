
#Libraries
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_edge_cut,minimum_edge_cut
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,AgglomerativeClustering 
from sklearn.metrics import silhouette_samples, silhouette_score
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from categorical_cluster import hierarchical_mixed
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
             discriminant_analysis, random_projection)
from collections import defaultdict

#Function to get adjacent cliques
def get_adjacent_cliques(clique,node_index):
	clique_adjacent = set()
	for node in clique:
		common_cliques = node_index[node]
		for common_clique in common_cliques:
			if common_clique!=clique:
				clique_adjacent.add(common_clique)

	return clique_adjacent


#Function to get overlapping communities
def get_communities(k,graph):
	#Pre-condition
	if k<2:
		print "Error! K has to be greater than or equal to 2"
		return 

	#Obtain the cliques
	cliques = nx.find_cliques(graph)

	#Select all the cliques of length k, frozenset => Immutable Set 
	cliques = [frozenset(c) for c in cliques if len(c)>=k]

	#Store the information about which nodes are in which cliques
	node_index = defaultdict(list)
    
    #For each node all the corresponding cliques which they are part of are stored
	for clique in cliques:
		for node in clique:
			node_index[node].append(clique)

	#Treat each clique as a node and create a new graph
	perlocation_graph = nx.Graph()

	#Create the nodes
	perlocation_graph.add_nodes_from(cliques)

	#Find the adjacent cliques
	for clique in cliques:
		probable_adjacents = get_adjacent_cliques(clique,node_index)
		for adjacent in probable_adjacents:
			#Now check if any of the adjacent ones have >= k-1 nodes in common -- If yes, combine them
			if len(clique.intersection(adjacent)) >= (k-1):
				perlocation_graph.add_edge(clique,adjacent)

	communities = []

	#Now obtain all the connected components
	for component in nx.connected_components(perlocation_graph):
		communities.append(list(component))


	return communities		


#Function for Test Graph Results
def test():
	#Initialise Test Graph
	test_graph = nx.Graph()

	#Add Edges to the Test Graph
	test_graph.add_edge(1,2,weight=1.0)
	test_graph.add_edge(1,3,weight=1.0)
	test_graph.add_edge(2,3,weight=1.0)
	test_graph.add_edge(3,4,weight=1.0)
	test_graph.add_edge(4,5,weight=1.0)
	test_graph.add_edge(4,6,weight=1.0)
	test_graph.add_edge(5,6,weight=1.0)

	#Write to File
	f_test = open('Tests/overlapping_communities.txt','w')
    
    #Get the communities with k set as 2
	communities = get_communities(2,test_graph)

	for community in communities:
		f_test.write(str(community))


	f_test.close()


def main(k):
	#Loading gene_interactions JSON file into a variable 
	with open('JSON rows/gene_interactions.json') as json_data:
		interactions = json.load(json_data)


	#Information about Graph Connectivity
	graph_info = interactions["results"]

	source = []
	target = []

	#test()
	

	#Creating NetworkX instance
	graph = nx.Graph()

	i = 0
	#Extracting the edge relationships
	for edge in graph_info:
		#temp = []
		source.append(edge[0])
		target.append(edge[2])
		#Don't add a self-loop
		if edge[0] == edge[2]:
			continue

		#Adding the edge in NetworkX
		graph.add_edge(edge[0],edge[2],weight=1.0)


	communities = get_communities(k,graph)

	print communities 



main(3)