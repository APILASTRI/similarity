"""  InterMine @ Open Genome Informatics 
	  -> Implementation of Markov Clustering to obtain the clusters in Gene Interaction Network              """
	  

#Libraries
from __future__ import division
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import math
import string
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.cluster import KMeans,AgglomerativeClustering, MiniBatchKMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import normalize
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
			 discriminant_analysis, random_projection)

""" Algorithm Steps : Markov Clustering

	 1. Find Adjacency Matrix
	 2. Add Self Self-Loops
	 3. Normalize the Matrix
	 4. Do matrix multiplication according to the given power value
	 5. Do the process of inflation for each column
	 6. Repeat Steps 4 and 5 until convergence 
	 7. Inspect the Matrix to get clusters                                         """


#Function to find holes in the feature set for the genes
def MCL(graph,inflation,e):
	#Nodes of the Graph
	nodes = graph.nodes()

	#Convert into adjacency matrix
	adjacency = nx.adjacency_matrix(graph)

	adjacency = adjacency.todense()

	#Normalization of Matrix
	normalized_matrix = adjacency / adjacency.sum(axis=0)	
    
    """ Step for Convergence  """
	#Iterate until convergence is reached
	for j in range(0,2):
		
		""" Power Operation """
		#Matrix multiplication
		for i in range(0,e-1):
			normalized_matrix = np.matmul(normalized_matrix,normalized_matrix)

		#Create Empty Numpy array
		new_array = np.zeros(shape=(len(adjacency),len(adjacency)))

		""" Inflation Operation """
		#Inflate each column with the given power
		for i in range(0,len(new_array)):
			new_array[i] = np.power(normalized_matrix[i],inflation)

		#Normalize the matrix
		normalized_matrix = new_array / new_array.sum(axis=0)
		
		

	
    
	return normalized_matrix		




def main():
	#Loading gene interactions JSON file into a variable 
	with open('JSON rows/gene_interactions.json') as json_data:
		interactions = json.load(json_data)


	#Information about Graph Connectivity
	graph_info = interactions["results"]

	#Creating NetworkX instance
	graph = nx.Graph()

	i = 0
	#Extracting the edge relationships
	for edge in graph_info:
		#Adding the edge in NetworkX
		graph.add_edge(edge[0],edge[2])
		if i == 5:
			break
		i += 1

	inflation = 2
	e = 2


	test_graph = nx.Graph()
	test_graph.add_edge(1,2)
	test_graph.add_edge(1,4)
	test_graph.add_edge(1,3)
	test_graph.add_edge(2,4)
	test_graph.add_edge(1,1)
	test_graph.add_edge(2,2)
	test_graph.add_edge(3,3)
	test_graph.add_edge(4,4)

	#Perform MCL Algorithm
	MCL_matrix = MCL(test_graph,inflation,e)	

	#print MCL_matrix




main()