import csv
import numpy as np
import pandas as pd
from igraph import *
import igraph as ig
import time
import sys
import argparse


def main():
	# Head of output console
	print("-----------------------------------------------------")
	print("-------------------Program Start---------------------")
	print("-----------------------------------------------------")
	# Record time
	startTime = time.time()

	# Get args
	args = get_args(sys.argv[1:], "CSV to GML script")

	# Set input and output steam
	filename_in = args.filename_in
	filename_out = args.filename_out
	
	# Import data
	toleratio = args.toleratio
	data,data_trimmed = readin(filename_in,toleratio)
	print("Data imported.")

	# Get the names
	dimension_names = [data[1:]] #unused in this program
	node_names = [vector[0] for vector in data][1:]  #remove 'gene_id'

	# Get adjacency matrix
	print("-----------------------------------------------------")
	print("Computing adj_matrix...")
	sim_meas = args.sim_meas
	edge_ratio = args.edge_ratio
	print("Similarity measure: " + sim_meas)
	adjMatrix = getAdjacency(data_trimmed,sim_meas,edge_ratio)
	print("Done.")

	# Now we can convert adjacency matrix to graph
	vertices = node_names
	edges = []
	n = len(adjMatrix)
	for i in range(0,n):
		for j in range(0,n):
			if adjMatrix[i][j]==1:
				edge = tuple([i,j])
				edges.insert(len(edges),edge)
	print("-----------------------------------------------------")
	print("Start Generating Graph...")
	g = Graph(vertex_attrs={"label":vertices}, edges=edges, directed=False)
	print("Graph Generated. There are " + str(len(vertices)) + " vertices and " + str(len(edges)) + " edges.")

	# Export to gml file
	print("-----------------------------------------------------")
	print("Start Writing to GML file...")
	ig.write(g,format="gml",filename=filename_out)
	print("Writing Completed. Output file is named " + filename_out)
	
	# Tail of output console
	# Report time
	print("Finished in %0.4f seconds" % (time.time() - startTime))
	print("-----------------------------------------------------")
	print("------------------Program Finish---------------------")
	print("-----------------------------------------------------")

def get_args(argList, name):

	parser = argparse.ArgumentParser(description=name)
	parser.add_argument("filename_in",type=str,help="File name of input csv.")
	parser.add_argument("filename_out",type=str,help="File name of output gml.")
	parser.add_argument("-t","--toleratio", type=float,action="store",dest="toleratio",help="Ratio of bad entires tolerated.", default=0.1)
	parser.add_argument("-e","--edge_ratio",type=float,action="store",dest="edge_ratio",help="Top percentage of edges.",default=0.1)
	parser.add_argument("-s","--sim_meas",type=str, action="store",dest="sim_meas",choices=["PCC","Euclid"],help="Similarity measure", default="PCC")
	
	return parser.parse_args(argList)


def readin(filename,toleratio):

	# Import from csv
	print("Reading files " + filename)
	data = np.array(pd.read_csv(filename, sep=',',header=None, error_bad_lines=False, index_col=False, dtype='unicode',low_memory=False))

	# Get the shape of data
	rowNum = len(data)
	colNum = len(data[0])
	sample_count = colNum-1

	# Remove genes with bad entries
	flags = np.zeros(rowNum)
	flags[0] = 1
	for i in range(1,len(data)):
		row = data[i]
		zero_count = 0
		for j in range(1,len(row)):
			#print(row[j])
			if not isinstance(float(row[j]),float):
				#print("Not digit")
				zero_count += 1
			else:
				if float(row[j])== 0:
					#print("Zero found")
					zero_count += 1
			row[j] = float(row[j]) # change to float for later computing
			#print("Normal input")
		if zero_count <= toleratio * sample_count:
			flags[i] = 1
	fixed_data = []
	for k in range(0,rowNum):
		if flags[k] == 1:
			fixed_data.append(data[k])

	# Trim the first row and column
	data_trimmed = np.delete(fixed_data,0,0)
	data_trimmed = np.delete(data_trimmed,0,axis=1)

	return fixed_data,data_trimmed

def getAdjacency(data,sim_meas,edge_ratio):

	# Get thbe shape of data
	n = len(data)  #node number
	m = len(data[0])  #dimension number

	# Initialize the adj_matrix
	size = (n,n)
	adj_matrix = np.zeros(size)
	pair_adj_count = int((np.square(n,dtype='int64')-n)/2)
	print("We have " + str(pair_adj_count) + " pairs of adjacencies.")

	# Initialize the adj_set
	adj_set = np.zeros(pair_adj_count)
	adj_index = 0

	# Start computing adj_matrix
	print("Computing adjacencies...")
	for i in range(0,n):
		for j in range(i+1,n):
			if isinstance(data[i], basestring) or isinstance(data[j], basestring):
				adj = 0 # If data is bad, treat as 0 adjacency
			else:
				if sim_meas=='Euclid':
					adj = 1/(1+np.linalg.norm(data[i]-data[j]))
				if sim_meas=='PCC':
					adj = pearson_correlation(data[i],data[j])
			#print(adj)
			adj_matrix[i][j] = adj
			adj_set[adj_index] = adj
			adj_index += 1
	print("Done.")

	# Get threshold
	print("Sorting adjacency list...")
	adj_set_sorted = np.sort(adj_set,kind='mergesort')
	#print(adj_set_sorted)
	print("Done.")
	edge_ratio = edge_ratio # percentage of edges
	print("Take top " + str(edge_ratio) + " adjacency as an edge.")
	adj_thresh = adj_set_sorted[int(pair_adj_count*(1-edge_ratio))]
	print("The threshold of similarity is " + str(adj_thresh))

	# Covert similarity to 0 or 1
	edge_count = 0
	for i in range(0,n):
		for j in range(i+1,n):
			if (adj_matrix[i][j]-adj_thresh)>0:
				adj_matrix[i][j] = 1
				edge_count += 1
			else:
				adj_matrix[i][j] = 0
	print(str(edge_count) + " edges formed.")

	return adj_matrix

def pearson_correlation(object1, object2):

    values = range(len(object1))
    
    # Summation over all attributes for both objects
    sum_object1 = sum([float(object1[i]) for i in values]) 
    sum_object2 = sum([float(object2[i]) for i in values])

    # Sum the squares
    square_sum1 = sum([pow(object1[i],2) for i in values])
    square_sum2 = sum([pow(object2[i],2) for i in values])

    # Add up the products
    product = sum([object1[i]*object2[i] for i in values])

    #Calculate Pearson Correlation score
    numerator = product - (sum_object1*sum_object2/len(object1))
    denominator = ((square_sum1 - pow(sum_object1,2)/len(object1)) * (square_sum2 - 
    	pow(sum_object2,2)/len(object1))) ** 0.5
		
    return (numerator/denominator if denominator != 0 else 0)

if __name__ == '__main__':
	main()
