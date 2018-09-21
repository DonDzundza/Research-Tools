# csv_to_gml.py

This tool uses python-igraph to construct a graph with network data, and exports the graph in "gml" format. 

## Usage

```
[Packages]

pandas,numpy,igraph

[Command line]

$python csv_to_gml.py [-h] [-t TOLERATIO] [-e EDGE_RATIO] [-s {PCC,Euclid}] filename_in filename_out

[Arguments]

Must-have arguments:     

    filename_in                 File name of input csv
    filename_out                File name of output gml

Optional arguments:

    -t,--toleratio              There could be some bad entries in the input file. 
                                We need to remove nodes with too many bad entries and keep the rest 
                                that we can tolerate(the percentage of bad entries is lower than toleratio).
                                Default value for this argument is 0.1
                                
    -e,--edge_ratio             When we build the graph from adjacency matrix, we can pick top 10% (by 
                                default) of the adjacencies and treat them as edges.
    
    -s,--sim_meas               Two similarity measurements can be chosen from. 
                                "Euclid" --> Euclidean Similarity (1/1+euclidean_dist)
                                "PCC" --> Pearson Correlation
```

## Example

```
You can find "testBlock.csv" in this directory. It is a sample dataset of gene expression network that has 
been normalized already. The output "testBlock.gml" is a graph that can be used for community detection.

[Command line]

$python csv_to_gml.py testBlock.csv testBlock.gml

[Console Output]

-----------------------------------------------------
-------------------Program Start---------------------
-----------------------------------------------------
Reading files testBlock.csv
Data imported.
-----------------------------------------------------
Computing adj_matrix...
Similarity measure: PCC
We have 91 pairs of adjacencies.
Computing adjacencies...
Done.
Sorting adjacency list...
Done.
Take top 0.1 adjacency as an edge.
The threshold of similarity is 0.403995643101493
9 edges formed.
Done.
-----------------------------------------------------
Start Generating Graph...
Graph Generated. There are 14 vertices and 9 edges.
-----------------------------------------------------
Start Writing to GML file...
Writing Completed. Output file is named testBlock.gml
Finished in 0.0127 seconds
-----------------------------------------------------
------------------Program Finish---------------------
-----------------------------------------------------

Output file is testBlock.gml
```
