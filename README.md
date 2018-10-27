# geneGraph

Tool for genome complexity computing

## Dependencies

* Python 3.5 or later
* Graphviz and pygraphviz libraries
* gene_graph_lib

    You can install gene_graph_lib with `pip install gene_graph_lib`

## Usage

### Generating of graph structure

If you don't have sif file with graph stucture, the first step is parsing of Orthofinder outputs.
To do this open `source` directory and type in terminal:
` python orthofinder_parse.py -i [path to txt file] -o [path and name prefix for output files] `

Output files:
* **graph.sif** - all edges list of genome graph
* **graph_freq.sif** - all edges frequency
* **graph.db** - database with all parsed informmation
* **graph_context.sif** - number of contexts, computed for each node in graph


### Complexity computing

Next step is computing of genome complexity.
To do this type in terminal:

`python start_computing.py -i graph.sif -o [path to output folder] --reference [name of reference genome]`

Additional parameters:
* ` --window ` - window size (default 20)
* ` --iterations ` - number of iterations in probabilistic method(default 500)
* ` --genomes_list ` - path to file with names list (default all names from *.sif will be used)
* ` --min_depth, --max_depth ` - minimum and maximum depth of generated paths in graph (default from 0 to inf)
* ` --save_db ` - path to database, created by orthfinder_parse.py (default data dont saved to db, only to txt)

Output files for each contig in reference genome:
* all_bridges_contig_n.txt
* window_variability_contig_n.txt
* main_chain_contig_n.txt
* params.txt
* IO_variability_table_contig_n.txt
* prob_IO_variability_table_contig_n.txt
* prob_window_variability_contig_n.txt
* prob_all_bridges_contig_n.txt

### Generating of subgraph

Okay, we computed complexity for each gene in our reference genome. Let's suppose that we found some interesting node and we want to observe its context. We developed script which allows us to draw small part of genomes graph. But current ` graph.sif ` file content information about full graph with too many nodes and edges to draw. 
So, to generate small part of graph you need `generate_subgraph.py` script.

Common usage is

`python generate_subgraph.py -i graph.sif -o subgraph --reference [name of reference genome] --start [name of start node] --end [name of end node]`

This command generates subgraph `subgraph.sif`, built from START......END simple chain of nodes in reference.

Additional parameters:
* ` --window ` - number of nodes, added to left and right side of refernce chain (default 20)
* ` --depth ` - maximum length of derivating path, which will be added to subgraph {default is length of reference chain)
* ` --tails ` - if derivating path too long, it will be replaced by left and right "tails". This parameter is tails length (default 5)
* ` --names_list ` - path to file with names list (default all names from *.sif will be used)

### Graph drawing

Now we can to run drawing of our mini-graph. Let's go to the `recombinatin_draw` directory in geneGraph folder and type in terminal:

`python run_drawing.py -i subgraph.sif -o subgrah_img`

This script generate:
* subgraph_img.ps file as image and 
* subgraph_img.dot file with DOT description of subgraph (DOT is popular graph structure description language)

Additional parameters:
* ` --freq_min` - minimal edge frequency to draw. Edge frequency is number of genomes with this edge.
* ` --da` - legacy parameter. Draws all edges in any case, but if edge frequency < ` --freq_min`, edge doesnt influence to graph layout.
