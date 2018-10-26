# geneGraph

Tool for genome complexity computing

## Dependencies

* Python 3.5 or later
* gene_graph_lib

    You can install gene_graph_lib with `pip install gene_graph_lib`

## Usage

If you have not sif file with graph stucture, the first step is parsing of Orthofinder outputs.
To do this open `source` directory and type in terminal:
` python orthofinder_parse.py -i [path to txt file] -o [path and name prefix for output files] `

Output files:
* **prefix.sif** - all edges list of genome graph
* **prefix_freq.sif** - all edges frequency
* **prefix.db** - database with all parsed informmation
* **prefix_context.sif** - number of contexts, computed for each node in graph

Next step is computing of genome complexity.
To do this type in terminal:

`python start_computing.py -i prefix.sif -o [path to output folder] --reference [name of reference genome]`

Additional parameters are:
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

## Generating of subgraph

Okay, we computed complexity for each gene in our reference genome. Let's suppose that we found some interesting node and we want to observe its context. We developed script which allows us to draw small part of genomes graph. But current ` prefix.sif ` file content information about full graph with too many nodes and edges to draw. 
So, to generate small part of graph you need `generate_subgraph.py` script.

Common usage is

`python generate_subgraph.py -i graph.sif -o subgraph --reference [name of reference genome] --start [name of start node] --end [name of end node]`

This command generate subgraph `subgraph.sif`, built from START......END simple chain of nodes in reference.

Additional parameters:
* ` --window ` - number of nodes, added to left and right side of refernce chain (default 20)
* ` --depth ` - maximum length of derivating path, which will be added to subgraph {default is length of reference chain)
* ` --tails ` - if derivating path too long, it will be replaced by left and right "tails". This parameter is tails length (default 5)
* ` --names_list ` - path to file with names list (default all names from *.sif will be used)

## Graph drawing

Now we can to run drawing of our mini-graph. Let's go to the `recombinatin_draw` directory in geneGraph folder and type in terminal:

`python run_drawing.py -i subgraph.sif -o subgrah_img`

This script generate:
* subgraph_img.ps file as image and 
* subgraph_img.dot file with DOT description of subgraph (DOT is popular graph structure description language)

