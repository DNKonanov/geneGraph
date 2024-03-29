# geneGraph

Command-line tool for genome complexity computing

## Dependencies

* Python 3.4 or later
* Graphviz and pygraphviz libraries (required only for subgraph plotting)
    `conda install -c anaconda graphviz`
    `conda install -c conda-forge pygraphviz==1.7`

* python3 modules
    `pip3 install -r requirements.txt`

## Usage

The most simple way to generate graph structure and estimate complexity is `gg.py` script. The first step is the inference of orthology groups, which is described [here](https://github.com/paraslonic/orthosnake). Output is `OrthoGroups.txt` file.

Now we can use `gg.py`:
`python gg.py -i [path to OrthoGroups.txt] -o [output dir]`

Output dir will be created automatically.
NB! This command calculates complexity profiles for ALL genomes in the dataset, and may spend a lot of time. To specify only one genome you should use `--reference` parameter.

Other parameters:
* ` --window ` - sliding window size (default 20)
* ` --iterations ` - number of iterations in probabilistic method(default 500)
* ` --genomes_list ` - path to file with names list (default all names from *.sif will be used)
* ` --min_depth, --max_depth ` - minimum and maximum depth of generated paths in graph (default from 0 to inf)
* ` --coalign ` - True/False value. If True, enable maximization of genomes collinerity. Quadratically increses analysis time. Default is False. 

Next sections describes parsing and complexity computing separately.

### Generating of graph structure

If you don't have a sif file with graph stucture, the first step is parsing of Orthofinder outputs.
To do this run in terminal:
` python parse_og.py -i [path to txt file with orthogroups] -o [path and name prefix for output files] `

For example:
`python parse_og.py -i ~/data/Mycoplasma/Results/Orthogroups.txt -o ~/data/outputs/Mycoplasma/graph`

Output files:
* **graph.sif** - all edges list of the genomes graph
* **graph.db** - SQLite database with all parsed informmation
* **graph_context.sif** - number of unique contexts, computed for each node in graph
* **graph_genes.sif** - list of all genes (nodes) from all genomes, with coordinates and Prokka annotations



### Complexity computing

To compute complexity type in terminal:

`python estimate_complexity.py -i graph.sif -o [path to output folder] --reference [name of reference genome]`

Additional parameters:
* ` --window ` - sliding window size (default 20)
* ` --iterations ` - number of iterations in probabilistic method(default 500)
* ` --genomes_list ` - path to file with names list (default all names from *.sif will be used)
* ` --min_depth, --max_depth ` - minimum and maximum depth of generated paths in graph (default from 0 to inf)
* ` --save_db ` - path to database, created by orthfinder_parse.py (default data will not be saved to db, only to txt)

Output folder contains number of files:
* `window_complexity_contig_*.txt` - by strain calculated complexity profile
* `prob_window_complexity_contig_*.txt` - complexity profile calculated by probabilistic approach (recommended to use)

Some additional results or information is located in `extended_info` subfolder:

* `extended_info/(prob_)main_chain_contig_*.txt` - simple chain of nodes in the reference
* `extended_info/(prob_)IO_complexity_contig_*.txt` - for each node in the reference chain a number of distinct deviating paths which started from this node or returned to this node (ignore window parameter)
* `extended_info/(prob_)all_bridges_contig_*.txt` - number af deviating paths between all pairs of nodes from the reference chain(ignore window paramter)

### Generating of subgraph

Okay, we computed complexity for each gene in the reference genome. Let's suppose that we found some interesting node and we want to observe its context. We developed script which allows us to draw genomes graph structure. But current `graph.sif` file content information about full graph with too many nodes and edges to draw. 
So, to generate a small part of graph you should use `generate_subgraph.py` script.

Common usage is

`python generate_subgraph.py -i graph.sif -o subgraph --reference [name of reference genome] --contig [name of contig] --start [start position] --end [end position]`

It's possible to set start and end by node names directly. For this add `--inputtype byNodeID` to command:

`python generate_subgraph.py -i graph.sif -o subgraph --reference [name of reference genome] --inputtype byNodeID --start [name of start node] --end [name of end node]`

These commands generate subgraph `subgraph.sif`, which is connected with START......END simple chain of nodes in the reference genome.

Additional parameters:
* ` --neighborhood ` - number of nodes, added to left and right side of refernce chain (default 20)
* ` --depth ` - maximum length of deviating paths, which will be added to the subgraph {default is the length of the reference chain)
* ` --tails ` - if deviating path too long, it will be replaced by left and right "tails". This parameter is tails length (default 5)
* ` --names_list ` - path to file with list of names for subgraph generating (default all names from *.sif will be used)

### Graph drawing

Now we can run drawing of our mini-graph. Just type in terminal:

`python plot_subgraph.py -i subgraph.sif -o subgrah_img`

This script generate:
* subgraph_img.ps file as image and 
* subgraph_img.dot file with DOT description of subgraph (DOT is popular graph structure description language)

Additional parameters:
* ` --freq_min` - minimal edge frequency to draw. Edge frequency is number of genomes with this edge.
* ` --da` - legacy parameter, it's not recommended to use. Draws all subgraph edges in any case, but edges wih frequency < `freq_min` do not influence to subgraph layout.


### Method validation

To verify our method we implemented a random rearrangements algorithm.
You can carry out the validation pipeline by one command: just open `geneGraph/source/validation` folder and run `bash validate.sh` script in terminal. This short pipeline generates sinus-distribution based random rearrangements models and creates the plot for the input distribution and the resulting complexity profile.


## Links

This tool is available as a web-service [Genome Complexity Browser](http://gcb.rcpcm.org) ([link to github](https://github.com/DNKonanov/Genome-Complexity-Browser)) and as a stand-alone app [GCB package](https://github.com/DNKonanov/GCB).

Genome Complexity Browser contains pre-computed complexity profiles and the graph structure for more than 140 prokariotic species.

GCB_package does not contain pre-computed complexity profiles, but there are easy-to-use scripts to add your own organisms in the package.

## References

Manolov, A., Konanov, D., Fedorov, D., Osmolovsky, I., Vereshchagin, R., & Ilina, E. (2020). [Genome Complexity Browser: Visualization and quantification of genome variability](https://doi.org/10.1371/journal.pcbi.1008222). PLoS computational biology, 16(10), e1008222. 

