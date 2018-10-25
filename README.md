# geneGraph

Tool for genome complexity computing

##Dependencies

* Python 3.5 or later
* gene_graph_lib

    You can install gene_graph_lib with <pip install gene_graph_lib>

##Usage

If you have not sif file with graph stucture, the first step is parsing of Orthofinder outputs.
To do this type in terminal:
<python orthofinder_parse.py -i [path to txt file] -o [path and name prefix for output files]>

Output files:
* **prefix.sif** - all edges list of genome graph
* **prefix_freq.sif** - all edges frequency
* **prefix.db** - database with all parsed informmation
* **prefix_context.sif** - number of contexts, computed for each node in graph

Next step is computing of genome complexity.
To do this type in terminal:

<python start_computing.py -i prefix.sif -o [path to output folder] --reference_stamm [name of reference strain]>

Additional parameters are:
* < --window> - window size (default 20)
* < --iterations> - number of iterations in probabilistic method(default is 500)
* < --names> - path to file with names list (default all names from *.sif will be used)
* < --min_depth, --max_depth> - minimum and maximum depth of generated paths in graph (default from 0 to inf)
* < --save_db> - path to database, created by orthfinder_parse.py (default data dont saved to db, only to txt)

Output files for each contig in reference genome:
* all_bridges_contig_n.txt
* window_variability_contig_n.txt
* main_chain_contig_n.txt
* params.txt
* IO_variability_table_contig_n.txt
* prob_IO_variability_table_contig_n.txt
* prob_window_variability_contig_n.txt
* prob_all_bridges_contig_n.txt