from gene_graph_lib.compute_complexity import GenomeGraph
from pickle import dump, load
from time import time


def find_uniq_paths(g, gene, ref, direction, depth=100):


    paths = []

    if direction == '+':
        step = 1
    elif direction == '-':
        step = -1

    for strain in g. list_graph:
        for contig in g.list_graph[strain]:
            c = g.list_graph[strain][contig]

            try:
                index = c.index(gene)
            except: continue
            
            path = [gene]
            index += step
            while len(path) <= depth:
                
                try:
                    test = c[index]
                except:
                    break
                if c[index] in ref:
                    path.append(c[index])
                    break

                path.append(c[index])
                index += step

            if path[-1] in ref:
                continue
            paths.append(path)


    used_genes = set([])
    for p in paths:
        for gene in p[:-1]:
            used_genes.add(gene)
        
    cleared_paths = []
    for p in paths:
        if p[-1] not in used_genes:
            if tuple(p) not in cleared_paths:
                cleared_paths.append(tuple(p))

    
    return len(cleared_paths)            


def check_edge_weight(g, gene_1, gene_2):
    
    value1 = value2 = 0
    for gene in g.dict_graph_freq[gene_1]:
        if gene[0] == gene_2:
            value1 = gene[1]
            break

    for gene in g.dict_graph_freq[gene_2]:
        if gene[0] == gene_1:
            value2 = gene[1]
            break

    return(max(value1, value2))

def check_insertion(g, gene_1, gene_2, insertion_len):
    
    insert_count = 0

    for strain in g.list_graph:
        for contig in g.list_graph[strain]:
            c = g.list_graph[strain][contig]

            try:
                d = abs(c.index(gene_1) - c.index(gene_2))
                if d < insertion_len:
                    insert_count += 1
            
            except:
                continue

    if insert_count > 0:
        return True

    else:
        return False

def find_template(g, depth=100, insertion_len=20, weight=0.05, uniq_paths=1):

    for name in g.list_graph:
        print(name)
        for contig in g.list_graph[name]:
            c = g.list_graph[name][contig]

            for i  in range(len(c)):
                insert = check_insertion(g, c[i], c[i + 1], insertion_len)
                if insert == False:
                    continue
                
                ref = c[max(0, i - depth):min(i + depth, len(c))]

                f_hits = find_uniq_paths(g, c[i+1], ref, '+',depth=depth)
                r_hits = find_uniq_paths(g, c[i], ref, '-', depth=depth)

                if f_hits >=uniq_paths and r_hits >= uniq_paths and check_edge_weight(g, c[i], c[i+1]) >= weight*len(g.list_graph):
                    print(name, g.genes_decode[c[i+1]], f_hits, g.genes_decode[c[i]], r_hits)


g = GenomeGraph()
g.read_graph('/home/dmitry/projects/Genome-Complexity-Browser-master/gcb_server/data/Escherichia_coli/Escherichia_coli.sif', generate_freq=True)
find_template(g, depth=20, insertion_len=50, uniq_paths=1)













