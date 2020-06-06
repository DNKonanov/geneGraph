from gene_graph_lib.compute_complexity import GenomeGraph
import time 
import sys
import argparse


global BIG_VALUE
BIG_VALUE = 999999999999

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', default='no', type=str, help='input_file')
parser.add_argument('-o', '--outfile', default='subgraph', help='out file prefix')
parser.add_argument('--inputtype', type=str, default='byposition', help='type of start and end parameters. Can be "byposition" or "byNodeID". Default is "byposition"')
parser.add_argument('--reference', type=str, default='auto', help='name of reference genome')
parser.add_argument('--contig', type=str, default=None, help='Contig name. Requiered, if there are more than one contigs in the reference genome')
parser.add_argument('--neighborhood', type=int, default=5, help='Number of neighborhood (default is 5)')
parser.add_argument('--start', type=str, default=None, help='start og for subgraph generating')
parser.add_argument('--end', type=str, default=None, help='end og for subgraph generating')
parser.add_argument('--depth', type=int, default=-1, help='max paths depth (default len of baseline)')
parser.add_argument('--tails', type=int, default=5, help='length of tails (default 5)')
parser.add_argument('--names_list', type=str, default='all', help='names list txt file')


args = parser.parse_args()

graph = GenomeGraph()
graph.read_graph(args.input_file)


contigs = list(graph.list_graph[args.reference].keys())
if len(contigs) > 1 and args.contig == None:
	print('There is {} contigs in the reference genome! Please, select contig with --contig parameter')
	sys.exit()

elif args.contig not in contigs and args.contig is not None:
	print('Selected contig is not in the reference genome! Please, check contig and reference names.')
	sys.exit()

elif len(contigs) == 1:
	contig = contigs[0]

else:
	contig = args.contig



def extract_nodes(graph, contig, startcoord, endcoord):
	
	global BIG_VALUE

	start_gene = None
	start_distance = BIG_VALUE
	end_gene = None
	end_distance = BIG_VALUE

	for gene in graph.genes_info[contig]:
		print(gene)
		cur_st_dist = abs(int(graph.genes_info[contig][gene]) - start_distance)
		cur_end_dist = abs(int(graph.genes_info[contig][gene]) - end_distance)

		if cur_st_dist < start_distance:
			start_gene = gene
			start_distance = cur_st_dist

		if cur_end_dist < end_distance:
			end_gene = gene
			end_distance = cur_end_dist

	return start_gene, end_gene

if args.inputtype == 'byposition':
	start, end = extract_nodes(graph, contig, int(args.start), int(args.end))
elif args.inputtype == 'byNodeID':
	start, end = args.start, args.end




subgraph, aim_chain = graph.generate_subgraph(start, end, reference=args.reference, window=args.window, tails=args.tails, depth=args.depth)
f_out = open(args.outfile + '.sif', '+a')
f_freq = open(args.outfile + '_freq.sif', '+a')

in_aim = 0
is_ref = 0
for stamm in subgraph:
	for contig in subgraph[stamm]:
		for i in range(len(contig) - 1):

			try:
				if aim_chain.index(contig[i + 1]) - aim_chain.index(contig[i]) == 1:
					in_aim = 1
			except ValueError:
				pass

			try:
				if subgraph[args.reference][0].index(contig[i + 1]) - subgraph[args.reference][0].index(contig[i]) == 1:
					is_ref = 1
			except ValueError:
				pass


			line = graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]] + ' ' + stamm + ' ' + str(in_aim) + ' ' + str(is_ref) + '\n'
			f_out.write(line)

			in_aim = 0
			is_ref = 0

f_out.close()
freq = {}

for name in subgraph:
	for contig in subgraph[name]:
		for i in range(len(contig) - 1):
			if graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]] not in freq:
				freq[graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]]] = [name]

			else:
				freq[graph.genes_decode[contig[i]] + ' ' + graph.genes_decode[contig[i + 1]]].append(name)


for pair in freq:
	f_freq.write(pair + ' ' + str(len(freq[pair])) + ' ')
	for name in freq[pair]:
		f_freq.write(name)
		if name == freq[pair][-1]:
			f_freq.write('\n')
			continue
		f_freq.write('|')

f_freq.close()