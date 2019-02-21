import argparse
from collections import OrderedDict
import sqlite3
import reverse
from find_context import find_context

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', default='no', type=str, help='input_file, generated by Orthofinder')
parser.add_argument('-o', '--out_file', default='paths', type=str, help='output file prefix (default paths)')
args = parser.parse_args()


graph = OrderedDict()
length_table = OrderedDict()
coord_table = OrderedDict()


for line in open(args.input_file, 'r'):

	OG, string = line.split(': ')[0], line.split(': ')[1][:-1]


	stamms = string.split(' ')
	for stamm in stamms:
		name = stamm.split('|')[0]
		start_coord = int(stamm.split('|')[-2])
		end_coord = int(stamm.split('|')[-1])
		coord = int(stamm.split('|')[1])
		description = stamm.split('|')[2]
		contig_name = stamm.split('|')[3]

		if name not in graph:
			graph[name] = {contig_name:[(start_coord, OG)]}
			coord_table[name] = {}
			length_table[name] = {}

		elif contig_name not in graph[name]:
			graph[name][contig_name] = [(start_coord, OG)]

		else:
			graph[name][contig_name].append((start_coord,OG))
		
		coord_table[name].update([(OG, [start_coord, end_coord, description])])
		length_table[name].update([(OG, end_coord - start_coord)])


for name in graph:
	 for contig in graph[name]:
		 graph[name][contig].sort()
		 graph[name][contig] = [graph[name][contig][i][1] for i in range(len(graph[name][contig]))]

out_context = open(args.out_file + '_context.txt', 'w')
og_context = find_context(graph)
for og in og_context:
	out_context.write(og + '\t' + str(og_context[og]) + '\n')


out = open(args.out_file + '.sif', 'w')


coord_list = {}
for name in graph:
	G = set([])
	S = set([])
	for contig in graph[name]:
		for gene in graph[name][contig]:
			if gene in S:
				G.add(gene)
				continue

			S.add(gene)
	coord_list[name] = {}
	for contig in graph[name]:
		coord_list[name][contig] = []
		fixed_contig = []
		for gene in graph[name][contig]:
			if gene in G:
				continue
			fixed_contig.append(gene)
			coord_list[name][contig].append([coord_table[name][gene][0], coord_table[name][gene][1]])

		graph[name][contig] = fixed_contig.copy()

graph, reversed_chains = reverse.reverse(graph, length_table)

for stamm in graph:
	for contig in graph[stamm]:
		if contig in reversed_chains[stamm]:
			print(contig)
			coord_list[stamm][contig].reverse()



db_ = sqlite3.connect(args.out_file + '.db')
c = db_.cursor()

c.execute('drop table if exists genomes_table')
c.execute('drop table if exists contigs_table')
c.execute('drop table if exists nodekey')
c.execute('drop table if exists nodes_table')
c.execute('drop table if exists complexity_table')
c.execute('drop table if exists edges_table')

c.execute('''
create table if not exists genomes_table(
	genome_id integer primary key,
	genome_code text,
	genome_name text

);''')

c.execute('''
create table if not exists contigs_table(
	contig_id integer primary key,
	contig_code text,
	genome_id integer not null,
	foreign key (genome_id) references genomes_table (genome_id)
);''')

c.execute('''
create table if not exists nodekey (
	node_id integer primary key,
	node_name text
);''')
c.execute('''
create table if not exists nodes_table (
	node_id integer primary key,
	node_name text,
	contig_id integer not null, 
	description text,
	start_coord integer,
	end_coord integer,
	foreign key (contig_id) references contigs_table (contig_id)
);''')
c.execute('''
create table if not exists complexity_table (
	node_id integer not null,
	contig_id integer not null,
	window_complexity float,
	prob_window_complexity float,
	io_complexity integer,
	prob_io_complexity integer,
	window integer,
	foreign key (node_id) references nodes_table (node_id),
	foreign key (contig_id) references contigs_table (contig_id)
);''')
c.execute('''
create table if not exists edges_table (
	source integer not null,
	target integer not null,
	frequency integer,
	genomes text,
	foreign key (source) references nodekey (node_id),
	foreign key (target) references nodekey (node_id)
);''')




nodes_set = set([])
for genome in graph:
	for contig in graph[genome]:
		chain = graph[genome][contig]

		nodes_set.update(set(chain))

nodes_codes_query = []
nodes_codes = {}
node_id = 0
for node in nodes_set:
	c.execute('insert into nodekey values (' + str(node_id) + ', "' + node + '")')
	nodes_codes[node] = node_id
	node_id += 1

db_.commit()

out_coord = open(args.out_file + '_genes.sif', 'w')
out_coord.write('genome\tcontig\tgene\tstart\tend\tdescription\n')

genome_key = 0
contig_key = 0
node_key = 0
for name in graph:
	print(name)
	print('---')

	c.execute('insert into genomes_table values(' + str(genome_key) + ', "' + name + '", "none")')
	
	for contig in graph[name]:
		if len(graph[name][contig]) < 2:
			continue
			
		c.execute('insert into contigs_table values(' + str(contig_key) + ', "' + contig + '", ' + str(genome_key) +')')

		print(contig)

		for i in range(len(graph[name][contig])):
			gene = graph[name][contig][i]
			c.execute('insert into nodes_table values(' + str(node_key) + ',"' + gene + '", ' +str(contig_key) + ', "' + coord_table[name][gene][2] + '", ' + str(coord_list[name][contig][i][0]) + ', ' + str(coord_list[name][contig][i][1]) + ')')
			out_coord.write(name + '\t' +contig + '\t' + gene + '\t' + str(coord_list[name][contig][i][0]) + '\t' + str(coord_list[name][contig][i][1]) + '\t' + coord_table[name][gene][2] + '\n')

			if i == len(graph[name][contig]) - 1:
				node_key += 1
				continue
			line = graph[name][contig][i] + ' ' + graph[name][contig][i + 1] + ' ' + name + ' ' + contig + '\n'
			out.write(line)
			node_key += 1
		contig_key += 1		
	genome_key += 1

out_coord.close()
out.close()
freq = {}
for genome in graph:
	for contig in graph[genome]:
		chain = graph[genome][contig]

		for i in range(len(chain) - 1):
			if (chain[i], chain[i+1]) not in freq:
				freq[(chain[i], chain[i+1])] = [genome]

			else:
				freq[(chain[i], chain[i+1])].append(genome)

edges_query = []


for edge in freq:
	edges_query.append({'source': nodes_codes[edge[0]], 'target': nodes_codes[edge[1]], 'frequency': len(freq[edge]), 'genomes': '\n'.join(freq[edge])})
	c.execute('insert into edges_table values(' + str(nodes_codes[edge[0]]) + ', ' + str(nodes_codes[edge[1]]) + ', ' + str(len(freq[edge])) + ', "' + '\n'.join(freq[edge]) + '")')


db_.commit()
db_.close()
