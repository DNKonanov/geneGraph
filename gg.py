import argparse
from collections import OrderedDict
import sqlite3
from gene_graph_lib.reverse import reverse
from gene_graph_lib.find_context import find_context
import numpy as np
from gene_graph_lib.compute_complexity import GenomeGraph


print('\ngeneGraph - tool to estimate genome variability\n')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', default='no', type=str, help='input_file (Orthofinder file)', required=True)
parser.add_argument('-o', '--out_dir', type=str, help='Output directory', required=True)
parser.add_argument('--reference', type=str, default='all', help='name of the reference genome. Default all')
parser.add_argument('--window', type=int, default=20, help='Size of window (default is 20)')
parser.add_argument('--iterations', type=int, default=500, help='number of iterations in stat computing (default is 500)')
parser.add_argument('--genomes_list', type=str, default='all', help='genomes list txt file')
parser.add_argument('--min_depth', type=int, default=0, help='min length of deviating path (default is 0)')
parser.add_argument('--max_depth', type=int, default=-1, help='max length of deviating path (default is inf)')

args = parser.parse_args()



graph = OrderedDict()
length_table = OrderedDict()
coord_table = OrderedDict()


import os
try:
    os.stat(args.out_dir)

except FileNotFoundError:
    os.mkdir(args.out_dir)

print('Parsing... It may take a few minutes...')


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

out_context = open('{name}/{name}_context.txt'.format(name=args.out_dir), 'w')
og_context = find_context(graph)
for og in og_context:
	out_context.write(og + '\t' + str(og_context[og]) + '\n')


out = open('{name}/{name}.sif'.format(name=args.out_dir), 'w')


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


print('Reversing...')
graph, reversed_chains = reverse(graph, length_table)

for stamm in graph:
	for contig in graph[stamm]:
		if contig in reversed_chains[stamm]:
			coord_list[stamm][contig].reverse()


print('Database filling...')
db_ = sqlite3.connect('{name}/{name}'.format(name=args.out_dir) + '.db')

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

out_coord = open('{name}/{name}_genes.sif'.format(name=args.out_dir), 'w')
out_coord.write('genome\tcontig\tgene\tstart\tend\tdescription\n')

genome_key = 0
contig_key = 0
node_key = 0
for name in graph:
	c.execute('insert into genomes_table values(' + str(genome_key) + ', "' + name + '", "none")')
	
	for contig in graph[name]:
		if len(graph[name][contig]) < 2:
			continue
			
		c.execute('insert into contigs_table values(' + str(contig_key) + ', "' + contig + '", ' + str(genome_key) +')')

		for i in range(len(graph[name][contig])):
			gene = graph[name][contig][i]
			c.execute('insert into nodes_table values(' + str(node_key) + ',"' + gene + '", ' +str(contig_key) + ', "' + coord_table[name][gene][2] + '", ' + str(coord_list[name][contig][i][0]) + ', ' + str(coord_list[name][contig][i][1]) + ')')
			out_coord.write(name + '\t' +contig + '\t' + gene + '\t' + str(coord_list[name][contig][i][0]) + '\t' + str(coord_list[name][contig][i][1]) + '\t' + coord_table[name][gene][2] + '\n')

			if i == len(graph[name][contig]) - 1:
				node_key += 1
				continue
			
			first_coord = int(np.mean(coord_list[name][contig][i]))
			second_coord = int(np.mean(coord_list[name][contig][i+1]))


			line = '{first_gene} {second_gene} {genome} {contig} {first_coord} {second_coord}\n'.format(
				first_gene=graph[name][contig][i],
				second_gene=graph[name][contig][i+1],
				genome=name,
				contig=contig, 
				first_coord=first_coord,
				second_coord=second_coord
			)

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

print()




f_params = open(args.out_dir + '/params.txt', 'w')


for arg in args.__dict__:

	f_params.write('{}: {}\n'.format(arg, args.__dict__[arg]))
f_params.close()



full_graph = OrderedDict()
length_table = OrderedDict()
coord_table = OrderedDict()

print('    ---    \nParsing with paralogues... It may take a few minutes...')


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

		if name not in full_graph:
			full_graph[name] = {contig_name:[(start_coord, end_coord, OG)]}
			coord_table[name] = {}
			length_table[name] = {}

		elif contig_name not in full_graph[name]:
			full_graph[name][contig_name] = [(start_coord, end_coord, OG)]

		else:
			full_graph[name][contig_name].append((start_coord, end_coord, OG))
		

		coord_table[name].update([(OG, (start_coord, end_coord, description))])
		length_table[name].update([(OG, end_coord - start_coord)])

graph = OrderedDict()
for name in full_graph:
	graph[name] = full_graph[name].copy()
	for contig in graph[name]:
		full_graph[name][contig].sort()
		graph[name][contig].sort()
		graph[name][contig] = [graph[name][contig][i][2] for i in range(len(graph[name][contig]))]




out = open('{name}/{name}_pars.sif'.format(name=args.out_dir), 'w')


edges = []



G = set([])
for name in graph:

	S = set([])
	for contig in graph[name]:
		for gene in graph[name][contig]:
			if gene in S:
				G.add(gene)
				continue

			S.add(gene)


print('Orthologization...')

paralogues_table = {pg:[] for pg in G}
coord_list = {}
for name in graph:
	coord_list[name] = {}
	for contig in graph[name]:
		
		
		if len(graph[name][contig]) < 2:
			continue

		coord_list[name][contig] = []
		fixed_contig = []
		gene = 1
		graph[name][contig] = ['>'] + graph[name][contig] + ['>']
		while gene <= len(graph[name][contig]) - 2:
			context = False

			if (graph[name][contig][gene]) in paralogues_table:

				length = 1
				try:
					
					if gene - 2 < 0:
						f = 1/0

					context_left = graph[name][contig][gene - 1]
					start_par_coord = full_graph[name][contig][gene - 2][1]

				except:
					context_left = '>'
					start_par_coord = full_graph[name][contig][gene-1][0]
				
				current_gene = graph[name][contig][gene]
				while graph[name][contig][gene + 1] == current_gene:
					gene += 1
					length += 1

				try:
					context_right = graph[name][contig][gene + 1]
					end_par_coord = full_graph[name][contig][gene][0]
				
				except:
					context_right = '>'
					end_par_coord = full_graph[name][contig][gene-1][1]
					

				for i in paralogues_table[graph[name][contig][gene]]:
					if context_left in i and context_right in i:

						new_name = 'P' + graph[name][contig][gene][1:] + '_' + str(i[2])

						if new_name in fixed_contig:
							copy = 0
							while new_name + '_copy' + str(copy) in fixed_contig:
								copy += 1

							new_name = new_name + '_copy' + str(copy)

						fixed_contig.append(new_name)
						coord_list[name][contig].append((start_par_coord, end_par_coord))
						length_table[name].update([(new_name, length_table[name][graph[name][contig][gene]])])
						
						coord_table[name].update([(new_name, (start_par_coord, end_par_coord, coord_table[name][graph[name][contig][gene]][2]))])

						context = True
						break
				if context != True:
					par = len(paralogues_table[graph[name][contig][gene]])
					paralogues_table[graph[name][contig][gene]].append([context_left, context_right, par])

					
					new_name  = 'P' + graph[name][contig][gene][1:] + '_' + str(paralogues_table[graph[name][contig][gene]][-1][2])

					fixed_contig.append(new_name)
					coord_list[name][contig].append((start_par_coord, end_par_coord))

					length_table[name].update([(new_name, length_table[name][graph[name][contig][gene]])])
					coord_table[name].update([(new_name, (start_par_coord, end_par_coord, coord_table[name][graph[name][contig][gene]][2]))])
					

				gene += 1
				continue
			


			fixed_contig.append(graph[name][contig][gene])
			coord_list[name][contig].append((coord_table[name][graph[name][contig][gene]][0], coord_table[name][graph[name][contig][gene]][1]))
			gene += 1

		graph[name][contig] = fixed_contig


print('Reversing...')
graph, reversed_chains = reverse(graph, length_table)




for name in graph:
	for contig in graph[name]:
		if contig in reversed_chains[name]:
			coord_list[name][contig].reverse()



print('Database filling...')

db_ = sqlite3.connect('{name}/{name}_pars.db'.format(name=args.out_dir))
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

out_coord = open('{name}/{name}_pars_genes.sif'.format(name=args.out_dir), 'w')
out_coord.write('genome\tcontig\tgene\tstart\tend\tdescription\n')

genome_key = 0
contig_key = 0
node_key = 0
for name in graph:

	c.execute('insert into genomes_table values(' + str(genome_key) + ', "' + name + '", "none")')
	
	for contig in graph[name]:
		if len(graph[name][contig]) < 2:
			continue
			
		c.execute('insert into contigs_table values(' + str(contig_key) + ', "' + contig + '", ' + str(genome_key) +')')

		for i in range(len(graph[name][contig])):
			gene = graph[name][contig][i]
			c.execute('insert into nodes_table values(' + str(node_key) + ',"' + gene + '", ' +str(contig_key) + ', "' + coord_table[name][gene][2] + '", ' + str(coord_list[name][contig][i][0]) + ', ' + str(coord_list[name][contig][i][1]) + ')')
			out_coord.write(name + '\t' +contig + '\t' + gene + '\t' + str(coord_list[name][contig][i][0]) + '\t' + str(coord_list[name][contig][i][1]) + '\t' + coord_table[name][gene][2] + '\n')

			if i == len(graph[name][contig]) - 1:
				node_key += 1
				continue


			first_coord = int(np.mean(coord_list[name][contig][i]))
			second_coord = int(np.mean(coord_list[name][contig][i+1]))


			line = '{first_gene} {second_gene} {genome} {contig} {first_coord} {second_coord}\n'.format(
				first_gene=graph[name][contig][i],
				second_gene=graph[name][contig][i+1],
				genome=name,
				contig=contig, 
				first_coord=first_coord,
				second_coord=second_coord
			)
			
			
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

print('Complete!')



from pickle import dump
from gene_graph_lib.compute_complexity import GenomeGraph
import os
import sqlite3


f = open('scripts/strains_decode.txt')

codes = {}

for line in f:
    code = line.split(' ')[-1][:-3]
    name = ''.join(line.split(' ')[:-1])
    codes[code] = name

org = args.out_dir

files = os.listdir(org + '/')
if '.dump' not in ' '.join(files):
	print(org, end=' ')
	print('dumping...')
	
	
	try:
		g = GenomeGraph()
		g.read_graph(org + '/' + org + '.sif')
		dump_file = open(org + '/' + org + '.dump', 'wb')
		dump(g, dump_file)
	except:
		pass

	try:
		g = GenomeGraph()
		g.read_graph(org + '/' + org + '_pars.sif')
		dump_file = open(org + '/' + org + '_pars.dump', 'wb')
		dump(g, dump_file)
	except:
		pass

	skip = os.path.isfile(org + '/' + org + '.db')

	if skip == True:
		connect = sqlite3.connect(org + '/' + org + '.db')
		c = connect.cursor()
		
		genome_codes = [q for q in c.execute('select genome_id,genome_code from genomes_table')]

		for code in genome_codes:
			for ref_code in codes:
				if ref_code in code[1]:
					c.execute('update genomes_table set genome_name="' + codes[ref_code] + '" where genome_id=' + str(code[0]))

		
		connect.commit()
		connect.close()

	#pars table
	skip = os.path.isfile( org + '/' + org + '_pars.db')

	if skip == True:
			
		connect = sqlite3.connect(org + '/' + org + '_pars.db')
		c = connect.cursor()
		
		genome_codes = [q for q in c.execute('select genome_id,genome_code from genomes_table')]

		for code in genome_codes:
			for ref_code in codes:
				if ref_code in code[1]:
					c.execute('update genomes_table set genome_name="' + codes[ref_code] + '" where genome_id=' + str(code[0]))

		connect.commit()
		connect.close()
	



references_list = list(graph.keys())
if args.reference == 'all':
	pass

else:
	if args.reference not in references_list:
		print('{} is not found in the selected sif file!'.format(args.reference))
		import sys
		sys.exit()
		
	references_list = [args.reference]

graph = GenomeGraph()


imput_file = '{name}/{name}.sif'.format(name=args.out_dir)




graph.read_graph('{name}/{name}.sif'.format(name=args.out_dir), names_list=args.genomes_list)
for reference in references_list:

    outdir = args.out_dir + '/{}'.format(reference.replace('/', '_'))

    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass


    graph.compute_complexity(
        args.out_dir + '/{}'.format(reference.replace('/', '_')), 
        reference, 
        window=args.window, 
        iterations=args.iterations, 
        min_depth=args.min_depth, 
        max_depth=args.max_depth, 
        save_db='{name}/{name}.db'.format(name=args.out_dir)
        )




