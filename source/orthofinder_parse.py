import argparse
import reverse
from collections import OrderedDict
from find_context import find_context
import manage_db

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

out_context = open(args.out_file + '_context.txt', 'a+')
og_context = find_context(graph)
for og in og_context:
	out_context.write(og + '\t' + str(og_context[og]) + '\n')


out = open(args.out_file + '.sif', 'a+')
out_freq = open(args.out_file + '_freq.sif', 'a+')


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


db = manage_db.create_connection(args.out_file + '.db')
sql_exec = """ CREATE TABLE IF NOT EXISTS stamms_table (
					id integer PRIMARY KEY,
					stamm text
				); """
manage_db.create_table(db, sql_exec)

sql_exec = """ CREATE TABLE IF NOT EXISTS contigs_table (
					id integer PRIMARY KEY,
					stamm_key text,
					contig text
				); """
manage_db.create_table(db, sql_exec)

sql_exec = """ CREATE TABLE IF NOT EXISTS og_table (
					id integer PRIMARY KEY,
					contig integer,
					og text,
					description text,
					start_coord integer,
					end_coord integer
				); """
manage_db.create_table(db, sql_exec)

sql_exec = """ CREATE TABLE IF NOT EXISTS og_complexity_table (
					contig integer,
					og text,
					win_var float,
					prob_win_var float,
					io float,
					prob_io float
				); """
manage_db.create_table(db, sql_exec)

sql_exec = """ CREATE TABLE IF NOT EXISTS freq_table (
					id integer PRIMARY KEY,
					edge_freq int,
					edge text,
					stamms_list text
				);"""

manage_db.create_table(db, sql_exec)


cursor = db.cursor()
stamm_key = 0
contig_key = 0
og_key = 0
for name in graph:
	print(name)
	print('---')
	cursor.execute('INSERT INTO stamms_table VALUES (' + str(stamm_key) + ',"' + name + '")')
	for contig in graph[name]:
		
		if len(graph[name][contig]) < 2:
			continue
			
		cursor.execute('INSERT INTO contigs_table VALUES (' + str(contig_key) + ',' + str(stamm_key) + ',"' + str(contig) + '")')
		
		for i in range(len(graph[name][contig])):
			gene = graph[name][contig][i]
			cursor.execute('INSERT INTO og_table VALUES('+ str(og_key) + ',' + str(contig_key) + ',"' + gene + '","' + coord_table[name][gene][2] + '",' + str(coord_list[name][contig][i][0]) + ',' + str(coord_list[name][contig][i][1]) + ')')

			og_key += 1
					

			if i == len(graph[name][contig]) - 1:
				continue
			line = graph[name][contig][i] + ' ' + graph[name][contig][i + 1] + ' ' + name + ' ' + contig + '\n'
			out.write(line)
		contig_key += 1	
	stamm_key += 1

freq = {}
for name in graph:
	for contig in graph[name]:
		for i in range(len(graph[name][contig]) - 1):
			if graph[name][contig][i] + ' ' + graph[name][contig][i + 1] not in freq:
				freq.update([(graph[name][contig][i] + ' ' + graph[name][contig][i + 1], [name])])

			else:
				freq[graph[name][contig][i] + ' ' + graph[name][contig][i + 1]].append(name)


freq_key = 0
for pair in freq:
	stamms_list = ''

	out_freq.write(pair + ' ' + str(len(freq[pair])) + ' ')
	for name in freq[pair]:
		out_freq.write(name)
		stamms_list += name
		if name is freq[pair][-1]:
			out_freq.write('\n')
			continue
		out_freq.write('|')
		stamms_list += ('\n')
	cursor.execute('INSERT INTO freq_table VALUES (' + str(freq_key) + ',' + str(len(freq[pair])) + ',"' + pair + '", ' + '"' + stamms_list + '")')
	freq_key += 1




db.commit()

db.close()
