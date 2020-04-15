import argparse
from gene_graph_lib.compute_complexity import GenomeGraph
import os

print('geneGraph - tool to estimate genome variability\n')

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', default='no', type=str, help='input_file (Orthofinder file)', required=True)
parser.add_argument('-o', '--output_dir', default='STDOUT', type=str, help='Output directory')
parser.add_argument('--reference', type=str, default='auto', help='name of reference genome', required=True)
parser.add_argument('--window', type=int, default=20, help='Size of window (default is 20)')
parser.add_argument('--iterations', type=int, default=500, help='number of iterations in stat computing (default is 500)')
parser.add_argument('--genomes_list', type=str, default='all', help='genomes list txt file')
parser.add_argument('--min_depth', type=int, default=0, help='min length of deviating path (default is 0)')
parser.add_argument('--max_depth', type=int, default=-1, help='max length of deviating path (default is inf)')
parser.add_argument('--fill_db', type=str, default=None, help='db path')

args = parser.parse_args()


try:
    os.stat(args.output_dir)

except FileNotFoundError:
    os.mkdir(args.output_dir)

f_params = open(args.output_dir + '/params.txt', 'w')
f_params.write('reference: ' + args.reference + '\n')
f_params.write('window: ' + str(args.window) + '\n')
f_params.write('iterations: ' + str(args.iterations) + '\n')

graph = GenomeGraph()

graph.read_graph(args.input_file, names_list=args.genomes_list)
graph.compute_complexity(args.output_dir, args.reference, window=args.window, iterations=args.iterations, min_depth=args.min_depth, max_depth=args.max_depth, save_db=args.save_db)

