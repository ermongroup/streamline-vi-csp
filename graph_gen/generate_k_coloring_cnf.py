#!/usr/bin/env python

import subprocess as sp
import argparse as ap
import tempfile
import numpy as np
import sys
import os
import shutil

parser = ap.ArgumentParser()
parser.add_argument('--num_nodes', type=int, default=5)
parser.add_argument('--edge_density', type=float, default=1)
parser.add_argument('--num_colors', type=int, default=3)
parser.add_argument('--cnf_file_path', type=str, default='kgraph.cnf')
args = parser.parse_args()

assert args.cnf_file_path.endswith('.cnf'), ('Filename must end in .cnf', args.cnf_file_path)
assert not os.path.exists(args.cnf_file_path), ('File already exists:', args.cnf_file_path)

print 'EDGE DENSITY SHOULD SCALE AS k ln k (https://arxiv.org/pdf/0911.2322.pdf)'

assert args.edge_density <= args.num_nodes - 1

directory = tempfile.mkdtemp()
graph_txt_path = os.path.join(directory, 'graph.txt')

num_edges = int(args.num_nodes * args.edge_density)
edges = set()
while len(edges) < num_edges:
    i = np.random.randint(1, args.num_nodes+1)
    j = np.random.randint(1, args.num_nodes+1)
    i, j = min(i, j), max(i, j)

    if i == j: continue
    if (i, j) in edges: continue
    edges.add((i, j))

with open(graph_txt_path, 'w') as f:
    print >>f, 'p %d %d %d' % (args.num_colors, args.num_nodes, num_edges)
    for i, j in sorted(edges):
        print >>f, 'e %d %d' % (i, j)

with open(args.cnf_file_path, 'w') as f:
    graph2cnf_path = os.path.join(os.path.dirname(__file__), 'graph2cnf.py')
    cmd = [sys.executable, graph2cnf_path, graph_txt_path]
    sp.call(cmd, stdout=f)

shutil.rmtree(directory)

print 'Finished writing to', args.cnf_file_path
