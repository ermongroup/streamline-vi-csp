#!/usr/bin/env python

import subprocess as sp
import argparse as ap
import tempfile
import numpy as np
import sys
import random
import os
import shutil

parser = ap.ArgumentParser()
parser.add_argument('--num_variables', type=int, default=5)
parser.add_argument('--xor_density', type=float, default=1)
parser.add_argument('--xor_num_vars', type=int, default=3)
parser.add_argument('--cnf_file_path', type=str, default='xor.cnf')
args = parser.parse_args()

assert args.cnf_file_path.endswith('.cnf'), ('Filename must end in .cnf', args.cnf_file_path)
assert not os.path.exists(args.cnf_file_path), ('File already exists:', args.cnf_file_path)

from sympy.logic.boolalg import ITE, And, Xor, Or, Not
import sympy
from sympy.logic.boolalg import to_cnf

num_xors = int(args.num_variables * args.xor_density)

variables = sympy.symbols(','.join(map(str, range(1, args.num_variables+1))))
xors = []
for _ in range(num_xors):
    xor_vars = random.sample(variables, args.xor_num_vars)
    for i in range(len(xor_vars)):
        if np.random.random() < 0.5:
            xor_vars[i] = Not(xor_vars[i])
    xor = Xor(*xor_vars)
    if np.random.random() < 0.5:
        xor = Not(xor)
    xors.append(xor)
formula = And(*xors)
cnf_formula = to_cnf(formula)

with open(args.cnf_file_path, 'w') as f:
    print >>f, 'p cnf %d %d' % (args.num_variables, len(cnf_formula.args))

    for clause in cnf_formula.args:
        for variable in clause.args:
            if isinstance(variable, sympy.Symbol):
                print >>f, variable.name,
            else:
                print >>f, '-%s' % (variable.args[0].name),
        print >>f, '0'

print 'Finished writing to', args.cnf_file_path
