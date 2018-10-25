from subprocess import STDOUT, check_output, Popen, PIPE, TimeoutExpired
import argparse
from tqdm import tqdm
import collections
import math
import signal
import tempfile
import shutil
import os
from contextlib import contextmanager
import numpy as np

@contextmanager
def cwd(path):
    oldpwd=os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

def create_env():
    d = tempfile.mkdtemp()
    for executable in ['sp', 'walksat', 'merge', 'verify']:
        assert os.path.exists(executable), executable
        shutil.copy(executable, d)
    return d

def run_sat(path, timeout, streamlining_rounds=1):
    with Popen('./sp -%%1 -t%d -l%s' % (streamlining_rounds, path), shell=True, stdout=PIPE, stderr=PIPE, preexec_fn=os.setsid) as process:
        try:
            sat_output = process.communicate(timeout=timeout)[0]
        except TimeoutExpired:
            os.killpg(process.pid, signal.SIGINT) # send signal to the process group
            os.killpg(process.pid, signal.SIGKILL) # send signal to the process group
            sat_output = process.communicate()[0]

    sat_output = str(sat_output)
    is_contradiction = 'contradiction' in sat_output
    is_sat = 'ASSIGNMENT FOUND' in sat_output
    return is_sat, is_contradiction

def run_xor_trial(streamlining_rounds, num_vars, density, xor_num_vars=2, timeout=60):
    d = create_env()
    output = check_output('''python2 graph_gen/generate_xor_cnf.py \
        --num_variables=%d \
        --xor_density=%0.2f \
        --xor_num_vars=%d \
        --cnf_file_path=%s/xor.cnf \
        ''' % (num_vars, density, xor_num_vars, d), shell=True)

    with cwd(d):
        sat_output = run_sat('xor.cnf', timeout, streamlining_rounds=streamlining_rounds)
        shutil.rmtree(d)
        return sat_output

def run_k_color_trial(streamlining_rounds, num_nodes, edge_density, num_colors, timeout=60):
    d = create_env()
    output = check_output('''python2 graph_gen/generate_k_coloring_cnf.py \
        --num_nodes=%d \
        --edge_density=%0.2f \
        --num_colors=%d \
        --cnf_file_path=%s/graph.cnf \
        ''' % (num_nodes, edge_density, num_colors, d), shell=True)

    with cwd(d):
        sat_output = run_sat('graph.cnf', timeout, streamlining_rounds=streamlining_rounds)
        shutil.rmtree(d)
        return sat_output

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('num_trials', type=int, default=100)
    args = ap.parse_args()

    for density in np.linspace(2.5, 3.2, 10):
        print('density: %0.2f' % density)
        results_dict = collections.defaultdict(list)
        for _ in tqdm(range(args.num_trials)):
            for rounds in [0, 10, 50, 100]:
                is_sat, is_contradiction = run_k_color_trial(streamlining_rounds=rounds, num_nodes=800, edge_density=density, num_colors=5, timeout=30)
                results_dict[rounds].append(is_sat)

        for rounds, results in sorted(results_dict.items()):
            print(rounds, 'constraints success rate:', np.mean(results))

if __name__ == '__main__':
    main()
