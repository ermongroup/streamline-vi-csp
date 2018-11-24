Streamlined Survey Propogation
============================================

This repository provides a reference implementation for solving constraint satisfaction problems via streamlined survey propogation as described in the paper:


> Streamlining Variational Inference for Constraint Satisfaction Problems  
Aditya Grover, Tudor Achim, Stefano Ermon  
Advances in Neural Information Processing Systems (NIPS), 2018  
Paper: arXiv link TODO

## General

The codebase has been built on top of the survey propogation implementation of A. Braunstein, M. Mezard, and R. Zecchina as described in the paper "Survey propagation: an algorithm for satisfiability". It is implemented in C/C++ and tested on Ubuntu 16.04.

## Setup

To compile the binaries run the following command from the root directory

```
make all
```

This will create a binary file for `sp` in the root directory (and others which will be directly accessed by `sp`).

## Options

For a full list of options, run:

```
./sp -h
```

Key options are described below:

```
  -l CSP in CNF representation (if none provided, random k-SAT instance is generated)
  -k length of each clause 
  -n number of variables 
  -m number of clauses 
  -a clause/variable ratio
  -s seed for reproducibility
  -% percentage of paired disjunctions (denoted as R in the paper)
  -t number of streamlining iterations (denoted as T in the paper)
  -d limit on the streamlined disjunctions per variable
  -p prefix path where all the generated files (cnf formula, streamlined formula etc.) are dumped
  
```

## Examples

Baseline *survey inspired decimation* on a random 3-SAT instance with 50,000 variables and clause to variable ratio of 4.235:

```
./sp -n50000 -a4.235 -k3 -%1 -t0 -d2 -s1
```

*Survey inspired streamlining* for the same problem instance:

```
./sp -n50000 -a4.235 -k3 -%1 -t90 -d2 -s1
```

*Survey inspired streamlining* for an arbitrary CSP accessed via the filepath `csp/1.cnf`:

```
./sp -%1 -lcsp/1.cnf -t80
```


## Citing

If you find this codebase useful in your research, please consider citing the following paper:


>@inproceedings{grover2018streamlining,  
  title={Streamlining Variational Inference for Constraint Satisfaction Problems},  
  author={Grover, Aditya and Achim, Tudor and Ermon, Stefano},  
  booktitle={Advances in Neural Information Processing Systems},  
  year={2018}}