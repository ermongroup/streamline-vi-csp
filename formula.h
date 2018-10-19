/*
    Copyright 2002 Alfredo Braunstein, Michele Leone, Marc Mézard, 
                   Martin Weigt and Riccardo Zecchina

    This file is part of SP (Survey Propagation).

    SP is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or 
    (at your option) any later version.

    SP is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SP; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <stdio.h>

#define STREAMLINEFORMULA "streamlineformula.tmp.cnf"
#define FORMULA "formula.tmp.cnf"
#define SUBFORMULA "subformula.tmp.cnf"
#define NOCONVERGENCE "noconvergence.tmp.cnf"
#define SPSOL "spsol.tmp.lst"
#define SOLUTION "solution.tmp.lst"

#define WSATSOL "wsatsol.tmp.lst"
#define WSATOUT "wsat.tmp.out"
#define WSATCUTOFF 100000000

//buffer for executing extern shell commands
#define BUFFERSIZE 1000
struct pistruct {
	double p; //product of (1+eta) of literals of same sign
	double m; //product of (1+eta) of literals of opposite sign
#ifndef FAST_ITERATION
	int pzero;
	int mzero;
#endif
};

struct weightstruct {
        double p; //plus
        double z; //zero
        double m; //minus
};

struct literalstruct {
	int var;   //varnum=1,...,N
	unsigned char bar;  //bar=0,1
	double eta; //eta of the literal Q(u)=delta(u)+eta*delta(u+bar?-1:1)
};

struct clausestruct {
	struct literalstruct *literal; //list of literals
	int type;     //type=0,1,... actual number of lits
	int lits;     //lits=0,1,... initial number of lits
};

struct clauselist {
	struct clausestruct *clause; //in which clause
	int lit; //in which literal on such clause
};

struct vstruct {
	int clauses; //in how many clauses the var appears
	struct clauselist *clauselist; //list of clauses
	struct pistruct pi; //product of (1+eta) of these clauses
	int spin; //spin of the var, 0=unfixed
        int artificialclauses; //in how many artificial implied clauses the var appears
    int proxy_spin;
};

void addClause(int num_vars, int* vars, int* polarities);
bool check_unsat();
void initformula();
int readformula(FILE * source);
int readvarformula(char *);
void randomformula(int);
int writeformula(FILE *sink);
int writeformulariccardo(FILE *sink);
int writesolution(FILE *sink);
int writepartialformula(FILE *sink);
int simplify(int var);
int fix(int var, int spin);
int fixTwo(int var, int var2, int spin, int spin2);
int fixThree(int var, int var2, int var3, int spin, int spin2, int spin3);
int callwsat(int cutoff, char *formula, char *output);
int solvesubsystem();
void print_stats();
struct weightstruct normalize(struct weightstruct H);
