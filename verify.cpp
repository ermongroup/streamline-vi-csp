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
#include <stdlib.h>

int *solution;
int unfixed=0;

FILE * sol;
int maxvar=0;


int readsolution(char *name)
{
	int num,i;
	sol=fopen(name,"r");
	while(fscanf(sol, "%i",&num)==1) {
		if(maxvar<abs(num))
			maxvar=abs(num);
	}
	fclose(sol);
	printf("%i variables found\n",maxvar);
	solution=new int[maxvar]; //calloc(maxvar,sizeof(int));
	if(!solution) {
		fprintf(stderr, "verify: not enough memory\n");
		exit(-1);
	}
	sol=fopen(name,"r");
	while(fscanf(sol, "%i",&num)==1) {
		solution[abs(num)]=num;
	}
	for(i=1; i<=maxvar; i++) {
		if(solution[i]==0)
			unfixed++;
	}
	if(unfixed) {
		printf("%i unfixed variables\n",unfixed);
		exit (-1);
	}
	return 0;
}

int readformula()
{
	int num,is_true=0,clause=0,unsat=0,M,N,aux,sat=0;
	while((aux=getc(stdin))=='c') {
		while(getc(stdin)!='\n');
	}
	ungetc(aux,stdin);
	scanf("p cnf %i %i", &N, &M);
	while(scanf("%i",&num)==1) {
		if(num==0) {
			if(!is_true) {
				printf("clause %i not satisfied\n",clause);
				unsat++;
			} else {
				sat++;
			}
			is_true=0;
			clause++;
		} else {
			if(solution[abs(num)]==num)
				is_true=1;
			if(abs(num)<1 || abs(num)>maxvar) {
				printf("variable %i out of range\n", abs(num));
				exit(-1);
			}
		}
	}
	if(M!=clause) {
		printf("wrong number of clauses %i != %i\n", M , clause);
		return -1;
	}
	printf("%i clauses sat, %i unsat (out of %i read)\n",sat,unsat,clause);
	if(unsat)
		exit (-1);
	return 0;
}


int main(int argc, char **argv)
{
	if(argc<2) {
		fprintf(stderr, "Usage: %s solution < formula\n",argv[0]);
		return -1;
	}
	printf("verifying solution %s...", argv[1]);
	readsolution(argv[1]);
	readformula();
	return 0;
}
