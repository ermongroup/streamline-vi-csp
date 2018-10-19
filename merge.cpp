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

int main(int argc, char ** argv)
{	
	int i,num,max=0;
	int *sol;
	FILE * WSS, *KPS;
	if(argc<3) {
		fprintf(stderr, "%s spsolution wsatsolution\n",argv[0]);
		exit(-1);
	} 
	// else {
	// 	fprintf(stdout, "merging %s over %s\n", argv[1],argv[2]);
	// }
	KPS=fopen(argv[1],"r");
	WSS=fopen(argv[2],"r");
	if(!KPS || !WSS) {
		printf("file not found\n");
	}
	while(fscanf(WSS,"%i",&num)==1) {
		if(abs(num)>max)
			max=abs(num);
	}
	while(fscanf(KPS,"%i",&num)==1) {
		if(abs(num)>max)
			max=abs(num);
	}
	sol=new int[max+1]; //calloc(max+1,sizeof(int));
	if(!sol) {
		fprintf(stderr, "merge: not enough memory\n");
		exit(-1);
	}
	fclose(KPS);
	fclose(WSS);
	KPS=fopen(argv[1],"r");
	WSS=fopen(argv[2],"r");
	while(fscanf(WSS,"%i",&num)==1) {
		sol[abs(num)]=num;
	}
	while(fscanf(KPS,"%i",&num)==1) {
		sol[abs(num)]=num;
	}
	for(i=1; i<=max; i++) {
		printf("%i\n",sol[i]);
	}
	
	return 0;

}
