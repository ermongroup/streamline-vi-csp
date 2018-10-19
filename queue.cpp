#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "queue.h"

int *queue;
int *inqueue;
int Qlength=0;
int Qsize;
int Qhead;
int Qtail;


void queue_reset()
{
	memset(inqueue,0,Qsize*sizeof(int));
	Qhead=Qtail=Qlength=0;
}

void queue_init(int size)
{
	queue = new int[size+1]; //calloc(size+1, sizeof(int));
	inqueue = new int[size+1]; //calloc(size+1, sizeof(int));
	Qsize=size+1;
	Qlength=Qhead=Qtail=0;
}

void queue_add(int i)
{
	if(inqueue[i])
		return;
	inqueue[i]=1;
	queue[Qtail++]=i;
	Qlength++;
	if(Qtail==Qhead) {
		fprintf(stderr, "queue overflow... you have found a bug\n");
		exit(-1);
	}
	if(Qtail==Qsize)
		Qtail=0;
}

int queue_get()
{
	int ret;
	if(Qhead==Qtail)
		return -1;
	Qlength--;
	ret=queue[Qhead];
	inqueue[ret]=0;
	if(++Qhead==Qsize)
		Qhead=0;
	return ret;;
}
