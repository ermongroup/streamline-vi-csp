CURDIR="$(shell pwd)"

FLAGS=-D__DIR__='${CURDIR}' -Wall -Winline -O3 -Wno-unused-result -std=c++11
#FLAGS=-D__DIR__='${CURDIR}' -Wall -Winline -g -O0 -Wno-unused-result
#FLAGS=-Wall -Winline -O3 -funroll-loops
#FLAGS=-Wall -Winline -O3 -funroll-loops
#FLAGS=-g -Wall -pg -O0 
#FLAGS=-g
CC=g++

LDFLAGS=-lboost_iostreams -lm -lboost_system -lboost_filesystem

all:	sp merge verify walksat
sp: sp.cpp sp.h formula.h formula.o random.o queue.o
	${CC} ${FLAGS} sp.cpp random.o formula.o queue.o -o sp ${LDFLAGS}
formula.o: formula.h formula.cpp
	${CC} ${FLAGS} -c formula.cpp -o formula.o ${LDFLAGS}
random.o: random.c random.h
	${CC} ${FLAGS} -c random.c -o random.o
queue.o: queue.cpp queue.h
	${CC} ${FLAGS} -c queue.cpp -o queue.o
clean:
	rm -f core *~ *.o  gmon.out *.tmp.* *.out *.dat sp walksat verify merge
binclean: clean
	rm -f walksat sp spprof verify merge
walksat: walksat.c
	gcc -w ${LDFLAGS} -O6 walksat.c -o walksat -lm
merge: merge.cpp
	${CC} -Wall merge.cpp -o merge
verify: verify.cpp
	${CC} -Wall verify.cpp -o verify
