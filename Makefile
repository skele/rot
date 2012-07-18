# other D's: -DLINUX -DCKQ -DDUMPNEIGH -DDOCLOCK -O0 -g

GPP= g++
CFLAGS= -Wall -DSORTQ -pipe
LDLIBS= -lm
LINTFLAGS=-DSORTQ -DDOCLOCK
OBJS= rot.o distance.o misc.o mkshell.o procopt.o readdbl.o spawnit.o usage.o linuxlib.o
TARS= Makefile README marc.1 run rot.h rot.c misc.c mkshell.c procopt.c readdbl.c spawnit.c usage.c linuxlib.c distance.cu
NVCC= nvcc
NVCCFLAGS= -lgomp
NVCCSOURCE= distance.cu

all: rot

distance.o: $(NVCCSOURCE)
	$(NVCC) $(NVCCFLAGS) -c $(NVCCSOURCE)
	
rot: $(OBJS)
	$(NVCC) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS) $(NVCCFLAGS)

clean:
	rm -f rot *.o *.core *.BAK

tar: $(TARS)
	tar -s /^/rot-0.1a\\//p -czvf rot-0.1a.tgz $(TARS)

