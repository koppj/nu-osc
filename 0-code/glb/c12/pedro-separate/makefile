# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1


prefix = /home/jkopp/software/installed/globes-devel
exec_prefix = ${prefix}
libdir = ${exec_prefix}/lib
globesconf= $(exec_prefix)/bin/globes-config

CFLAGS = -std=c99 -Wall -g -O2 -D_BSD_SOURCE -fgnu89-inline -ffast-math $(INCFLAGS)

INCFLAGS:=$(shell $(globesconf) --include) -I${HOME}/projects/globes/Distribution/source
LDFLAGS:=$(shell $(globesconf) --libs)
LTLDFLAGS:=$(shell $(globesconf) --ltlibs)

karmen: karmen.c myprob.c degfinder.c chiNPwrapper.c
	gcc -c karmen.c  myprob.c  degfinder.c chiNPwrapper.c ${CFLAGS}
	gcc karmen.o myprob.o degfinder.o chiNPwrapper.o -o karmen  ${CFLAGS} ${INCFLAGS} $(LDFLAGS)
	./karmen

lsnd: lsnd.c myprob-lsnd.c degfinder.c chiNPwrapper.c
	gcc -c lsnd.c  myprob-lsnd.c  degfinder.c chiNPwrapper.c ${CFLAGS}
	gcc lsnd.o myprob-lsnd.o degfinder.o chiNPwrapper.o -o lsnd  ${CFLAGS} ${INCFLAGS} $(LDFLAGS)
	./lsnd

lsnd-pc: lsnd.c myprob-lsnd.c degfinder.c chiNPwrapper.c
	gcc -c lsnd.c  myprob-lsnd.c  degfinder.c chiNPwrapper.c -std=c99 -ffast-math -O3  -m64
	gcc lsnd.o myprob-lsnd.o degfinder.o chiNPwrapper.o -o lsnd  $(LDFLAGS) $(local_LDFLAGS) -std=c99 -ffast-math -O3 -m64
	./lsnd

spectrum:
	globes  KARMEN-nuecarbon.glb -s --parameters="0,0,0,0,0,0" -o nue-carbon.dat

spec:
	globes  LSND-nuecarbon.glb -s --parameters="0,0,0,0,0,0" -o LSND-nue-carbon.dat


