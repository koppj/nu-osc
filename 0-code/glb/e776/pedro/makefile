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

#e776: E776.c myprob.c degfinder.c chiNPwrapper.c
#	gcc -c E776.c myprob.c  degfinder.c chiNPwrapper.c -std=c99 -ffast-math -O3 -fnested-functions -m32
#	gcc E776.o myprob.o  degfinder.o chiNPwrapper.o -o e776  $(LDFLAGS) $(local_LDFLAGS) -std=c99 -fnested-functions -ffast-math -O3 -m32
#	time ./e776
##	rm *.o
##	rm e776

e776-pc: E776.c myprob.c degfinder.c chiNPwrapper.c
	gcc -c E776.c myprob.c  degfinder.c chiNPwrapper.c ${CFLAGS}
	gcc E776.o myprob.o  degfinder.o chiNPwrapper.o -o e776  ${CFLAGS} $(INCFLAGS) $(LDFLAGS)
	./e776

dis: E776-dis.glb
	globes -s E776-dis.glb -o dis.dat

app: E776.glb
	globes -s E776.glb -o app.dat
