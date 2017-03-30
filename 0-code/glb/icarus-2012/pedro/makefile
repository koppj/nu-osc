# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1



prefix = /usr/local
#prefixfma = /home/accioly/installedsoftware
exec_prefix = ${prefix}
#exec_prefixfma = ${prefixfma}
libdir = ${exec_prefix}/lib
#libdirfma = ${exec_prefixfma}/lib
globesconf= $(exec_prefix)/bin/globes-config
#globesconffma= $(exec_prefixfma)/bin/globes-config

#local_CFLAGS = -g -O4 
local_CFLAGS = -g -O4 -fnested-functions

INCFLAGS:=$(shell $(globesconf) --include)
#INCFLAGSfma:=$(shell $(globesconffma) --include)
#local_LDFLAGS:=$(shell $(globesconf) --libs)
local_LDFLAGS:=-L/usr/local/lib -lglobes -L/sw/lib -lgsl -lgslcblas -lm -lm
local_LTLDFLAGS:=$(shell $(globesconf) --ltlibs)
#local_LTLDFLAGSfma:=$(shell $(globesconffma) --ltlibs)

icarus: icarus.c snu.c degfinder.c chiNPwrapper.c
	gcc -c icarus.c    snu.c degfinder.c chiNPwrapper.c -std=c99 -ffast-math -O3 -fnested-functions #-fbounds-check  #-m32 
	gcc icarus.o    snu.o degfinder.o chiNPwrapper.o -o icarus  $(LDFLAGS) $(local_LDFLAGS) -std=c99 -fnested-functions -ffast-math -O3 #-fbounds-check #-m32
	./icarus


