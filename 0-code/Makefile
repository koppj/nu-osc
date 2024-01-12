# GLoBES -- General LOng Baseline Experiment Simulator
# (C) 2002 - 2007,  The GLoBES Team
#
# GLoBES is mainly intended for academic purposes. Proper
# credit must be given if you use GLoBES or parts of it. Please
# read the section 'Credit' in the README file.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# The configure script should have replaced anything as required to
# obtain a working Makefile, which may serve as a template for own
# applications.
#
# This simple Makefile is for the GLoBES examples
#
# Compile example files with ``make example1'' etc.
#
# This Makefile assumes dynamic libraries, installed to either
# the default prefix /usr/local/ or to a user-defined directory 
# called ${prefix}.
#
# For linking against a specific version of GLoBES, libglobes.so can be 
# replaced by the respective library, such as libglobes.so.0.0.1


# Note: - To compile with MPI support, add -lmpich option to LDFLAGS
#         and use mpicc instead of gcc
#       - To include profiling information, use -pg option in CFLAGS
#         and LDFLAGS and use -static option in LDFLAGS,
#         and link against -lc_p if available (use GLB_LDFLAGS_PROF instead
#         of GLB_LDFLAGS)
#       - don't use -ffast-math as this prevents consistency checks on osc. params in snu.c

# Profiling?
#export PROF_FLAGS=-pg -O0 -no-pie
export EXTRA_FLAGS=$(PROF_FLAGS)

# Determine operating system
UNAME_S := $(shell uname -s)

# Set platform-specific variables
export GCC = gcc
export GPP = g++
export prefix = ${HOME}/Dropbox/nu-rec/globes/installed/
ifeq ($(UNAME_S),Linux)
  INCFLAGS = $(shell $(globesconf) --include) \
             -I../../globes/Distribution/source
  export CFLAGS   = -std=c99 -Wall -g -O2 -D_DEFAULT_SOURCE -fPIC -fgnu89-inline \
                    -Wno-write-strings $(INCFLAGS) $(EXTRA_FLAGS)
  export CXXFLAGS = -std=c++11 -Wall -g -O2 -Wno-write-strings \
    		      $(INCFLAGS) $(EXTRA_FLAGS)
endif
ifeq ($(UNAME_S),Darwin)
  INCFLAGS = $(shell $(globesconf) --include) \
	     -I../../globes/Distribution/source
  export CFLAGS   = -std=c99 -Wall -g -O2 -D_BSD_SOURCE -fgnu89-inline \
                    -Wno-write-strings $(INCFLAGS) $(EXTRA_FLAGS)
  export CXXFLAGS = -std=c++11 -Wall -g -O2 -Wno-write-strings \
                    $(INCFLAGS) $(EXTRA_FLAGS)
endif

export exec_prefix = ${prefix}
export libdir = ${exec_prefix}/lib
export globesconf = $(exec_prefix)/bin/globes-config

ifeq ($(PROF_FLAGS),"")
  GLB_LDFLAGS = $(shell $(globesconf) --libs)
else
  GLB_LDFLAGS = -g -pg -L${prefix}/lib \
	        -Wl,-rpath -Wl,${prefix}/lib \
  	        -lglobes -L/home/jkopp/software/installed/lib -lgsl -lgslcblas \
	        -static-libgcc -lgcc -lstdc++ -lm -lc
#	         -nodefaultlibs -static-libgcc -lgcc -lstdc++ -lm_p -lc
#	         -nodefaultlibs -static -static-libgcc -lgcc -lstdc++ -lm -lc
endif

LDFLAGS = $(GLB_LDFLAGS)

BIN     = nu
OBJ     = nu.o const.o sys.o snu.o degfinder.o sensitivities.o prem.o \
	  mb.o mb-jk.o muboone.o e776.o c12.o icarus-2012.o icarus-2014.o \
          opera-2013.o minos.o prior.o  util.o
HEADER  = nu.h const.h snu.h LSND.h KARMEN.h constraints.h \
	  sbl/definitions.h reactors/definitions.h db-neos/db-neos.h osc_decay/osc_decay.h
SUBDIRS = sbl reactors db-neos

all: $(BIN)

$(BIN): $(OBJ) $(SUBDIRS)
	${GPP} $(CFLAGS) $(INCFLAGS) $(shell ls *.o sbl/*.o reactors/*.o db-neos/*.o) \
	 	 $(LDFLAGS) -o $(BIN)

%.o : %.c $(HEADER)
	${GCC} $(CFLAGS) $(INCFLAGS) -c $<

%.o : %.cc $(HEADER)
	${GPP} $(CXXFLAGS) $(INCFLAGS) -c $<

.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)
clean:
	rm -f $(BIN) $(OBJ)
	@for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

