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


export prefix = ${HOME}/software/installed/globes-devel
export exec_prefix = ${prefix}
export libdir = ${exec_prefix}/lib
export globesconf= $(exec_prefix)/bin/globes-config

export solar_dir=solar/

export GCC = gcc
export GPP = g++

GLB_LDFLAGS = $(shell $(globesconf) --libs)
GLB_LDFLAGS_PROF = -g -pg -L/home/jkopp/software/installed/globes-devel/lib \
		   -Wl,-rpath -Wl,/home/jkopp/software/installed/globes-devel/lib \
		   -lglobes -L/home/jkopp/software/installed/lib -lgsl -lgslcblas \
#		   -nodefaultlibs -static-libgcc -lgcc -lstdc++ -lm_p -lc
#		   -nodefaultlibs -static -static-libgcc -lgcc -lstdc++ -lm -lc
INCFLAGS = $(shell $(globesconf) --include) -I${HOME}/projects/globes/Distribution/source
LDFLAGS = $(GLB_LDFLAGS) \
	  -L./atm/LibMath -L./atm/LibPhys -L./atm/LibWrap \
	  -L./$(solar_dir)/LibMath -L./$(solar_dir)/LibPhys \
	  -L./$(solar_dir)/LibWrap -L./$(solar_dir)/LibProb \
	  -latmwrapper -latmphysics -latmmathres \
      	  -lsolarwrapper -lsolarphysics -lsolarmathres -lprobsun -static
          #-lmontecubes

# For some reason, -D_BSD_SOURCE is required for mathematical constants
# such as M_PI, to be defined
export CFLAGS   = -std=c99 -Wall -g -O2 -D_BSD_SOURCE -fgnu89-inline -ffast-math -Wno-write-strings $(INCFLAGS)
export CXXFLAGS = -Wall -g -O2 -Wno-write-strings $(INCFLAGS)

BIN     = nu
OBJ     = nu.o const.o sys.o snu.o degfinder.o sensitivities.o prem.o \
	  mb.o e776.o c12.o icarus.o prior.o
HEADER  = nu.h const.h snu.h sbl/definitions.h reactors/definitions.h
SUBDIRS = sbl reactors

all: $(BIN)

$(BIN): $(OBJ) $(SUBDIRS)
#	${GPP} $(CFLAGS) $$(ls *.o sbl/*.o reactors/*.o) -o $(BIN) $(INCFLAGS) $(LDFLAGS)
	${GPP} $(CFLAGS) $(shell ls *.o sbl/*.o reactors/*.o) -o $(BIN) $(INCFLAGS) $(LDFLAGS)

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

