include $(ROOTSYS)/etc/Makefile.arch

# Configuration
# -------------

ROOTCONFIG = root-config
GENIECONFIG = genie-config

CFLAGS     = $(ROOTCFLAGS) -fopenmp -Wall -g -O2

INCFLAGS   = -I$(shell $(ROOTCONFIG) --incdir) -I$(MYSOFTWARE)/include/GENIE
LIBFLAGS   = -lstdc++ -lMathMore -lMinuit -lGeom -lEG -lEGPythia6 -lxml2 -llog4cpp \
	     -Wl,-rpath -Wl,${HOME}/local/software/lhapdf-5.9.1/lib/.libs -lLHAPDF \
	     $(shell $(ROOTCONFIG) --glibs) $(shell $(GENIECONFIG) --libs)


# Files included in project
# -------------------------

BIN    = minos
OBJ    = minos.o
HEADER =


# Rules
# -----

$(BIN): $(OBJ)
	@echo "Linking $@"
	@g++ $(CFLAGS) $(INCFLAGS) $(OBJ) $(LIBFLAGS) $(LDFLAGS) -o $(BIN)

%.o : %.c $(HEADER)
	@echo "Compling $@"
	@g++ $(CFLAGS) $(INCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(BIN) $(OBJ)


