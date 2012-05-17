# --------------------------------------------------------------------
# Makefile for mT2 code
# --------------------------------------------------------------------

include $(ROOTSYS)/test/Makefile.arch

# Configuration
# -------------

ROOTCONFIG = root-config
GENIECONFIG = genie-config

CFLAGS     = $(ROOTCFLAGS) -Wall -g -O2

INCFLAGS   = -I$(shell $(ROOTCONFIG) --incdir)
LIBFLAGS   = -lstdc++ -lMathMore -lMinuit -lGeom -lEG -lEGPythia6 -lxml2 -llog4cpp \
	     $(shell $(ROOTCONFIG) --glibs) $(shell $(GENIECONFIG) --libs)


# Files included in project
# -------------------------

BIN    = nc
OBJ    = nc.o
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


