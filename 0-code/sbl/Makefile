OBJ     = linear_algebra.o mynrfunc.o read.o miniboone.o minib-anti.o flux_Xsection.o \
	  nomad.o karmen.o lsnd.o cdhs.o atm-de-dm.o
HEADER  = def-reactors.h definitions.h minimize-simplex.h minimization.h mynrfunc.h
SUBDIRS = 

all: $(OBJ)

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

