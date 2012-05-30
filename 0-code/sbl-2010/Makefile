OBJ     = linear_algebra.o minimize-simplex.o mynrfunc.o read.o miniboone.o minib-anti.o flux_Xsection.o \
	  nomad.o karmen.o lsnd.o class_react-fit.o class_react-flux.o bugey.o palo-verde.o \
	  sbl_reactors.o chooz.o cdhs.o atm-de-dm.o
HEADER  = def-reactors.h definitions.h minimize-simplex.h minimization.h mynrfunc.h
SUBDIRS = 

all: $(OBJ)

%.o : %.c $(HEADER)
	gcc $(CFLAGS) $(INCFLAGS) -c $<

%.o : %.cc $(HEADER)
	g++ $(CXXFLAGS) $(INCFLAGS) -c $<

.PHONY: subdirs $(SUBDIRS)
subdirs: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)
clean:
	rm -f $(BIN) $(OBJ)

