OBJ     = bugey.o danss.o gallium.o class_fit.o dayabay.o db-flux.o \
          kamland.o dc.o main.o qrom.o reno.o chooz.o class_flux.o class_minimize.o \
	  linear_algebra.o palo-verde.o read.o sbl_reactors.o class_read_2dim.o
HEADER  = definitions.h qrom.h
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

