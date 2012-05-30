OBJ     = bugey.o gallium.o class_fit.o dayabay.o kamland.o main.o \
   	  qrom.o reno.o chooz.o class_flux.o dc.o linear_algebra.o \
	  palo-verde.o read.o sbl_reactors.o
HEADER  = definitions.h qrom.h
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

