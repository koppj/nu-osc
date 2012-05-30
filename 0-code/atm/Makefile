OBJ     = 
HEADER  =
SUBDIRS = 

all: $(OBJ)

%.o : %.c $(HEADER)
	gcc $(CFLAGS) $(INCFLAGS) -c $<

%.o : %.cc $(HEADER)
	g++ $(CXXFLAGS) $(INCFLAGS) -c $<

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)

clean:
	rm -f $(BIN) $(OBJ)

