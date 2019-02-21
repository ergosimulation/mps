UNAME_S := $(shell uname -s)
export CPPFLAGS+= -O3 -std=c++11
ifeq ($(UNAME_S),Linux)
	export CPPFLAGS += -static -Wl,--no-as-needed
endif
ifeq ($(UNAME_S),MINGW64_NT-10.0)
	export CPPFLAGS += -static -Wl,--no-as-needed
endif

# NAME OF LIBRARY
MPSLIB = mpslib/mpslib.a

# LINK LIBRARIES
LDLIBS =  -lstdc++ -lpthread

all: mps_genesim mps_snesim_list mps_snesim_tree

.PHONY: mpslib
mpslib:
	make -C mpslib

mps_genesim: mpslib
	$(CXX) $(CPPFLAGS) mps_genesim.cpp ENESIM_GENERAL.cpp $(MPSLIB) -o $@ -I mpslib/ $(LDLIBS)

mps_snesim_tree: mpslib
	$(CXX) $(CPPFLAGS) mps_snesim_tree.cpp SNESIMTree.cpp $(MPSLIB) -o $@ -I mpslib/ $(LDLIBS)

mps_snesim_list: mpslib
	$(CXX) $(CPPFLAGS) mps_snesim_list.cpp SNESIMList.cpp $(MPSLIB) -o $@ -I mpslib/ $(LDLIBS)


.PHONY: clean
clean:
	rm -f *.o mps mpslib/*.o $(MPSLIB)

cleanexe:
	rm -f *.o mps *.exe mps_genesim mps_snesim_tree mps_snesim_list mpslib/*.o $(MPSLIB)

.PHONY: cleano
cleano:
	rm -f *.o mpslib/*.o $(MPSLIB)
