CC=g++
CXX=g++
CXXLD=g++
CCLD=gcc

CXXFLAGS= -g -Wall -I -I$  -pthread -fPIC -DPIC
CCFLAGS=$(CXXFLAGS)

ROOTSYS = /opt64/root/new/

ROOTCFLAGS   := $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTGLIBS    := $(shell ${ROOTSYS}/bin/root-config --glibs)
ROOTLDFLAGS  := $(shell ${ROOTSYS}/bin/root-config --ldflags)
CXXFLAGS += $(ROOTCFLAGS)

INCLUDES:= analyzer.h main.h sort_functions.h structures.h unpacker.h

OBJECTS:= analyzer.o main.o sort_functions.o unpacker.o

LIBS  = -lm $(ROOTGLIBS) -lz

SRCS:= analyzer.cpp main.cpp sort_functions.cpp unpacker.cpp	

all: main

main:  $(SRCS) $(OBJECTS) $(INCLUDES) 
	$(CXX) -o DANCE_Analysis $(OBJECTS) $(CXXFLAGS) $(ROOTGLIBS) $(LIBS)

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -c $< 
%.o: %.cxx 
	$(CXX) $(CXXFLAGS) -c $< 
%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o DANCE_Analysis
# DO NOT DELETE

print_env:
	@echo $(ROOTSYS)
	@echo $(ROOTCFLAGS)
	@echo $(ROOTGLIBS)
	@echo $(ROOTLDFLAGS)
