##*************************##
##  Christopher J. Prokop  ##
##  cprokop@lanl.gov       ##
##  Makefile               ## 
##  Last Edit: 04/15/19    ##  
##*************************##

CC=g++
CXX=g++
CXXLD=g++
CCLD=gcc

CXXFLAGS= -g -Wall -I -I$  -pthread -fPIC -DPIC
CCFLAGS=$(CXXFLAGS)

ROOTCFLAGS   := $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOTGLIBS    := $(shell ${ROOTSYS}/bin/root-config --glibs)
ROOTLDFLAGS  := $(shell ${ROOTSYS}/bin/root-config --ldflags)
CXXFLAGS += $(ROOTCFLAGS)

INCLUDES:= message.h calibrator.h validator.h eventbuilder.h analyzer.h main.h sort_functions.h unpacker.h unpack_vx725_vx730.h structures.h global.h 

OBJECTS:= message.o calibrator.o validator.o eventbuilder.o analyzer.o main.o sort_functions.o unpacker.o unpack_vx725_vx730.o

LIBS  = -lm $(ROOTGLIBS) -lz -lbz2

SRCS:= message.cpp calibrator.cpp validator.cpp eventbuilder.cpp analyzer.cpp main.cpp sort_functions.cpp unpacker.cpp unpack_vx725_vx730.cpp 

all: main

main:  $(SRCS) $(OBJECTS) $(INCLUDES) 
	$(CXX) -o DANCE_Analysis $(OBJECTS) $(CXXFLAGS) $(ROOTGLIBS) $(LIBS)

%.o: %.cpp ${INCLUDES} 
	$(CXX) $(CXXFLAGS) -c $< 
%.o: %.cxx ${INCLUDES}
	$(CXX) $(CXXFLAGS) -c $< 
%.o: %.c ${INCLUDES}
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -f *.o DANCE_Analysis
# DO NOT DELETE

print_env:
	@echo $(ROOTSYS)
	@echo $(ROOTCFLAGS)
	@echo $(ROOTGLIBS)
	@echo $(ROOTLDFLAGS)
