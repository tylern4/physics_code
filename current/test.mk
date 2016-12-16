UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -std=c++11
endif

ROOTLIBS	= $(shell root-config --libs)
CXXFLAGS =      -O2 -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
LDFLAGS = $(shell root-config --libs)
CXX = g++
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

TARGET  = e1d

.PHONY: all
all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CXX) $(CXXFLAGS) $< -c -o $@