# Pedro G. K. Bertella | 05/07/2021

CXX=g++
CXXFLAGS=-O2 -fPIC -Wall -MMD -MP -MF -nostartfiles 
LDFLAGS = -shared
LIBS=-lCGAL -lgmp -lgmpxx
MYEXEC=build/libApprox.so
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
MYOBJECTS=$(SOURCES:.cpp=.o)
.PHONY: all clean

$(shell mkdir -p build)

all: $(TARGET)

all: $(MYEXEC)
$(MYEXEC): $(MYOBJECTS)
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

clean:
	rm -f $(MYEXEC) $(MYOBJECTS)
