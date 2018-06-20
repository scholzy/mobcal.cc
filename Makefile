appname := mobcal

# CXX := g++-8
# CXXFLAGS := -std=c++11 -g -O2 -march=native -Wall -pedantic -ffast-math -fopenmp

# CXX := clang++
# CXXFLAGS := -std=c++11 -g -Wall -pedantic -O2 -ffast-math

CXX := icc
CXXFLAGS := -std=c++11 -g -Wall -pedantic -Ofast -qopenmp

srcfiles := $(shell find . -name "*.cc")
objects := $(patsubst %.cc, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	-rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

test: $(appname)
	@-./$(appname) tests/n2o.inp

clean:
	-rm -f $(objects)

dist-clean: clean
	-rm -f *~ .depend
