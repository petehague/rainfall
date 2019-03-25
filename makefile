CXX = g++
CPPFLAGS = -pg -g -std=c++11

COREOBJS = bin/main.o bin/number.o bin/queue.o bin/model.o bin/atom.o

all: rainfall basc hmc

rainfall: $(COREOBJS) bin/testfunc.o
	$(CXX) $(CPPFLAGS) $(COREOBJS) bin/testfunc.o -orainfall

basc: $(COREOBJS) bin/bascfunc.o
	$(CXX) $(CPPFLAGS) $(COREOBJS) bin/bascfunc.o -obasc

hmc: bin/pick.o bin/hmc.o
	$(CXX) $(CPPFLAGS) bin/pick.o bin/hmc.o -ohmctest

bin/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $< -c -o$@

clean:
	rm -rf bin/*.o
	rm -f hmctest
	rm -f rainfall
	rm -f basc
