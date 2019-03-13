CXX = g++
CPPFLAGS = -g -std=c++11

rainfall: bin/number.o bin/main.o bin/testfunc.o
	$(CXX) $(CPPFLAGS) bin/main.o bin/number.o bin/testfunc.o -orainfall

basc: bin/number.o bin/main.o bin/bascfunc.o
	$(CXX) $(CPPFLAGS) bin/main.o bin/number.o bin/bascfunc.o -obasc

hmc: bin/pick.o bin/hmc.o
	$(CXX) $(CPPFLAGS) bin/pick.o bin/hmc.o -ohmctest

bin/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $< -c -o$@

clean:
	rm -rf bin/*.o
	rm -f hmctest
	rm -f rainfall
	rm -f basc
